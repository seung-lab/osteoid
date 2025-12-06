#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <limits>
#include <span>
#include <vector>

#include "chunk.hpp"
#include "unordered_dense.hpp"

namespace py = pybind11;

template <typename T>
py::array_t<T> pairs_to_numpy(
	const std::span<std::pair<T, T>>& pairs
) {
	size_t n = pairs.size();
	py::array_t<T> arr({n, size_t(2)});
	auto buf = arr.template mutable_unchecked<2>();

	for (size_t i = 0; i < n; ++i) {
			buf(i, 0) = pairs[i].first;
			buf(i, 1) = pairs[i].second;
	}

	return arr;
}

template <typename T>
std::vector<std::pair<T,T>> unique(
	std::vector<std::pair<T,T>>& input
) {
	if (input.size() == 0) {
		return std::vector<std::pair<T,T>>();
	}

	std::sort(input.begin(), input.end());
	
	const uint64_t N = input.size();

	std::vector<std::pair<T,T>> uniq;
	uniq.reserve(N / 10);

	uniq.push_back(input[0]);
	
	for (uint64_t i = 1; i < N; i++) {
		if (input[i] != input[i-1]) {
			uniq.push_back(input[i]);
		}
	}

	return uniq;
}

template <typename EDGE_T>
py::list compute_components_impl(
	const py::array_t<EDGE_T> &edges_arr,
	const uint64_t Nv
) {
	const uint64_t Ne = edges_arr.size() >> 1;

	py::buffer_info buf = edges_arr.request();

	if (buf.ndim != 2 || buf.strides[1] != sizeof(EDGE_T)) {
		throw std::runtime_error("Array must be 2D and C-contiguous");
	}

	EDGE_T* edges = static_cast<EDGE_T*>(buf.ptr);

	std::vector<bool> visited(Nv);
	std::vector< ankerl::unordered_dense::set<EDGE_T> > index(Nv);

	for (size_t i = 0; i < 2 * Ne; i += 2) {
		EDGE_T e1 = edges[i];
		EDGE_T e2 = edges[i+1];

		index[e1].insert(e2);
		index[e2].insert(e1);
	}

	auto extract_component = [&](EDGE_T start) -> py::array_t<EDGE_T> {
		std::vector<std::pair<EDGE_T, EDGE_T>> edge_list;
		std::vector<EDGE_T> stack = {};
		std::vector<EDGE_T> parents = {};

		visited[start] = true;
		for (EDGE_T child : index[start]) {
			stack.push_back(child);
			parents.push_back(start);
		}

		while (stack.size()) {
			EDGE_T node = stack.back();
			EDGE_T parent = parents.back();
 
			stack.pop_back();
			parents.pop_back();

			if (node == parent) {
				continue;
			}
			else if (node < parent) {
				edge_list.emplace_back(node, parent);
			}
			else {
				edge_list.emplace_back(parent, node);
			}

			// check visited after because you can visit a node 
			// multiple times for different parents, but you don't
			// want to keep re-adding it to the stack.
			if (visited[node]) {
				continue;
			}

			visited[node] = true;

			for (EDGE_T child : index[node]) {
				if (child == parent) {
					continue;
				}
				stack.push_back(child);
				parents.push_back(node);
			}
		}

		if (edge_list.size() == 0) {
			return py::array(py::dtype("uint64"), {0, 2});
		}

		auto uniq = unique<EDGE_T>(edge_list);
		return pairs_to_numpy<EDGE_T>(uniq);
	};

	py::list forest;
	for (size_t i = 0; i < Ne * 2; i += 2) {
		EDGE_T edge = edges[i];

		if (visited[edge]) {
			continue;
		}

		forest.append(
			extract_component(edge)
		);
	}

	return forest;
}	

py::list compute_components(
	const py::array &edges_arr,
	const uint64_t Nv
) {
	py::buffer_info buf = edges_arr.request();

	if (buf.ndim != 2) {
		throw std::runtime_error("Array must be 2D and C-contiguous");
	}

	int data_width = edges_arr.dtype().itemsize();

	if (data_width == 1) {
		if (buf.strides[1] != sizeof(uint8_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return compute_components_impl<uint8_t>(edges_arr, Nv);
	}
	else if (data_width == 2) {
		if (buf.strides[1] != sizeof(uint16_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return compute_components_impl<uint16_t>(edges_arr, Nv);
	}
	else if (data_width == 4) {
		if (buf.strides[1] != sizeof(uint32_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return compute_components_impl<uint32_t>(edges_arr, Nv);
	}
	else {
		if (buf.strides[1] != sizeof(uint64_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return compute_components_impl<uint64_t>(edges_arr, Nv);
	}
}

template <typename VERT_T, typename EDGE_T>
py::dict chunk_skeleton_impl(
	const py::array_t<VERT_T>& vertex_arr,
	const py::array_t<EDGE_T>& edges_arr,
	const float cx, const float cy, const float cz,
	const std::optional<float> origin_x = std::nullopt,
	const std::optional<float> origin_y = std::nullopt,
	const std::optional<float> origin_z = std::nullopt
) {
	py::buffer_info vbuf = vertex_arr.request();
	if (vbuf.ndim != 2 || vbuf.strides[1] != sizeof(VERT_T) || vbuf.shape[1] != 3) {
		throw std::runtime_error("Array must be 3D and C-contiguous");
	}
	VERT_T* vertices = static_cast<VERT_T*>(vbuf.ptr);

	py::buffer_info ebuf = edges_arr.request();
	if (ebuf.ndim != 2 || ebuf.strides[1] != sizeof(EDGE_T) || ebuf.shape[1] != 2) {
		throw std::runtime_error("Array must be 2D and C-contiguous");
	}
	EDGE_T* edges = static_cast<EDGE_T*>(ebuf.ptr);

	const uint64_t Nv = vbuf.shape[0];
	const uint64_t Ne = ebuf.shape[0];

	auto line_grid = fastosteoid::chunk::chunk_line(
		vertices, Nv,
		edges, Ne,
		cx, cy, cz,
		origin_x, origin_y, origin_z
	);

	VERT_T minx = vertices[0], miny = vertices[1], minz = vertices[2];
	VERT_T maxx = vertices[0], maxy = vertices[1], maxz = vertices[2];
	
	for (uint64_t i = 0; i < 3 * Nv; i += 3) {
		minx = std::min(minx, vertices[i + 0]);
		miny = std::min(miny, vertices[i + 1]);
		minz = std::min(minz, vertices[i + 2]);
		maxx = std::max(maxx, vertices[i + 0]);
		maxy = std::max(maxy, vertices[i + 1]);
		maxz = std::max(maxz, vertices[i + 2]);
	}

	if (origin_x.has_value()) {
		minx = *origin_x;
	}
	if (origin_y.has_value()) {
		miny = *origin_y;
	}
	if (origin_z.has_value()) {
		minz = *origin_z;
	}

	int64_t gsx = std::max(static_cast<int64_t>(1), static_cast<int64_t>(std::ceil((maxx - minx) / cx)));
	int64_t gsy = std::max(static_cast<int64_t>(1), static_cast<int64_t>(std::ceil((maxy - miny) / cy)));
	int64_t gsz = std::max(static_cast<int64_t>(1), static_cast<int64_t>(std::ceil((maxz - minz) / cz)));

	py::dict chunked_skeletons;
	int64_t idx = 0;
	for (int64_t gz = 0; gz < gsz; gz++) {
		for (int64_t gy = 0; gy < gsy; gy++) {
			for (int64_t gx = 0; gx < gsx; gx++, idx++) {
				const auto& line_obj = line_grid[idx];

				if (line_obj.points.empty() || line_obj.edges.empty()) {
					continue;
				}

				size_t num_verts = line_obj.points.size() / 3;
				size_t num_edges = line_obj.edges.size() / 2;
				
				auto chunk_vertices = py::array_t<VERT_T>({
					static_cast<py::ssize_t>(num_verts), 
					static_cast<py::ssize_t>(3)
				});
				auto chunk_edges = py::array_t<EDGE_T>({
					static_cast<py::ssize_t>(num_edges), 
					static_cast<py::ssize_t>(2)
				});

				auto verts_buf = chunk_vertices.request();
				auto edges_buf = chunk_edges.request();
				
				VERT_T* verts_ptr = static_cast<VERT_T*>(verts_buf.ptr);
				EDGE_T* edges_ptr = static_cast<EDGE_T*>(edges_buf.ptr);
				
				std::memcpy(verts_ptr, line_obj.points.data(), line_obj.points.size() * sizeof(VERT_T));
				std::memcpy(edges_ptr, line_obj.edges.data(), line_obj.edges.size() * sizeof(EDGE_T));
				
				py::tuple key = py::make_tuple(gx, gy, gz);
				chunked_skeletons[key] = py::make_tuple(chunk_vertices, chunk_edges);
			}
		}
	}
   
	return chunked_skeletons;
}

py::dict chunk_skeleton(
	const py::array_t<float>& vertex_arr,
	const py::array& edges_arr,
	const float cx, const float cy, const float cz,
	const std::optional<float> origin_x = std::nullopt,
	const std::optional<float> origin_y = std::nullopt,
	const std::optional<float> origin_z = std::nullopt
) {
	py::buffer_info vbuf = vertex_arr.request();
	if (vbuf.ndim != 2) {
		throw std::runtime_error("Vertex array must be 2D");
	}
	if (vbuf.strides[1] != sizeof(float)) {
		throw std::runtime_error("Vertex array must be C-contiguous");
	}
	
	py::buffer_info ebuf = edges_arr.request();
	if (ebuf.ndim != 2) {
		throw std::runtime_error("Edge array must be 2D");
	}
	
	int edge_width = edges_arr.dtype().itemsize();
	
	if (edge_width == 1) {
		if (ebuf.strides[1] != sizeof(uint8_t)) {
			throw std::runtime_error("Edge array must be C-contiguous");
		}
		return chunk_skeleton_impl<float, uint8_t>(vertex_arr, edges_arr, cx, cy, cz, origin_x, origin_y, origin_z);
	}
	else if (edge_width == 2) {
		if (ebuf.strides[1] != sizeof(uint16_t)) {
			throw std::runtime_error("Edge array must be C-contiguous");
		}
		return chunk_skeleton_impl<float, uint16_t>(vertex_arr, edges_arr, cx, cy, cz, origin_x, origin_y, origin_z);
	}
	else if (edge_width == 4) {
		if (ebuf.strides[1] != sizeof(uint32_t)) {
			throw std::runtime_error("Edge array must be C-contiguous");
		}
		return chunk_skeleton_impl<float, uint32_t>(vertex_arr, edges_arr, cx, cy, cz, origin_x, origin_y, origin_z);
	}
	else { // 8 bytes
		if (ebuf.strides[1] != sizeof(uint64_t)) {
			throw std::runtime_error("Edge array must be C-contiguous");
		}
		return chunk_skeleton_impl<float, uint64_t>(vertex_arr, edges_arr, cx, cy, cz, origin_x, origin_y, origin_z);
	}
}

PYBIND11_MODULE(fastosteoid, m) {
	m.doc() = "Accelerated osteoid functions."; 
	m.def("compute_components", &compute_components, "Find skeleton components.");
	m.def("chunk_skeleton", &chunk_skeleton, "Cut a skeleton into a grid of chunks.");
}
