#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <vector>
#include <limits>

// #include "builtins.hpp"
#include "unordered_dense.hpp"

namespace py = pybind11;

template <typename EDGE_T>
py::list paths_to_pylist(const std::vector<std::vector<EDGE_T>>& paths) {
    py::list out(paths.size());

    for (size_t i = 0; i < paths.size(); i++) {
        const auto& p = paths[i];
        py::array_t<EDGE_T> arr(p.size());
        auto* dst = static_cast<EDGE_T*>(arr.mutable_data());
        std::memcpy(dst, p.data(), p.size() * sizeof(EDGE_T));
        out[i] = arr;
    }

    return out;
}

template <typename T>
py::array_t<T> pairs_to_numpy(
	const std::vector<std::pair<T, T>>& pairs
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
py::tuple linked_paths_impl(const py::array_t<EDGE_T>& edges_arr) {
	auto buf = edges_arr.request();
	if (buf.ndim != 2 || buf.shape[1] != 2) {
		throw std::runtime_error("edges must be an (N,2) array");
	}
	const uint64_t Ne = edges_arr.size() >> 1;
	EDGE_T* edges = static_cast<EDGE_T*>(buf.ptr);


	uint64_t Nv = 0;
	for (size_t i = 0; i < 2 * Ne; i++) {
		Nv = std::max(Nv, static_cast<uint64_t>(edges[i]));
	}
	Nv++;

	std::vector<bool> visited(Nv);
	std::vector< ankerl::unordered_dense::set<EDGE_T> > index(Nv);

	for (size_t i = 0; i < 2 * Ne; i += 2) {
		EDGE_T e1 = edges[i];
		EDGE_T e2 = edges[i+1];

		index[e1].insert(e2);
		index[e2].insert(e1);
	}

	auto extract_branches = [&](EDGE_T start) -> py::tuple {
		std::vector<EDGE_T> stack;
		std::vector<EDGE_T> parents;

		std::vector<std::pair<EDGE_T, EDGE_T>> edge_list;
		std::vector<std::vector<EDGE_T>> paths;
 		std::vector<EDGE_T> path;

 		path.push_back(start);

		visited[start] = true;
		for (EDGE_T child : index[start]) {
			stack.push_back(child);
			parents.push_back(start);
		}

		bool cycle_detected = false;

		while (stack.size()) {
			EDGE_T node = stack.back();
			EDGE_T parent = parents.back();

			stack.pop_back();
			parents.pop_back();

			if (node == parent) {
				continue;
			}

 			path.push_back(node);

			// check visited after because you can visit a node 
			// multiple times for different parents, but you don't
			// want to keep re-adding it to the stack.
			if (visited[node]) {
				continue;
			}

			while (true) {
				visited[node] = true;
				
				bool picked = false;
				EDGE_T next = 0;
				
				for (EDGE_T child : index[node]) {
					if (child == parent) {
						continue;
					}
					else if (visited[child]) {
						cycle_detected = true;
						continue;
					}
					else if (picked) {
						if (node < child) {
							edge_list.emplace_back(node, child);
						}
						else {
							edge_list.emplace_back(child, node);
						}
						stack.push_back(child);
						parents.push_back(node);
					}
					else {
						picked = true;
						next = child;
						path.push_back(child);
					}
				}

				if (picked) {
					parent = node;
					node = next;
				}
				else {
					break;
				}
			}

			paths.push_back(path);
			path.clear();
		}

		auto uniq = unique<EDGE_T>(edge_list);
		py::tuple out(3);
		out[0] = paths_to_pylist<EDGE_T>(paths);
		out[1] = pairs_to_numpy<EDGE_T>(uniq);
		out[2] = cycle_detected;

		return out;
	};

	py::tuple ret(4);

	py::list all_paths;
	py::list all_edges;

	uint64_t N = 0;
	bool has_cycle = false;

	for (size_t i = 0; i < Ne * 2; i += 2) {
		EDGE_T edge = edges[i];

		if (visited[edge]) {
			continue;
		}

		py::tuple tup = extract_branches(edge);
		all_paths.append(tup[0]);
		all_edges.append(tup[1]);
		has_cycle |= tup[2].cast<bool>();;
		N++;
	}

	ret[0] = all_paths;
	ret[1] = all_edges;
	ret[2] = has_cycle;
	ret[3] = N;

	return ret;
}

py::tuple linked_paths(const py::array &edges_arr) {
	py::buffer_info buf = edges_arr.request();

	if (buf.ndim != 2 || buf.shape[1] != 2) {
		throw std::runtime_error("edges must be an (N,2) array");
	}

	int data_width = edges_arr.dtype().itemsize();

	if (data_width == 1) {
		if (buf.strides[1] != sizeof(uint8_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return linked_paths_impl<uint8_t>(edges_arr);
	}
	else if (data_width == 2) {
		if (buf.strides[1] != sizeof(uint16_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return linked_paths_impl<uint16_t>(edges_arr);
	}
	else if (data_width == 4) {
		if (buf.strides[1] != sizeof(uint32_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return linked_paths_impl<uint32_t>(edges_arr);
	}
	else {
		if (buf.strides[1] != sizeof(uint64_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return linked_paths_impl<uint64_t>(edges_arr);
	}
}

// Ne = size of edges / 2
// Nv = number of vertices (max of edge values)
template <typename EDGE_T>
bool has_cycle_impl(const py::array_t<EDGE_T>& edges) {
	auto buf = edges.request();
	if (buf.ndim != 2 || buf.shape[1] != 2) {
		throw std::runtime_error("edges must be an (N,2) array");
	}
	constexpr EDGE_T sentinel = std::numeric_limits<EDGE_T>::max();
	const size_t Ne = buf.shape[0];

	if (Ne == 0) {
		return false;
	}

	const EDGE_T* ptr = static_cast<const EDGE_T*>(buf.ptr);

	EDGE_T maxv = 0;
	for (size_t i = 0; i < 2 * Ne; i++) {
		maxv = std::max(maxv, ptr[i]);
	}
	
	if (maxv == sentinel) {
		throw std::runtime_error("has_cycle: edges contains its maximum dtype value, which is used as a sentinel.");
	}

	const size_t Nv = static_cast<size_t>(maxv) + 1;

	std::vector< ankerl::unordered_dense::set<EDGE_T> > index(Nv);

	for (size_t i = 0; i < 2 * Ne; i += 2) {
		const EDGE_T e1 = ptr[i];
		const EDGE_T e2 = ptr[i+1];

		if (e1 == e2) {
			continue;
		}

		index[e1].insert(e2);
		index[e2].insert(e1);
	}

	const EDGE_T root = ptr[0];
	EDGE_T node = sentinel;
	EDGE_T parent = sentinel;

	std::vector<EDGE_T> stack;
	std::vector<EDGE_T> parents;

	parents.push_back(sentinel);
	stack.push_back(root);
	
	std::vector<bool> visited(Nv, false);

	while (!stack.empty()) {
		node = stack.back();
		parent = parents.back();

		stack.pop_back();
		parents.pop_back();

		visited[node] = true;

		for (EDGE_T child : index[node]) {
			if (child == parent) {
				continue;
			}
			else if (visited[child]) {
				return true;
			}

			stack.push_back(child);
			parents.push_back(node);
		}
	}

	return false;
}

bool has_cycle(const py::array& edges) {
	py::buffer_info buf = edges.request();

	if (buf.ndim != 2) {
		throw std::runtime_error("Array must be 2D and C-contiguous");
	}

	int data_width = edges.dtype().itemsize();

	if (data_width == 1) {
		if (buf.strides[1] != sizeof(uint8_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return has_cycle_impl<uint8_t>(edges);
	}
	else if (data_width == 2) {
		if (buf.strides[1] != sizeof(uint16_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return has_cycle_impl<uint16_t>(edges);
	}
	else if (data_width == 4) {
		if (buf.strides[1] != sizeof(uint32_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return has_cycle_impl<uint32_t>(edges);
	}
	else {
		if (buf.strides[1] != sizeof(uint64_t)) {
			throw std::runtime_error("Array must be C-contiguous");
		}
		return has_cycle_impl<uint64_t>(edges);
	}
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

PYBIND11_MODULE(fastosteoid, m) {
	m.doc() = "Accelerated osteoid functions."; 
	m.def("linked_paths", &linked_paths, "Compute the linked paths edge representation.");
	m.def("has_cycle", &has_cycle, "Check whether this connected component has a cycle.");
	m.def("compute_components", &compute_components, "Find skeleton components.");
}
