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

// uint64_t c_eytzinger_binary_search(uint64_t x, uint64_t* array, size_t N) {
//     // int64_t block_size = 8; // two cache lines 64 * 2 / 8
//     uint64_t k = 1;
//     while (k <= (uint64_t)N) {
//         // __builtin_prefetch(array + k * block_size);
//         // multiply by 2 b/c index is [label, pos, label, pos]
//         k = 2 * k + (array[(k - 1) << 1] < x); 
//     }
//     k >>= mb_ffs(~k);
//     k -= 1;

//     if (k >= 0 && array[k << 1] == x) {
//         return k;
//     }

//     return -1;
// }

// uint32_t c_eytzinger_sort_indices(
//     uint32_t* array, uint32_t n, 
//     uint32_t i /*=0*/, uint32_t k /*=1*/
// ) {
//     if (k <= n) {
//         i = c_eytzinger_sort_indices(array, n, i, 2 * k);
//         array[k-1] = i;
//         i++;
//         i = c_eytzinger_sort_indices(array, n, i, 2 * k + 1);
//     }
//     return i;
// }

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

template <typename T>
py::array_t<T> pairs_to_numpy(
	const std::vector<std::pair<T, T>>& pairs
) {
		size_t n = pairs.size();

		// Create an uninitialized NumPy array of shape (n, 2)
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
	m.def("has_cycle", &has_cycle, "Check whether this connected component has a cycle.");
	m.def("compute_components", &compute_components, "Find skeleton components.");
}
