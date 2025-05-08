#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <vector>
#include <limits>
#include "unordered_dense.hpp"

namespace py = pybind11;

template <typename T>
py::array to_numpy(
	T* output,
	const uint64_t sx, const uint64_t sy
) {
	py::capsule capsule(output, [](void* ptr) {
		if (ptr) {
			delete[] static_cast<T*>(ptr);
		}
	});

	uint64_t width = sizeof(T);

	// C order
	return py::array_t<T>(
		{sx,sy},
		{sy * width, width},
		output,
		capsule
	);
}

py::list compute_components(
	const py::array_t<uint32_t> &edges_arr,
	const uint64_t Nv
) {
	const uint64_t Ne = edges_arr.size() >> 1;

    py::buffer_info buf = edges_arr.request();

    // Sanity check: 2D, C-contiguous, correct item size
    if (buf.ndim != 2 || buf.strides[1] != sizeof(uint32_t)) {
        throw std::runtime_error("Array must be 2D and C-contiguous");
    }

    uint32_t* edges = static_cast<uint32_t*>(buf.ptr);

	std::vector<bool> visited(Nv);
	std::vector< ankerl::unordered_dense::set<uint32_t> > index(Nv);

	for (size_t i = 0; i < 2 * Ne; i += 2) {
		uint32_t e1 = edges[i];
		uint32_t e2 = edges[i+1];

		index[e1].insert(e2);
		index[e2].insert(e1);
	}

	auto extract_component = [&](uint32_t start){
		std::vector<std::pair<uint32_t,uint32_t>> edge_list;
		std::vector<uint32_t> stack = { start };
		std::vector<uint32_t> parents = { std::numeric_limits<uint32_t>::max() };

		// printf("start: %d\n", start);

		while (stack.size()) {
			uint32_t node = stack.back();
			uint32_t parent = parents.back();
 
			stack.pop_back();
			parents.pop_back();

			if (node == parent) {
				continue;
			}
			else if (node < parent) {
				edge_list.emplace_back(node, parent);	
			}
			else {
				// printf("gt\n");
				edge_list.emplace_back(parent, node);
			}

			// check visited after because you can visit a node 
			// multiple times for different parents, but you don't
			// want to keep re-adding it to the stack.
			if (visited[node]) {
				// printf("visited.\n");
				continue;
			}

			visited[node] = true;

			for (uint32_t child : index[node]) {
				// printf("pb %d\n", child);
				stack.push_back(child);
				parents.push_back(node);
			}
		}

		// printf("done\n");

		if (edge_list.size() <= 1) {
			return py::array(py::dtype("uint32"), {0, 2});
		}

    	const size_t n = edge_list.size() - 1;

		uint32_t* flat_data = new uint32_t[n * 2]();
		for (uint64_t i = 1; i < edge_list.size(); i++) {
			flat_data[(i - 1) << 1] = edge_list[i].first;
			flat_data[((i - 1) << 1) + 1] = edge_list[i].second;
		}

		return to_numpy<uint32_t>(flat_data, n, 2);

		// uint64_t* flat_data_pairs = reinterpret_cast<uint64_t*>(flat_data);
		// std::sort(flat_data_pairs, flat_data_pairs + n);

		// uint64_t uniq = (n > 0) ? 1 : 0;
		// for (uint64_t i = 1; i < n; i++) {
		// 	if (flat_data_pairs[i] != flat_data_pairs[i-1]) {
		// 		uniq++;
		// 	}
		// }

		// uint64_t* uniq_data = new uint64_t[uniq];
		// uniq_data[0] = flat_data_pairs[0];
		// uint64_t j = 1;
		// for (uint64_t i = 1; i < n; i++) {
		// 	if (flat_data_pairs[i] != flat_data_pairs[i-1]) {
		// 		uniq_data[j] = flat_data_pairs[i];
		// 		j++;
		// 	}
		// }

		// delete[] flat_data_pairs;

		// return to_numpy<uint32_t>(reinterpret_cast<uint32_t*>(uniq_data), uniq, 2);
	};

	py::list forest;
	for (size_t i = 0; i < Ne * 2; i += 2) {
		uint32_t edge = edges[i];

		if (visited[edge]) {
			continue;
		}

		forest.append(
			extract_component(edge)
		);
	}

	return forest;
}

PYBIND11_MODULE(fastosteoid, m) {
	m.doc() = "Accelerated osteoid functions."; 
	m.def("compute_components", &compute_components, "Find skeleton components.");
}
