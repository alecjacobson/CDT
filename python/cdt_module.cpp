#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "cdt.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// Indices come in and go out as int64: unsigned indices are a foot-gun on the
// numpy side, where `T - 1` on a uint32 array silently wraps.
using IndexArray = nb::ndarray<const int64_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using PointArray = nb::ndarray<const double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

// Hand ownership of `data` to numpy without copying it. `cols == 0` asks for a
// one dimensional array of `rows` entries.
template <typename T>
nb::ndarray<nb::numpy, T> to_numpy(std::vector<T>&& data, size_t rows, size_t cols) {
    auto* held = new std::vector<T>(std::move(data));
    nb::capsule owner(held, [](void* p) noexcept { delete (std::vector<T>*)p; });
    const size_t shape[2] = {rows, cols};
    return nb::ndarray<nb::numpy, T>(held->data(), cols == 0 ? 1 : 2, shape, owner);
}

std::vector<uint32_t> to_indices(const IndexArray& A, int64_t expected_cols, const char* name) {
    if ((int64_t)A.shape(1) != expected_cols)
        throw std::invalid_argument(std::string(name) + " must be #" + name + " by " +
                                    std::to_string(expected_cols));
    const size_t n = A.shape(0) * A.shape(1);
    const int64_t* p = A.data();
    std::vector<uint32_t> out(n);
    for (size_t i = 0; i < n; i++) {
        if (p[i] < 0 || p[i] > (int64_t)UINT32_MAX)
            throw std::invalid_argument(std::string(name) + " contains an out of range index");
        out[i] = (uint32_t)p[i];
    }
    return out;
}

auto tetrahedralize(
    const PointArray& V,
    const IndexArray& F,
    const std::optional<IndexArray>& E,
    bool bounding_box,
    bool verbose,
    bool inner_only) {
    if (V.shape(1) != 3) throw std::invalid_argument("V must be #V by 3");

    const std::vector<uint32_t> f = to_indices(F, 3, "F");
    const std::vector<uint32_t> e = E ? to_indices(*E, 2, "E") : std::vector<uint32_t>();

    cdt::Options opts;
    opts.bounding_box = bounding_box;
    opts.verbose = verbose;
    opts.inner_only = inner_only;

    cdt::Result r = [&] {
        // Long, and it prints; let other threads run and let ^C be noticed.
        nb::gil_scoped_release release;
        return cdt::tetrahedralize(
            V.data(), (uint32_t)V.shape(0),
            f.data(), (uint32_t)F.shape(0),
            e.empty() ? nullptr : e.data(), (uint32_t)(E ? (*E).shape(0) : 0),
            opts);
    }();

    const size_t nv = r.num_vertices(), nt = r.num_tets();
    std::vector<int64_t> tets(r.tets.begin(), r.tets.end());
    std::vector<int64_t> vmap(r.vertex_map.begin(), r.vertex_map.end());
    const size_t n_in = vmap.size();

    return std::make_tuple(
        to_numpy(std::move(r.vertices), nv, 3),
        to_numpy(std::move(tets), nt, 4),
        to_numpy(std::move(r.labels), nt, 0),
        to_numpy(std::move(vmap), n_in, 0),
        (int64_t)r.num_steiner_vertices,
        r.is_polyhedron,
        r.face_recovery_ok);
}

} // namespace

NB_MODULE(_cdt, m) {
    m.doc() = "Constrained Delaunay tetrahedrization (CDT)";

    m.def("tetrahedralize", &tetrahedralize,
        "V"_a, "F"_a, "E"_a = nb::none(),
        "bounding_box"_a = false, "verbose"_a = false, "inner_only"_a = false,
        R"(Constrained Delaunay tetrahedrization of a PLC.

Returns (vertices, tets, labels, vertex_map, num_steiner_vertices,
is_polyhedron, face_recovery_ok). Prefer the `cdt.tetrahedralize` wrapper,
which names these.)");

    m.attr("OUTER") = (int)cdt::OUTER;
    m.attr("INNER") = (int)cdt::INNER;
}
