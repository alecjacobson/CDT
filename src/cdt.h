#ifndef CDT_H
#define CDT_H

// In-memory interface to the CDT algorithm.
//
// This is the entry point shared by the command line tool, the python
// bindings and the matlab mex bindings: it takes a PLC given as plain arrays
// and returns the tetrahedrization as plain arrays, without touching the file
// system.

#include <cstdint>
#include <cstddef>
#include <vector>

namespace cdt {

// Per-tetrahedron classification. Values match the DT_* macros in delaunay.h.
enum Label : uint8_t {
    UNKNOWN = 0,
    OUTER = 1,  // outside the input polyhedron
    INNER = 2   // inside the input polyhedron
};

struct Options {
    // Enclose the input in a box of eight extra vertices (-b). Without this
    // the tetrahedrization only covers the convex hull of the input.
    bool bounding_box = false;
    // Report progress on stdout (-v).
    bool verbose = false;
    // Drop the tets classified as OUTER from the output (-r). Vertices are
    // kept as they are, so the output may contain unreferenced vertices.
    bool inner_only = false;
};

struct Result {
    // Vertex positions, xyz interleaved. Steiner points are rounded to the
    // closest representable double, so a tet may be degenerate or inverted
    // even though the exact CDT is not.
    std::vector<double> vertices;
    // Tet corners, four indices into `vertices` per tet. INNER tets come
    // first, then OUTER ones; `labels` says which is which. Corners are
    // ordered so that det[b-a, c-a, d-a] > 0, i.e. positively oriented. Note
    // this differs from the corner order the command line tool writes to
    // .tet files.
    std::vector<uint32_t> tets;
    std::vector<uint8_t> labels; // one Label per tet

    // vertex_map[i] is the position of the i-th input vertex in `vertices`.
    // The input is deduplicated and reordered, so this is not the identity.
    std::vector<uint32_t> vertex_map;

    uint32_t num_input_vertices = 0;   // unique input vertices, plus 8 if bounding_box
    uint32_t num_steiner_vertices = 0; // vertices added to make the CDT exist
    // False if some edge of the input has an odd number of incident faces, in
    // which case inside/outside is meaningless and every tet is INNER.
    bool is_polyhedron = false;
    // False if face recovery had to fall back on the slower method.
    bool face_recovery_ok = false;

    size_t num_vertices() const { return vertices.size() / 3; }
    size_t num_tets() const { return tets.size() / 4; }
};

// Constrained Delaunay tetrahedrization of the PLC made of `nv` vertices
// (`V`, xyz interleaved), `nf` triangles (`F`, three indices into V each) and
// `ne` constraint segments (`E`, two indices into V each). `E` may be null
// when `ne` is 0. Validity of the PLC is assumed, not verified.
//
// Throws std::runtime_error on invalid input.
Result tetrahedralize(
    const double* V, uint32_t nv,
    const uint32_t* F, uint32_t nf,
    const uint32_t* E, uint32_t ne,
    const Options& opts);

} // namespace cdt

#endif
