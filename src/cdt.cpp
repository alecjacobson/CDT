#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#endif

#include "cdt.h"
#include "cdt_pipeline.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>

#include "PLC.h"
#include "logger.h"

TetMesh* createSteinerCDT(inputPLC& plc, const char* options, SteinerCDTStats* stats) {
    bool log = false, bbox = false, verbose = false, logscreen = false;

    for (size_t i = 0; i < strlen(options); i++) switch (options[i]) {
    case 'l':
        log = true; break;
    case 'b':
        bbox = true; break;
    case 'v':
        verbose = true; break;
    case 'w':
        logscreen = true; break;
    } // Just ignore unknown options

    if (bbox) plc.addBoundingBoxVertices();

    if (logscreen) {
        log = true;
        startLogging(NULL);
    }
    else if (log) startLogging(plc.input_file_name);

    // Build a delaunay tetrahedrization of the vertices
    TetMesh* tin = new TetMesh;
    tin->init_vertices(plc.coordinates.data(), plc.numVertices());
    tin->tetrahedrize();

    if (verbose) printf("DT of the vertices built\n");

    if (log) logTimeChunk();

    // Build a structured PLC linked to the Delaunay tetrahedrization
    PLCx Steiner_plc(
          *tin,
          plc.triangle_vertices.data(),
          plc.numTriangles(),
          plc.edge_vertices.data(),
          plc.numEdges());

    // Recover segments by inserting Steiner points in both the PLC and the tetrahedrization
    Steiner_plc.segmentRecovery_HSi(!verbose);

    if (log) logTimeChunk();

    // Recover PLC faces by locally remeshing the tetrahedrization
    bool sisMethodWorks = Steiner_plc.faceRecovery(!verbose);

    if (log) logTimeChunk();

    // Mark the tets which are bounded by the PLC.
    // If the PLC is not a valid polyhedron (i.e. it has odd-valency edges)
    // all the tets but the ghosts are marked as "internal".
    uint32_t num_inner_tets = (uint32_t)Steiner_plc.markInnerTets();

    if (log) logTimeChunk();

    if (log) {
        logMemInfo();
        logBoolean(Steiner_plc.is_polyhedron);
        logInteger(plc.numVertices());
        logInteger(Steiner_plc.input_nt);
        logInteger(Steiner_plc.numSteinerVertices());
        logInteger(tin->countNonGhostTets());
        logInteger(num_inner_tets);
        size_t nflip, nflat;
        tin->hasBadSnappedOrientations(nflip, nflat);
        logInteger((uint32_t)nflat);
        logInteger((uint32_t)nflip);
        logBoolean(sisMethodWorks);
        finishLogging();
    }

    if (stats) {
        stats->is_polyhedron = Steiner_plc.is_polyhedron;
        stats->face_recovery_ok = sisMethodWorks;
        stats->num_steiner_vertices = Steiner_plc.numSteinerVertices();
        stats->num_inner_tets = num_inner_tets;
    }

    return tin;
}

namespace cdt {

namespace {

// Same ordering as inputPLC.h's vertex_compare: it compares differences, so
// -0.0 and 0.0 land in the same bucket, exactly as the deduplication does.
int coord_compare(const double* a, const double* b) {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return 4 * ((dx > 0) - (dx < 0)) +
           2 * ((dy > 0) - (dy < 0)) +
           ((dz > 0) - (dz < 0));
}

// inputPLC deduplicates and spatially reorders its vertices without ever doing
// arithmetic on the coordinates, so every input vertex still appears verbatim
// in `plc.coordinates`. Recover the correspondence by looking the coordinates
// up rather than by threading an index through the reordering.
std::vector<uint32_t> build_vertex_map(
    const double* V, uint32_t nv, const std::vector<double>& out_coords) {
    const uint32_t n_out = (uint32_t)(out_coords.size() / 3);

    std::vector<uint32_t> sorted(n_out);
    std::iota(sorted.begin(), sorted.end(), 0u);
    const double* c = out_coords.data();
    std::sort(sorted.begin(), sorted.end(), [c](uint32_t a, uint32_t b) {
        return coord_compare(c + 3 * a, c + 3 * b) < 0;
    });

    std::vector<uint32_t> map(nv);
    for (uint32_t i = 0; i < nv; i++) {
        const double* p = V + 3 * i;
        const auto it = std::lower_bound(
            sorted.begin(), sorted.end(), p,
            [c](uint32_t a, const double* q) { return coord_compare(c + 3 * a, q) < 0; });
        if (it == sorted.end() || coord_compare(c + 3 * (*it), p) != 0)
            throw std::runtime_error("cdt: could not map input vertex " + std::to_string(i) +
                                     " onto an output vertex");
        map[i] = *it;
    }
    return map;
}

void validate(const double* V, uint32_t nv, const uint32_t* F, uint32_t nf,
              const uint32_t* E, uint32_t ne) {
    // Fewer than four vertices cannot span a tet, and inputPLC's deduplication
    // reads uninitialized memory when handed a single vertex.
    if (nv < 4) throw std::runtime_error("cdt: at least 4 input vertices are required");
    if (nf < 1) throw std::runtime_error("cdt: at least 1 input triangle is required");
    for (size_t i = 0; i < 3 * (size_t)nv; i++)
        if (!std::isfinite(V[i]))
            throw std::runtime_error("cdt: input vertex " + std::to_string(i / 3) +
                                     " is not finite");
    for (size_t i = 0; i < 3 * (size_t)nf; i++)
        if (F[i] >= nv)
            throw std::runtime_error("cdt: triangle " + std::to_string(i / 3) +
                                     " references out of range vertex " + std::to_string(F[i]));
    for (size_t i = 0; i < 2 * (size_t)ne; i++)
        if (E[i] >= nv)
            throw std::runtime_error("cdt: segment " + std::to_string(i / 2) +
                                     " references out of range vertex " + std::to_string(E[i]));
}

} // namespace

Result tetrahedralize(
    const double* V, uint32_t nv,
    const uint32_t* F, uint32_t nf,
    const uint32_t* E, uint32_t ne,
    const Options& opts) {
    initFPU();

    if (E == nullptr) ne = 0;
    validate(V, nv, F, nf, E, ne);

    // inputPLC::postProcess takes mutable pointers; keep the caller's arrays
    // untouched.
    std::vector<double> vcopy(V, V + 3 * (size_t)nv);
    std::vector<uint32_t> fcopy(F, F + 3 * (size_t)nf);
    std::vector<uint32_t> ecopy(E, E + 2 * (size_t)ne);

    inputPLC plc;
    plc.input_file_name = "";
    plc.postProcess(vcopy.data(), nv, fcopy.data(), nf, ecopy.data(), ne, opts.verbose);

    Result result;
    // Built before the bounding box vertices are appended, so that it only
    // ever refers to vertices the caller passed in.
    result.vertex_map = build_vertex_map(V, nv, plc.coordinates);

    std::string options;
    if (opts.bounding_box) options += 'b';
    if (opts.verbose) options += 'v';

    SteinerCDTStats stats;
    const std::unique_ptr<TetMesh> tin(createSteinerCDT(plc, options.c_str(), &stats));

    result.num_input_vertices = plc.numVertices();
    result.num_steiner_vertices = stats.num_steiner_vertices;
    result.is_polyhedron = stats.is_polyhedron;
    result.face_recovery_ok = stats.face_recovery_ok;

    const uint32_t n_out_v = tin->numVertices();
    result.vertices.resize(3 * (size_t)n_out_v);
    for (uint32_t i = 0; i < n_out_v; i++) {
        double* p = result.vertices.data() + 3 * (size_t)i;
        // apap: round to the closest representable double, matching saveTET.
        if (!tin->vertices[i]->getApproxXYZCoordinates(p[0], p[1], p[2], true))
            throw std::runtime_error("cdt: output vertex " + std::to_string(i) +
                                     " has no floating point approximation");
    }

    const uint32_t n_tets = tin->numTets();
    result.tets.reserve(4 * (size_t)n_tets);
    result.labels.reserve(n_tets);
    const auto emit = [&](uint32_t t) {
        const uint32_t* n = tin->tet_node.data() + 4 * (size_t)t;
        // CDT's internal corner order gives det[b-a, c-a, d-a] < 0. Swap the
        // last two so the output is positively oriented, which is what the
        // usual signed volume formulas (and the rest of the tooling around
        // these bindings) expect. The .tet files written by the command line
        // tool keep the internal order.
        const uint32_t flipped[4] = {n[0], n[1], n[3], n[2]};
        result.tets.insert(result.tets.end(), flipped, flipped + 4);
        result.labels.push_back((uint8_t)tin->mark_tetrahedra[t]);
    };
    for (uint32_t t = 0; t < n_tets; t++)
        if (tin->mark_tetrahedra[t] == DT_IN) emit(t);
    if (!opts.inner_only)
        for (uint32_t t = 0; t < n_tets; t++)
            if (!tin->isGhost(t) && tin->mark_tetrahedra[t] != DT_IN) emit(t);

    return result;
}

} // namespace cdt
