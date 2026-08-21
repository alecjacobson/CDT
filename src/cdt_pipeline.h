#ifndef CDT_PIPELINE_H
#define CDT_PIPELINE_H

// Internal entry point: the CDT pipeline at the TetMesh level. Callers that
// only need plain arrays should use cdt.h instead.

#include <cstdint>

#include "delaunay.h"
#include "inputPLC.h"

struct SteinerCDTStats {
    bool is_polyhedron = false;
    bool face_recovery_ok = false;
    uint32_t num_steiner_vertices = 0;
    uint32_t num_inner_tets = 0;
};

// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!
// 'options' is a (possibly empty) string of characters, each controlling
// one option as follows:
// l: log results to cdt_log.csv
// b: add eight vertices to enclose everything in a box
// v: verbose mode
// w: log to screen
//
// The returned mesh is owned by the caller.
TetMesh* createSteinerCDT(inputPLC& plc, const char* options, SteinerCDTStats* stats = nullptr);

#endif
