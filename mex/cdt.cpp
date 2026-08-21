// MATLAB mex interface to the CDT algorithm. See cdt.m for the help text.

#include "cdt.h"

#include <mex.h>

#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <ostream>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

namespace {

// Send std::cout / std::cerr to the MATLAB console instead of to whatever
// terminal MATLAB happens to have been launched from.
class MexStreamBuf : public std::streambuf {
protected:
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        mexPrintf("%.*s", (int)n, s);
        return n;
    }
    int overflow(int c = EOF) override {
        if (c != EOF) mexPrintf("%c", (char)c);
        return 1;
    }
};

// Redirects for the lifetime of the mex call.
class ConsoleRedirect {
public:
    ConsoleRedirect() : out_(std::cout.rdbuf(&buf_)), err_(std::cerr.rdbuf(&buf_)) {}
    ~ConsoleRedirect() { std::cout.rdbuf(out_); std::cerr.rdbuf(err_); }
private:
    MexStreamBuf buf_;
    std::streambuf* out_;
    std::streambuf* err_;
};

[[noreturn]] void fail(const std::string& msg) { throw std::runtime_error(msg); }

// Read any real numeric matrix as column major doubles, so that callers do not
// have to remember to cast their index arrays.
std::vector<double> read_numeric(const mxArray* a, const char* name) {
    if (mxIsComplex(a) || mxIsSparse(a) || !mxIsNumeric(a))
        fail(std::string(name) + " must be a real dense numeric matrix");
    if (mxGetNumberOfDimensions(a) != 2)
        fail(std::string(name) + " must be a matrix");

    const size_t n = mxGetNumberOfElements(a);
    std::vector<double> out(n);
    const void* p = mxGetData(a);
    switch (mxGetClassID(a)) {
#define CDT_READ_AS(id, type) \
    case id: for (size_t i = 0; i < n; i++) out[i] = (double)((const type*)p)[i]; break;
        CDT_READ_AS(mxDOUBLE_CLASS, double)
        CDT_READ_AS(mxSINGLE_CLASS, float)
        CDT_READ_AS(mxINT8_CLASS, int8_T)
        CDT_READ_AS(mxUINT8_CLASS, uint8_T)
        CDT_READ_AS(mxINT16_CLASS, int16_T)
        CDT_READ_AS(mxUINT16_CLASS, uint16_T)
        CDT_READ_AS(mxINT32_CLASS, int32_T)
        CDT_READ_AS(mxUINT32_CLASS, uint32_T)
        CDT_READ_AS(mxINT64_CLASS, int64_T)
        CDT_READ_AS(mxUINT64_CLASS, uint64_T)
#undef CDT_READ_AS
    default:
        fail(std::string(name) + " has an unsupported class");
    }
    return out;
}

// MATLAB is column major and 1-indexed; the algorithm wants row major and
// 0-indexed.
std::vector<uint32_t> read_indices(
    const mxArray* a, mwSize cols, uint32_t bound, const char* name) {
    if (mxGetN(a) != cols)
        fail(std::string(name) + " must be #" + name + " by " + std::to_string(cols));
    const std::vector<double> raw = read_numeric(a, name);
    const size_t rows = mxGetM(a);

    std::vector<uint32_t> out(rows * cols);
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++) {
            const double v = raw[i + j * rows];
            if (v != (double)(long long)v)
                fail(std::string(name) + " must contain integers");
            if (v < 1 || v > (double)bound)
                fail(std::string(name) + " contains index " + std::to_string((long long)v) +
                     " outside 1:" + std::to_string(bound));
            out[i * cols + j] = (uint32_t)v - 1;
        }
    return out;
}

bool read_flag(const mxArray* a, const char* name) {
    if (mxGetNumberOfElements(a) != 1 || (!mxIsLogicalScalar(a) && !mxIsNumeric(a)))
        fail(std::string(name) + " must be a logical scalar");
    return mxIsLogicalScalar(a) ? mxIsLogicalScalarTrue(a) : (mxGetScalar(a) != 0);
}

std::string read_string(const mxArray* a) {
    if (!mxIsChar(a)) fail("Parameter names must be strings");
    char* s = mxArrayToString(a);
    std::string out(s ? s : "");
    mxFree(s);
    return out;
}

// Column major mxArray out of row major data, with an optional offset applied
// to turn 0-indexed corners into 1-indexed ones.
template <typename T>
mxArray* to_mx(const std::vector<T>& data, size_t rows, size_t cols, double offset = 0) {
    mxArray* a = mxCreateDoubleMatrix((mwSize)rows, (mwSize)cols, mxREAL);
    double* p = mxGetPr(a);
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            p[i + j * rows] = (double)data[i * cols + j] + offset;
    return a;
}

void run(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 2) fail("Usage: [TV,TT,TL,I,S] = cdt(V,F,...)");
    if (nlhs > 5) fail("Too many output arguments");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0]) != 3)
        fail("V must be a real #V by 3 matrix of doubles");
    const size_t nv = mxGetM(prhs[0]);
    if (nv > UINT32_MAX) fail("V has too many rows");

    const double* vcm = mxGetPr(prhs[0]);
    std::vector<double> V(3 * nv);
    for (size_t i = 0; i < nv; i++)
        for (size_t k = 0; k < 3; k++) V[3 * i + k] = vcm[i + k * nv];

    const std::vector<uint32_t> F = read_indices(prhs[1], 3, (uint32_t)nv, "F");
    const size_t nf = mxGetM(prhs[1]);

    std::vector<uint32_t> E;
    size_t ne = 0;
    cdt::Options opts;

    for (int i = 2; i < nrhs; i++) {
        const std::string name = read_string(prhs[i]);
        if (++i >= nrhs) fail("Parameter '" + name + "' is missing a value");
        if (name == "Edges") {
            E = read_indices(prhs[i], 2, (uint32_t)nv, "Edges");
            ne = mxGetM(prhs[i]);
        } else if (name == "BoundingBox") {
            opts.bounding_box = read_flag(prhs[i], "BoundingBox");
        } else if (name == "InnerOnly") {
            opts.inner_only = read_flag(prhs[i], "InnerOnly");
        } else if (name == "Verbose") {
            opts.verbose = read_flag(prhs[i], "Verbose");
        } else {
            fail("Unknown parameter: " + name);
        }
    }

    const cdt::Result r = cdt::tetrahedralize(
        V.data(), (uint32_t)nv,
        F.data(), (uint32_t)nf,
        E.empty() ? nullptr : E.data(), (uint32_t)ne,
        opts);

    plhs[0] = to_mx(r.vertices, r.num_vertices(), 3);
    if (nlhs > 1) plhs[1] = to_mx(r.tets, r.num_tets(), 4, 1);
    if (nlhs > 2) plhs[2] = to_mx(r.labels, r.num_tets(), 1);
    if (nlhs > 3) plhs[3] = to_mx(r.vertex_map, r.vertex_map.size(), 1, 1);
    if (nlhs > 4) {
        const char* fields[] = {"num_steiner_vertices", "is_polyhedron", "face_recovery_ok"};
        plhs[4] = mxCreateStructMatrix(1, 1, 3, fields);
        mxSetField(plhs[4], 0, "num_steiner_vertices",
                   mxCreateDoubleScalar((double)r.num_steiner_vertices));
        mxSetField(plhs[4], 0, "is_polyhedron", mxCreateLogicalScalar(r.is_polyhedron));
        mxSetField(plhs[4], 0, "face_recovery_ok", mxCreateLogicalScalar(r.face_recovery_ok));
    }
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    std::string error;
    {
        ConsoleRedirect redirect;
        try {
            run(nlhs, plhs, nrhs, prhs);
        } catch (const std::exception& e) {
            error = e.what();
        } catch (...) {
            error = "cdt: unknown error";
        }
    }
    // Raised outside the try block so that every destructor has already run:
    // mexErrMsgIdAndTxt does not return.
    if (!error.empty()) mexErrMsgIdAndTxt("cdt:error", "%s", error.c_str());
}
