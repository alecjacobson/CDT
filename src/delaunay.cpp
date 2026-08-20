#include "delaunay.h"
#include <float.h>
#include <iomanip>

using namespace std;

void TetMesh::init_vertices(const double* coords, uint32_t num_v) {
    vertices.reserve(num_v);
    for (uint32_t i = 0; i < num_v; i++)
        vertices.push_back(new explicitPoint(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
    inc_tet.resize(num_v, UINT64_MAX);
    marked_vertex.resize(num_v, 0);
}

void TetMesh::init(uint32_t& unswap_k, uint32_t& unswap_l){
  const uint32_t n = numVertices();

  // Find non-coplanar vertices (we assume that no coincident vertices exist)
  int ori=0;
  uint32_t i=0, j=1, k=2, l=3;

  for (; ori == 0 && k < n - 1; k++)
      for (l = k + 1; ori == 0 && l < n; l++)
          ori = vOrient3D(i, j, k, l);

  l--; k--;

  if (ori == 0) {
      assert(0 && "TetMesh::init() - Input vertices do not define a volume");
      ip_error("TetMesh::init() - Input vertices do not define a volume.\n");
  }

  unswap_k = k;
  unswap_l = l;
  std::swap(vertices[k], vertices[2]); k=2;
  std::swap(vertices[l], vertices[3]); l=3;

  if(ori<0) std::swap(i, j); // Tets must have positive volume

  const uint32_t base_tet[] = { l, k, j, i, l, j, k, INFINITE_VERTEX, l, k, i, INFINITE_VERTEX, l, i, j, INFINITE_VERTEX, k, j, i, INFINITE_VERTEX };
  const uint64_t base_neigh[] = { 19, 15, 11, 7, 18, 10, 13, 3, 17, 14, 5, 2, 16, 6, 9, 1, 12, 8, 4, 0 };

  resizeTets(5);
  std::memcpy(getTetNodes(0), base_tet, 20 * sizeof(uint32_t));
  std::memcpy(getTetNeighs(0), base_neigh, 20 * sizeof(uint64_t));

  // set the vertex-(one_of_the)incident-tetrahedron relation
  inc_tet[i] = inc_tet[j] = inc_tet[k] = inc_tet[l] = 0;
}


void TetMesh::tetrahedrize() {
    uint32_t uk, ul;
    init(uk, ul); // First tet is made of vertices 0, 1, uk, ul

    // Need to unswap immediately to keep correct indexing and
    // ensure symbolic perturbation is coherent
    if (ul != 3) {
        std::swap(vertices[ul], vertices[3]);
        std::swap(inc_tet[ul], inc_tet[3]);
        for (uint32_t& tn : tet_node) if (tn == 3) tn = ul; else if (tn == ul) tn = 3;
    }

    if (uk != 2) {
        std::swap(vertices[uk], vertices[2]);
        std::swap(inc_tet[uk], inc_tet[2]);
        for (uint32_t& tn : tet_node) if (tn == 2) tn = uk; else if (tn == uk) tn = 2;
    }

    uint64_t ct = 0;
    for (uint32_t i = 2; i < numVertices(); i++) if (i != uk && i != ul) insertExistingVertex(i, ct);

    removeDelTets();
}


bool TetMesh::saveTET(const char* filename, bool inner_only) const
{
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveTET: Can't open file for writing.\n";
        return false;
    }

    f << numVertices() << " vertices\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    if (inner_only) {
        f << ngnt << " tets\n";
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    else {
        f << ngnt << " inner tets\n";
        f << countNonGhostTets()-ngnt << " outer tets\n";
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    
    f.close();

    return true;
}

bool TetMesh::saveVTU(const char* filename, bool inner_only) const
{
    ofstream f(filename);
    if (!f) {
        std::cerr << "\nTetMesh::saveVTU: Can't open file for writing.\n";
        return false;
    }

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    f << "  <UnstructuredGrid>\n";

    std::vector<uint32_t> tets_to_save;
    if (inner_only) {
        for (uint32_t i = 0; i < numTets(); i++) {
            if (mark_tetrahedra[i] == DT_IN) {
                tets_to_save.push_back(i);
            }
        }
    }
    else {
        for (uint32_t i = 0; i < numTets(); i++) {
            if (!isGhost(i)) {
                tets_to_save.push_back(i);
            }
        }
    }
    size_t num_t = tets_to_save.size();
    uint32_t num_v = numVertices();

    f << "    <Piece NumberOfPoints=\"" << num_v << "\" NumberOfCells=\"" << num_t << "\">\n";
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (uint32_t i = 0; i < num_v; i++) {
        f << "          " << *vertices[i] << "\n";
    }
    f << "        </DataArray>\n";
    f << "      </Points>\n";
    f << "      <Cells>\n";
    f << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (uint32_t i : tets_to_save) {
        const uint32_t* nodes = tet_node.data() + i * 4;
        f << "          " << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3] << "\n";
    }
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    f << "          ";
    for (uint32_t i = 1; i <= num_t; i++) {
        f << i * 4 << " ";
    }
    f << "\n        </DataArray>\n";
    f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    f << "          ";
    for (uint32_t i = 0; i < num_t; i++) {
        f << "10 "; // VTK_TETRA
    }
    f << "\n        </DataArray>\n";
    f << "      </Cells>\n";
    f << "    </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    f.close();
    return true;
}

bool TetMesh::saveMEDIT(const char* filename, bool inner_only) const
{
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveMEDIT: Can't open file for writing.\n";
        return false;
    }

    f << "MeshVersionFormatted 2\nDimension\n3\n";

    f << "Vertices\n" << numVertices() << "\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    if (inner_only) {
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << " 1\n";
        f << "Tetrahedra\n" << ngnt << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << tet_node[i * 4]+1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 1\n";
    }
    else {
        for (uint32_t i = 0; i < numVertices(); i++)
            f << *vertices[i] << " 1\n";
        f << "Tetrahedra\n" << countNonGhostTets() << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 1\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 2\n";
    }

    f.close();

    return true;
}


bool TetMesh::saveBinaryTET(const char* filename, bool inner_only) const
{
    ofstream f(filename, ios::binary);

    if (!f) {
        std::cerr << "\nTetMesh::saveBinaryTET: Can't open file for writing.\n";
        return false;
    }

    uint32_t num_v = numVertices(), num_t = 0;

    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) num_t++;

    f << num_v << " vertices\n";

    if (inner_only) {
        f << num_t << " tets\n";
    }
    else {
        f << num_t << " inner tets\n";
        f << countNonGhostTets() - num_t << " outer tets\n";
    }

    double c[3];
    for (uint32_t i = 0; i < numVertices(); i++) {
        vertices[i]->getApproxXYZCoordinates(c[0], c[1], c[2], true);
        f.write((const char*)(&c), sizeof(double) * 3);
    }

    const uint32_t* tnd = tet_node.data();

    if (inner_only) {
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
    }
    else {
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
    }

    f.close();

    return true;
}

bool TetMesh::saveBoundaryToOFF(const char* filename) const {
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveBoundaryToOFF: Can't open file for writing.\n";
        return false;
    }

    f << "OFF\n" << numVertices() << " ";

    size_t num_tris = 0;
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) num_tris++;

    f << num_tris << " 0\n";

    for (uint32_t i = 0; i < numVertices(); i++)
        f << *vertices[i] << "\n";

    uint32_t fv[3];
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) {
            getFaceVertices(i, fv);
            if (fv[0] == INFINITE_VERTEX || fv[1] == INFINITE_VERTEX || fv[2] == INFINITE_VERTEX) ip_error("Attempting to save skin of invalid In/Out classification.\n");
            f << "3 " << fv[0] << " " << fv[1] << " " << fv[2] << "\n";
        }
    f.close();

    return true;
}

bool TetMesh::saveRationalTET(const char* filename, bool inner_only)
{
#ifdef USE_INDIRECT_PREDS
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh::saveRationalTET: Can't open file for writing.\n";
        return false;
    }

    f << numVertices() << " vertices\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    if (inner_only) {
        f << ngnt << " tets\n";
        for (uint32_t i = 0; i < numVertices(); i++) {
            bigrational c[3];
            vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
            f << c[0] << " " << c[1] << " " << c[2] << "\n";
        }
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    else {
        f << ngnt << " inner tets\n";
        f << countNonGhostTets() - ngnt << " outer tets\n";
        for (uint32_t i = 0; i < numVertices(); i++) {
            bigrational c[3];
            vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
            f << c[0] << " " << c[1] << " " << c[2] << "\n";
        }
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }

    f.close();
#endif

    return true;
}

// Swap t and l while assuming that l is a valid tet (not to be deleted)
inline void TetMesh::moveDeletedToTail(uint64_t t, uint64_t l) {
    const uint64_t t2 = t >> 2, l2 = l >> 2;
    uint32_t* tnt = tet_node.data() + t, *lnt = tet_node.data() + l;
    uint64_t* tnn = tet_neigh.data() + t, *lnn = tet_neigh.data() + l;
    const uint32_t* te = tnt + 4;

    do {
        *tnt++ = *lnt;

        uint64_t neigh = *lnn++;
        *tnn++ = neigh;
        tet_neigh[neigh] = t++;

        if (*lnt != INFINITE_VERTEX && inc_tet[*lnt] == l2)
            inc_tet[*lnt] = t2;
        lnt++;
    } while (tnt != te);

    mark_tetrahedra[t2] = mark_tetrahedra[l2];
}

void TetMesh::removeDelTets() {
    uint64_t* dp = Del_deleted.data();
    uint64_t* de = dp + Del_deleted.size();
    uint64_t last_valid = numTets() - 1;
    while (isToDeleteSmall(last_valid)) last_valid--;

    while (dp != de) {
        if (*dp < (last_valid<<2)) {
            moveDeletedToTail(*dp, (last_valid << 2));
            last_valid--;
            while (isToDeleteSmall(last_valid)) last_valid--;
        }
        dp++;
    }

    resizeTets(++last_valid);
    Del_deleted.clear();
}

bool TetMesh::tetHasVertex(uint64_t t, uint32_t v) const {
    t <<= 2;
    return tet_node[t] == v || tet_node[t + 1] == v || tet_node[t + 2] == v || tet_node[t + 3] == v;
}

bool TetMesh::tetContainsPoint(uint64_t tet, const pointType* p) const {
    const uint32_t* Node = tet_node.data() + tet;
    if (Node[3] == INFINITE_VERTEX) return false;

    const pointType* v1 = vertices[Node[0]];
    const pointType* v2 = vertices[Node[1]];
    const pointType* v3 = vertices[Node[2]];
    const pointType* v4 = vertices[Node[3]];

    if (genericPoint::orient3D(*p, *v2, *v3, *v4) > 0) return false;
    if (genericPoint::orient3D(*v1, *p, *v3, *v4) > 0) return false;
    if (genericPoint::orient3D(*v1, *v2, *p, *v4) > 0) return false;
    if (genericPoint::orient3D(*v1, *v2, *v3, *p) > 0) return false;
    return true;
}

void TetMesh::oppositeTetEdge(const uint64_t tet, const uint32_t v[2], uint32_t ov[2]) const {
    int i = 0, j = 0;
    while (i < 4) {
        const uint32_t w = tet_node[tet + i];
        if (w != v[0] && w != v[1]) ov[j++] = w;
        i++;
    }
    assert(j == 2);
}

uint64_t TetMesh::getCornerFromOppositeTet(uint64_t t, uint64_t n) const {
    t <<= 2;
    for (int i = 0; i < 4; i++)
        if ((tet_neigh[t + i] >> 2) == n)
            return tet_neigh[t + i];
    assert(0);
    return UINT64_MAX;
}

void TetMesh::getFaceVertices(uint64_t t, uint32_t v[3]) const {
    uint64_t tv = t & 3;
    const uint32_t* Node = tet_node.data() + (t - tv);
    v[0] = Node[(++tv) & 3];
    v[1] = Node[(++tv) & 3];
    v[2] = Node[(++tv) & 3];
}

void TetMesh::getFaceSortedVertices(uint64_t t, uint32_t v[3]) const {
    getFaceVertices(t, v);
    if (v[0] > v[1]) std::swap(v[0], v[1]);
    if (v[1] > v[2]) std::swap(v[1], v[2]);
    if (v[0] > v[1]) std::swap(v[0], v[1]);
}

bool TetMesh::getTetsFromFaceVertices(uint32_t v1, uint32_t v2, uint32_t v3, uint64_t* nt) const {
    static std::vector<uint64_t> vt; // Static to avoid reallocation at each call
    VTfull(v1, vt);
    int i = 0;
    for (uint64_t t : vt) if (tetHasVertex(t, v2) && tetHasVertex(t, v3)) nt[i++] = t;
    vt.clear();
    return (i == 2);
}

uint64_t TetMesh::tetOppositeCorner(uint64_t t, uint32_t v1, uint32_t v2, uint32_t v3) const {
    const uint64_t tb = t << 2;
    const uint32_t* n = tet_node.data() + tb;
    for (int i = 0; i < 3; i++)
        if (n[i] != v1 && n[i] != v2 && n[i] != v3)
            return tet_neigh[tb + i];
    assert(n[3] != v1 && n[3] != v2 && n[3] != v3);
    return tet_neigh[tb + 3];
}

void TetMesh::resizeTets(uint64_t new_size) {
    mark_tetrahedra.resize(new_size, 0);
    new_size <<= 2;
    tet_node.resize(new_size);
    tet_neigh.resize(new_size);
}

void TetMesh::reserveTets(uint64_t new_capacity) {
    mark_tetrahedra.reserve(new_capacity);
    new_capacity <<= 2;
    tet_node.reserve(new_capacity);
    tet_neigh.reserve(new_capacity);
}

uint64_t TetMesh::searchTetrahedron(uint64_t tet, const uint32_t v_id)
{
    if (tet_node[tet + 3] == INFINITE_VERTEX)
        tet = getIthNeighbor(getTetNeighs(tet), 3);

    uint64_t i, f0 = 4;
    do {
        const uint32_t* Node = getTetNodes(tet);
        if (Node[3] == INFINITE_VERTEX) return tet;

        const uint64_t* Neigh = getTetNeighs(tet);
        for (i = 0; i < 4; i++)
            if (i != f0 && vOrient3D(Node[tetON1(i)], Node[tetON2(i)], Node[tetON3(i)], v_id) < 0) {
                tet = getIthNeighbor(Neigh, i);
                f0 = Neigh[i] & 3;
                break;
            }
    } while (i != 4);

    return tet;
}


int TetMesh::symbolicPerturbation(uint32_t indices[5]) const {
    int swaps = 0;
    int n = 5;
    int count;
    do {
        count = 0;
        n--;
        for (int i = 0; i < n; i++) {
            if (indices[i] > indices[i + 1]) {
                std::swap(indices[i], indices[i + 1]);
                count++;
            }
        }
        swaps += count;
    } while (count);

    n = vOrient3D(indices[1], indices[2], indices[3], indices[4]);
    if (n) return (swaps % 2) ? (-n) : n;

    n = vOrient3D(indices[0], indices[2], indices[3], indices[4]);
    return (swaps % 2) ? (n) : (-n);
}

int TetMesh::vertexInTetSphere(const uint32_t Node[4], uint32_t v_id) const {
    int det = vInSphere(Node[0], Node[1], Node[2], Node[3], v_id);
    if (det) return det;
    uint32_t nn[5] = { Node[0],Node[1],Node[2],Node[3],v_id };
    det = symbolicPerturbation(nn);
    if (det == 0.0) {
        std::cout << *vertices[Node[0]] << " (ID: " << Node[0] << ")\n";
        std::cout << *vertices[Node[1]] << " (ID: " << Node[1] << ")\n";
        std::cout << *vertices[Node[2]] << " (ID: " << Node[2] << ")\n";
        std::cout << *vertices[Node[3]] << " (ID: " << Node[3] << ")\n";
        std::cout << *vertices[v_id] << " (ID: " << v_id << ")\n";
        assert(0 && "Symbolic perturbation failed! Should not happen");
        ip_error("Symbolic perturbation failed! Should not happen.\n");
    }
    return det;
}

int TetMesh::vertexInTetSphere(uint64_t tet, uint32_t v_id) const
{
  const uint32_t* Node = getTetNodes(tet);
  int det;

  if (Node[3] == INFINITE_VERTEX) {
      if ((det = vOrient3D(Node[0], Node[1], Node[2], v_id)) != 0) return det;
      const uint32_t nn[4] = {Node[0], Node[1], Node[2], tet_node[tet_neigh[tet + 3]]};
      return -vertexInTetSphere(nn, v_id);
  }
  else return vertexInTetSphere(Node, v_id);
}

#ifdef USE_MAROTS_METHOD
void TetMesh::deleteInSphereTets(uint64_t tet, const uint32_t v_id)
{
  pushAndMarkDeletedTets(tet);

  for(uint64_t t = Del_deleted.size() - 1; t < Del_deleted.size(); t++) {
    uint64_t tet = Del_deleted[t];
    uint64_t* Neigh = getTetNeighs(tet);
    uint32_t* Node = getTetNodes(tet);

    uint64_t neigh = getIthNeighbor(Neigh, 0);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 1);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 2);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 3);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0){
        if(Node[1]<Node[2])
          bnd_push(v_id, Node[0], Node[2], Node[1], Neigh[3]);
        else
          bnd_push(v_id, Node[1], Node[0], Node[2], Neigh[3]);
      }
      else pushAndMarkDeletedTets(neigh);
    }
  }
}


void TetMesh::tetrahedrizeHole(uint64_t* tet){
  uint64_t clength = Del_deleted.size(); // Num tets removed
  uint64_t blength = numDelTmp(); // Num tets to insert

  uint64_t tn = numTets();

  if(blength > clength){
    for (uint64_t i = clength; i<blength; i++, tn++)
        Del_deleted.push_back(tn<<2);

    clength = blength;
    resizeTets(tn);
  }

  uint64_t start = clength - blength;

  for (uint64_t i=0; i<blength; i++)
  {
    const uint64_t tet = Del_deleted[i + start];
    uint32_t* Node = getTetNodes(tet);

    Node[0] = Del_tmp[i].node[0];
    Node[1] = Del_tmp[i].node[1];
    Node[2] = Del_tmp[i].node[2];
    Node[3] = Del_tmp[i].node[3];

    uint64_t bnd = Del_tmp[i].bnd;
    tet_neigh[tet] = bnd;
    tet_neigh[bnd] = tet;
    Del_tmp[i].bnd = tet;

    mark_tetrahedra[tet >> 2] = 0;

    if(tet_node[tet+3]!=INFINITE_VERTEX)
      for(uint32_t j=0; j<4; j++)
          inc_tet[tet_node[tet + j]] = tet>>2;
  }

  uint64_t tlength = 0;
  const uint64_t middle = blength * 3 / 2;

  uint64_t* Tmp = delTmpVec();
  const unsigned index[4] = { 2,3,1,2 };

  for (uint64_t i = 0; i < blength; i++)
  {
      uint64_t tet = Del_deleted[start + i];
      const uint32_t* Node = getTetNodes(tet);

      for (uint64_t j = 0; j < 3; j++)
      {
          uint64_t key = ((uint64_t)Node[index[j]] << 32) + Node[index[j + 1]];
          tet++;

          uint64_t k;
          for (k = 0; k < tlength; k++) if (Tmp[k] == key) break;

          if (k == tlength) {
              Tmp[tlength] = (key >> 32) + (key << 32);
              Tmp[middle + tlength] = tet;
              tlength++;
          }
          else {
              uint64_t pairValue = Tmp[middle + k];
              tet_neigh[tet] = pairValue;
              tet_neigh[pairValue] = tet;
              tlength--;
              if (k < tlength) {
                  Tmp[k] = Tmp[tlength];
                  Tmp[middle + k] = Tmp[middle + tlength];
              }
          }
      }
  }

  flushDelTmp();
  *tet = Del_deleted[start];
  Del_deleted.resize(start);
}

void TetMesh::insertExistingVertex(const uint32_t vi, uint64_t& ct)
{
    ct = searchTetrahedron(ct, vi);
    deleteInSphereTets(ct, vi);
    tetrahedrizeHole(&ct);
    uint64_t lt = ct;
    if (tet_node[lt + 3] == INFINITE_VERTEX) lt = tet_neigh[lt + 3];
    inc_tet[vi] = lt >> 2;
}

#else

// Expand by adjacencies to collect all tets whose circumsphere contains v_id
void TetMesh::getUnconstrainedCavity(const uint32_t v_id, const uint64_t tet, std::vector<uint64_t>& cavityCorners) {
    uint64_t ntet = searchTetrahedron(tet, v_id) >> 2;

    const uint32_t* tet_node_data = tet_node.data();
    const uint64_t* tet_neigh_data = tet_neigh.data();

    size_t first = Del_deleted.size();
    pushAndMarkDeletedTets(ntet << 2);

    for (size_t i = first; i < Del_deleted.size(); i++) {
        const uint64_t* nb = tet_neigh_data + Del_deleted[i];
        const uint64_t* nl = nb + 4;

        for (; nb < nl; nb++)
        {
            const uint64_t n0 = *nb >> 2;

            if (IEV_IS_UNVISITED(n0)) {
                if (vertexInTetSphere(n0 << 2, v_id) < 0) {
                    IEV_MARK_VISITED_TWICE(n0);
                    cavityCorners.push_back(*nb);
                }
                else pushAndMarkDeletedTets(n0 << 2);
            }
            else if (IEV_IS_VISITED_TWICE(n0)) cavityCorners.push_back(*nb);
        }
    }
}

// Expand by adjacencies to collect all tets whose circumsphere contains v_id.
// Epansion proceeds only across unconstrained facets.
void TetMesh::getConstrainedCavity(const uint32_t v_id, uint64_t& tet, std::vector<uint64_t>& cavity, uint32_t& cv1, uint32_t& cv2, uint32_t& cv3) {
    tet = searchTetrahedron(tet, v_id) & (~3);

    const uint32_t* tet_node_data = tet_node.data();
    const uint64_t* tet_neigh_data = tet_neigh.data();

    // Init cavity with all tets containing v_id (one, two or more)
    const pointType& v_p = *vertices[v_id];
    static std::vector<uint64_t> et;
    et.clear();
    const uint32_t* tn = tet_node_data + tet;
    size_t i;

    for (i = 0; i < 4; i++) {
        uint64_t opp = tet_neigh[tet + i] & (~3);
        uint32_t fv[3] = { tn[tetN1(i)] , tn[tetN2(i)] , tn[tetN3(i)] };
        if (fv[0] != INFINITE_VERTEX && fv[1] != INFINITE_VERTEX && fv[2] != INFINITE_VERTEX && vOrient3D(v_id, fv[0], fv[1], fv[2]) == 0) {
            if (!genericPoint::misaligned(*vertices[fv[0]], v_p, *vertices[fv[1]])) { cv1 = fv[0]; cv2 = fv[1]; ETfull(cv1, cv2, et); } // On edge v1-v2
            else if (!genericPoint::misaligned(*vertices[fv[1]], v_p, *vertices[fv[2]])) { cv1 = fv[1]; cv2 = fv[2]; ETfull(cv1, cv2, et); } // On edge v2-v3
            else if (!genericPoint::misaligned(*vertices[fv[2]], v_p, *vertices[fv[0]])) { cv1 = fv[2]; cv2 = fv[0]; ETfull(cv1, cv2, et); } // On edge v3-v1
            else { // On face
                cv1 = fv[0]; cv2 = fv[1]; cv3 = fv[2];
                cavity.push_back(tet); markToDelete(tet);
                cavity.push_back(opp); markToDelete(opp);
            }
            for (uint64_t y : et) { cavity.push_back(y<<2); markToDelete(y<<2); }
            break;
        }
    }
    if (i == 4) { cavity.push_back(tet); markToDelete(tet); } // Not on any face or edge

    // Expand cavity
    for (size_t i = 0; i < cavity.size(); i++) {
        const uint32_t mm = IEV_GET_IN_OR_OUT((cavity[i] >> 2));
        const uint64_t* nb = tet_neigh_data + cavity[i];
        const uint64_t* nl = nb + 4;
        for (; nb < nl; nb++)
        {
            uint64_t n0 = *nb >> 2;
            if (IEV_IS_UNVISITED(n0)) {
                if (mm != IEV_GET_IN_OR_OUT(n0) || vertexInTetSphere(n0 << 2, v_id) < 0) IEV_MARK_VISITED_TWICE(n0);
                else { n0 <<= 2; cavity.push_back(n0); markToDelete(n0); }
            }
        }
    }

    // Unmark outer corners
    for (uint64_t c : cavity) {
        const uint64_t* nb = tet_neigh_data + c;
        const uint64_t* nl = nb + 4;
        for (; nb < nl; nb++) if (!isToDelete(*nb)) IEV_MARK_UNVISITED((*nb >> 2));
    }

    for (uint64_t c : cavity) unmarkToDelete(c);
}

void TetMesh::undeleteTets(size_t num_deltets) {
    for (size_t i = Del_deleted.size() - num_deltets; i < Del_deleted.size(); i++)
        mark_tetrahedra[Del_deleted[i] >> 2] &= (~((uint32_t)1073741824)); // No longer to delete
    Del_deleted.resize(Del_deleted.size() - num_deltets);
}

void TetMesh::getCavityConnectivity(const std::vector<uint64_t>& cavityCorners, std::vector<uint64_t>& adjs) const {
    const uint32_t* tet_node_data = tet_node.data();
    const uint64_t* tet_neigh_data = tet_neigh.data();
    const uint32_t* tn;
    uint32_t v1, v2, v3, v[5];
    uint64_t b, c, c0, i, j = 0;

    adjs.resize(cavityCorners.size() * 3);

    for (const uint64_t t : cavityCorners) {
        c0 = tet_neigh_data[t];
        b = c0 & 3;
        tn = tet_node_data + c0 - b;
        v[0] = v[3] = tn[tetON1(b)];
        v[1] = v[4] = tn[tetON2(b)];
        v[2] = tn[tetON3(b)];

        for (i = 0; i < 3; i++) {
            v1 = v[i]; v2 = v[i + 1];
            if (v1 < v2) {
                c = c0;
                do {
                    v3 = tet_node_data[c];
                    c &= (~3);
                    tn = tet_node_data + c;
                    while (*tn == v1 || *tn == v2 || *tn == v3) tn++;
                    c = tet_neigh_data[tn - tet_node_data];
                } while (isToDelete(c));
                adjs[j++] = c;
                adjs[j++] = t;
            }
        }
    }
    assert(adjs.size() == j);
}

void TetMesh::retetrahedrizeCavity(const uint32_t v_id, const std::vector<uint64_t>& cavityCorners, const std::vector<uint64_t>& adjs, const uint32_t mark) {
    static const int fi[4][3] = { {2, 1, 3} ,{0, 2, 3} ,{1, 0, 3} ,{0, 1, 2} };
    uint32_t v1;
    uint32_t* tet_node_data = tet_node.data();
    uint64_t* tet_neigh_data = tet_neigh.data();

    // Resize the mesh to host the new tets
    uint64_t ntb, newpos = tet_node.size();
    if (cavityCorners.size() > Del_deleted.size()) {
        resizeTets(numTets() + (cavityCorners.size() - Del_deleted.size()));
        tet_node_data = tet_node.data();
        tet_neigh_data = tet_neigh.data();
    }

    // Create the new tets
    for (const uint64_t c : cavityCorners) {
        IEV_MARK_UNVISITED((c >> 2));
        if (Del_deleted.empty()) {
            ntb = newpos;
            newpos += 4;
        }
        else {
            ntb = Del_deleted.back();
            Del_deleted.pop_back();
        }
        const uint64_t cb = c & 3;
        const uint32_t* cr = tet_node_data + (c - cb);
        uint32_t* cn = tet_node_data + ntb;

        assert(cr[fi[cb][2]] == INFINITE_VERTEX || vOrient3D(v_id, cr[fi[cb][0]], cr[fi[cb][1]], cr[fi[cb][2]]) > 0);

        *cn++ = v_id;
        *cn++ = cr[fi[cb][0]];
        *cn++ = cr[fi[cb][1]];
        *cn++ = cr[fi[cb][2]];

        tet_neigh_data[ntb] = c; tet_neigh_data[c] = ntb;

        ntb >>= 2;
        if ((*(--cn)) != INFINITE_VERTEX) {
            inc_tet[*cn] = ntb;
            inc_tet[*(--cn)] = ntb;
            inc_tet[*(--cn)] = ntb;
            inc_tet[v_id] = ntb;
        }
        mark_tetrahedra[ntb] = mark;
    }

    // Restore inner connectivity
    for (size_t i = 0; i < adjs.size(); ) {
        uint64_t c = tet_neigh_data[adjs[i++]] & (~3); c++;
        uint64_t n = tet_neigh_data[adjs[i++]] & (~3); n++;
        uint64_t k = c;

        v1 = tet_node_data[c];
        if (v1 == tet_node_data[n + 1] || v1 == tet_node_data[n + 2] || v1 == tet_node_data[n]) v1 = tet_node_data[++c];
        if (v1 == tet_node_data[n + 1] || v1 == tet_node_data[n + 2] || v1 == tet_node_data[n]) ++c;
        v1 = tet_node_data[n];
        if (v1 == tet_node_data[k + 1] || v1 == tet_node_data[k + 2] || v1 == tet_node_data[k]) v1 = tet_node_data[++n];
        if (v1 == tet_node_data[k + 1] || v1 == tet_node_data[k + 2] || v1 == tet_node_data[k]) ++n;

        tet_neigh_data[c] = n; tet_neigh_data[n] = c;
    }
}

// Collect all tets whose circumsphere contains v_id and replace them
// with a star of new tets originating at v_id

void TetMesh::insertExistingVertex(const uint32_t v_id, uint64_t& tet)
{
    static std::vector<uint64_t> cavityCorners, adjs; // Static to avoid reallocation on each call
    getUnconstrainedCavity(v_id, tet, cavityCorners);
    getCavityConnectivity(cavityCorners, adjs);
    retetrahedrizeCavity(v_id, cavityCorners, adjs, 0);

    tet = tet_neigh[cavityCorners.back()];

    cavityCorners.clear();
}

#endif

void TetMesh::VT(uint32_t v, std::vector<uint64_t>& vt) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v];

    vt_queue.push_back(tetCornerAtVertex(t << 2, v));
    mark_Tet_31(t);

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb+j)&3];
            const uint64_t tbb = tb >> 2;
            if (tet_node[tb] != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                const uint64_t nt = tetCornerAtVertex(tb & (~3), v);
                vt_queue.push_back(nt); 
                mark_Tet_31(tbb); 
            }
        }
    }

    for (uint64_t t : vt_queue) {
        t >>= 2;
        unmark_Tet_31(t);
        vt.push_back(t);
    }

    vt_queue.clear();
}

void TetMesh::VV(uint32_t v, std::vector<uint32_t>& vv) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v];
    const uint64_t tb = t << 2;

    const uint64_t s = tetCornerAtVertex(tb, v);
    vt_queue.push_back(s);
    mark_Tet_31(t);

    const uint32_t* tn = tet_node.data() + tb;
    const uint64_t sb = s & 3;
    for (int j = 1; j < 4; j++) {
        const uint32_t w = tn[(sb + j) & 3];
        marked_vertex[w] |= 128;
        vv.push_back(w);
    }

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb + j) & 3];
            const uint64_t tbb = tb >> 2;
            const uint32_t w = tet_node[tb];
            if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                vt_queue.push_back(tetCornerAtVertex(tb & (~3), v));
                mark_Tet_31(tbb);
                if (!(marked_vertex[w] & 128)) {
                    marked_vertex[w] |= 128;
                    vv.push_back(w);
                }
            }
        }
    }

    for (uint64_t t : vt_queue) unmark_Tet_31(t>>2);
    vt_queue.clear();
    for (uint32_t w : vv) marked_vertex[w] &= 127;
}

void TetMesh::ET(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    VT(v1, et);
    for (size_t i = 0; i < et.size();)
        if (!tetHasVertex(et[i], v2)) {
            std::swap(et[i], et[et.size() - 1]);
            et.pop_back();
        }
        else i++;
}

void TetMesh::ETfull(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    VTfull(v1, et);
    for (size_t i = 0; i < et.size();)
        if (!tetHasVertex(et[i], v2)) {
            std::swap(et[i], et[et.size() - 1]);
            et.pop_back();
        }
        else i++;
}

void TetMesh::ETcorners(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    uint64_t t;
    VTfull(v1, et);
    for (uint64_t s : et) if (tetHasVertex(s, v2)) { t = (s<<2); break; }
    
    while (tet_node[t] == v1 || tet_node[t] == v2) t++;

    et.clear();

    uint64_t c0 = t;
    do {
        et.push_back(t); // Add tet
        uint64_t oc = tet_neigh[t] & (~3); // Get next base
        uint32_t cv = tet_node[t];
        t &= (~3);
        while (tet_node[t] == v1 || tet_node[t] == v2 || tet_node[t] == cv) t++;
        t = tetCornerAtVertex(oc, tet_node[t]); // Get corresp corner at opposite tet
    } while (t != c0);
}

void TetMesh::VTfull(uint32_t v, std::vector<uint64_t>& vt) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t s, t = inc_tet[v];
    vt_queue.push_back(t);
    mark_Tet_31(t);

    while (!vt_queue.empty()) {
        t = vt_queue.back();
        vt_queue.pop_back();
        vt.push_back(t);
        t <<= 2;
        s = tet_neigh[t] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 1] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 2] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 3] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
    }

    for (uint64_t t : vt) unmark_Tet_31(t);
}


bool TetMesh::hasEdge(uint32_t v1, uint32_t v2) const {
    static std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call

    //vt_queue.clear();
    //VT(v1, vt_queue);
    //for (uint64_t t : vt_queue) if (tetHasVertex(t, v2)) return true;
    //return false;

    uint64_t t = inc_tet[v1];
    const uint64_t tb = t << 2;
    if (tet_node[tb] == v2 || tet_node[tb + 1] == v2 || tet_node[tb + 2] == v2 || tet_node[tb + 3] == v2) return true;

    vt_queue.push_back(tetCornerAtVertex(tb, v1));
    mark_Tet_31(t);

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb + j) & 3];
            const uint64_t tbb = tb >> 2;
            const uint32_t w = tet_node[tb];
            if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                vt_queue.push_back(tetCornerAtVertex(tbb << 2, v1));
                mark_Tet_31(tbb);
                if (w == v2) {
                    for (uint64_t t : vt_queue) unmark_Tet_31(t >> 2);
                    vt_queue.clear();
                    return true;
                }
            }
        }
    }

    for (uint64_t t : vt_queue) unmark_Tet_31(t >> 2);
    vt_queue.clear();
    return false;
}


void TetMesh::swapTets(const uint64_t t1, const uint64_t t2) 
{
    if (t1 == t2) return;

    const uint64_t t1_id = t1<<2;
    const uint64_t t2_id = t2<<2;

    // update VT base relation
    for (int i = 0; i < 3; i++) if (inc_tet[tet_node[t1_id + i]] == t1) inc_tet[tet_node[t1_id + i]] = t2;
    if (tet_node[t1_id + 3] != INFINITE_VERTEX && inc_tet[tet_node[t1_id + 3]] == t1) inc_tet[tet_node[t1_id + 3]] = t2;

    for (int i = 0; i < 3; i++) if (inc_tet[tet_node[t2_id + i]] == t2) inc_tet[tet_node[t2_id + i]] = t1;
    if (tet_node[t2_id + 3] != INFINITE_VERTEX && inc_tet[tet_node[t2_id + 3]] == t2) inc_tet[tet_node[t2_id + 3]] = t1;

    // Update nodes and marks
    for (int i = 0; i < 4; i++) std::swap(tet_node[t1_id + i], tet_node[t2_id + i]);
    std::swap(mark_tetrahedra[t1], mark_tetrahedra[t2]);

    // update neigh-neigh relations
    const uint64_t ng1[] = { tet_neigh[t1_id + 0], tet_neigh[t1_id + 1], tet_neigh[t1_id + 2], tet_neigh[t1_id + 3] };
    const uint64_t ng2[] = { tet_neigh[t2_id + 0], tet_neigh[t2_id + 1], tet_neigh[t2_id + 2], tet_neigh[t2_id + 3] };

    for (int i = 0; i < 4; i++) if ((ng2[i] >> 2) != t1) tet_neigh[ng2[i]] = t1_id + i;
    for (int i = 0; i < 4; i++) if ((ng1[i] >> 2) != t2) tet_neigh[ng1[i]] = t2_id + i;

    for (int i = 0; i < 4; i++)
        if ((ng2[i] >> 2) != t1) tet_neigh[t1_id + i] = tet_neigh[t2_id + i];
        else tet_neigh[t1_id + i] = (tet_neigh[t2_id + i] & 3) + (t2 << 2);

    for (int i = 0; i < 4; i++)
        if ((ng1[i] >> 2) != t2) tet_neigh[t2_id + i] = ng1[i];
        else tet_neigh[t2_id + i] = (ng1[i] & 3) + (t1 << 2);
}

size_t TetMesh::markInnerTets(uint64_t single_start) {
    std::vector<uint64_t> C;

    // All ghosts are DT_OUT
    for (size_t i = 0; i < numTets(); i++)
        mark_tetrahedra[i] = (isGhost(i)) ? DT_OUT : DT_UNKNOWN;

    if (single_start != UINT64_MAX) C.push_back(single_start);
    else for (size_t i = 0; i < numTets(); i++)
        if (mark_tetrahedra[i] == DT_OUT) C.push_back(i);

    for (size_t i = 0; i < C.size(); i++) {
        uint64_t t = C[i];
        for (int j = 0; j < 4; j++) {
            const uint64_t n = tet_neigh[t * 4 + j];
            const uint64_t n2 = n >> 2;
            if (mark_tetrahedra[n2] == DT_UNKNOWN) {
                if (!cornerMask[n]) {
                    mark_tetrahedra[n2] = mark_tetrahedra[t];
                }
                else {
                    mark_tetrahedra[n2] = ((mark_tetrahedra[t] == DT_IN) ? (DT_OUT) : (DT_IN));
                }
                C.push_back(n2);
            }
        }
    }

    return std::count(mark_tetrahedra.begin(), mark_tetrahedra.end(), DT_IN);
}

void TetMesh::markInnerTetsNonManifold() {
    std::vector<uint64_t> C;

    // All ghosts are DT_OUT
    for (size_t i = 0; i < numTets(); i++)
        mark_tetrahedra[i] = (isGhost(i)) ? DT_OUT : DT_IN;

    for (size_t i = 0; i < numTets(); i++)
        if (mark_tetrahedra[i] == DT_OUT) C.push_back(i);

    for (size_t i = 0; i < C.size(); i++) {
        uint64_t t = C[i];
        for (int j = 0; j < 4; j++) {
            const uint64_t n = tet_neigh[t * 4 + j];
            const uint64_t n2 = n >> 2;
            if (mark_tetrahedra[n2] == DT_IN && !cornerMask[n]) {
                mark_tetrahedra[n2] = DT_OUT;
                C.push_back(n2);
            }
        }
    }
}

bool TetMesh::hasBadSnappedOrientations(size_t& num_flipped, size_t& num_flattened) const {
    const uint32_t* tn = tet_node.data();
    const uint32_t* end = tn + tet_node.size();
    num_flipped = num_flattened = 0;
    explicitPoint v[4];
    while (tn < end) {
        if (tn[3] != INFINITE_VERTEX) {
            for (int i = 0; i < 4; i++) {
                const pointType* p = vertices[tn[i]];
                if (p->isExplicit3D()) v[i] = p->toExplicit3D();
                else p->apapExplicit(v[i]);
            }
            const int o = pointType::orient3D(v[0], v[1], v[2], v[3]);
            if (o > 0) num_flipped++;
            else if (o == 0) num_flattened++;
        }
        tn += 4;
    }

    return (num_flipped || num_flattened);
}

void TetMesh::checkMesh(bool checkDelaunay) const {
    size_t i;
    const uint32_t num_vertices = (uint32_t)vertices.size();
    // Check tet nodes	
    for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
        const uint32_t* tn = tet_node.data() + i * 4;
        if (tn[0] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[1] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[2] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[3] != INFINITE_VERTEX && tet_node[i * 4 + 3] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[0] == tn[1] || tn[0] == tn[2] || tn[0] == tn[3]
            || tn[1] == tn[2] || tn[1] == tn[3] || tn[2] == tn[3]) 
            assert(0 && "Wrong tet node indexes!\n");
    }

    // Check neighbors	
    for (i = 0; i < numTets() * 4; i++) if (!isToDelete(i))
        if (tet_neigh[i] >= tet_neigh.size() || tet_neigh[tet_neigh[i]] != i)
            assert(0 && "Wrong neighbor!\n");

    // Check neighbor-node coherence
    for (i = 0; i < numTets() * 4; i++) if (!isToDelete(i)) {
        if (tetHasVertex(tet_neigh[i] >> 2, tet_node[i]))
            assert(0 && "Incoherent neighbor!\n");
        else {
            uint32_t v[3];
            getFaceVertices(i, v);
            if (!tetHasVertex(tet_neigh[i] >> 2, v[0])) assert(0 && "Incoherent face at neighbors!\n");
            if (!tetHasVertex(tet_neigh[i] >> 2, v[1])) assert(0 && "Incoherent face at neighbors!\n");
            if (!tetHasVertex(tet_neigh[i] >> 2, v[2])) assert(0 && "Incoherent face at neighbors!\n");
        }
    }

    // Check vt*	
    for (i = 0; i < num_vertices; i++) if (inc_tet[i]!=UINT64_MAX) {
        if (inc_tet[i] >= numTets())
            assert(0 && "Wrong vt* (out of range)!\n");
        if (isGhost(inc_tet[i]))
            assert(0 && "Wrong vt* (ghost tet)!\n");
        if (isToDeleteSmall(inc_tet[i]))
            assert(0 && "Wrong vt* (deleted tet)!\n");
        const uint32_t* tn = tet_node.data() + inc_tet[i] * 4;
        if (tn[0] != i && tn[1] != i && tn[2] != i && tn[3] != i)
            assert(0 && "Wrong vt*!\n");
    }

    // Check marks
    //for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2))
    //    if (mark_tetrahedra[i])
    //        assert(0 && "Marked tet\n");

    // Check geometry
    for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
        const uint32_t* tn = tet_node.data() + i * 4;
        if (tn[3] != INFINITE_VERTEX && vOrient3D(tn[0], tn[1], tn[2], tn[3]) <= 0) assert(0 && "Inverted/degn tet\n");
    }

    if (checkDelaunay) {
        for (size_t i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
            const uint32_t* n = tet_node.data() + (i * 4);
            if (n[3] == INFINITE_VERTEX) continue;
            for (int j = 0; j < 4; j++) {
                uint64_t oppc = tet_neigh[i * 4 + j];
                uint32_t ov = tet_node[oppc];
                if (ov != INFINITE_VERTEX && ((mark_tetrahedra[i] & 3) == (mark_tetrahedra[oppc>>2] & 3)) && vertexInTetSphere(n, ov) > 0) assert(0 && "Non delaunay\n");
            }
        }
    }

    printf("checkMesh passed\n");
}

uint32_t TetMesh::findEncroachingPoint_inexact(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    static std::vector<uint64_t> enc_queue; // Static to avoid reallocation upon each call

    // Start collecting tetrahedra incident at the endpoints
    VT(ep0, enc_queue);

    for (uint64_t j : enc_queue) mark_Tet_1(j);

    const vector3d p0 = vertices[ep0];
    const vector3d p1 = vertices[ep1];
    const double eslen = (p0 - p1).sq_length();

    vector3d ep;
    uint32_t enc_pt_i = UINT32_MAX;

    marked_vertex[ep0] = marked_vertex[ep1] = 1;

    // Collect all encroaching points while expanding around insphere vertices
    for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
        const uint64_t tet = enc_queue[ti];
        const uint64_t tb = tet << 2;

        // Check each tet vertex for 'insphereness' and keep track of the one with largest sphere
        const uint32_t* tn = tet_node.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint32_t ui = tn[i];
            if (!marked_vertex[ui]) {           
                const vector3d& pui = vertices[ui];
                if (((pui - p0).sq_length() + (pui - p1).sq_length()) <= eslen) {
                    marked_vertex[ui] = 1;
                    if (enc_pt_i == UINT32_MAX || vector3d::hasLargerSphere(p0, p1, pui, ep)) {
                        ep = pui; enc_pt_i = ui;
                        tet_e = tb;
                    }
                } 
                else marked_vertex[ui] = 2;
            }
        }

        const int nvmask[] = { (marked_vertex[tn[0]] == 1), (marked_vertex[tn[1]] == 1), (marked_vertex[tn[2]] == 1), (marked_vertex[tn[3]] == 1) };
        const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

        // Expand on adjacent tets if at least one common vertex is insphere
        const uint64_t* tg = tet_neigh.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint64_t nc = tg[i];
            const uint64_t n = nc >> 2;
            if (is_marked_Tet_1(n)==2 || tet_node[nc] == INFINITE_VERTEX) continue;
            const int totmarked = totmarkeda - nvmask[i];
            if (totmarked) {
                mark_Tet_1(n);
                enc_queue.push_back(n);
            }
        }
    }

    // Clear all marks
    marked_vertex[ep0] = marked_vertex[ep1] = 0;
    for (uint64_t j : enc_queue) {
        unmark_Tet_1(j);
        j <<= 2;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j]] = 0;
    }
    enc_queue.clear();

    return enc_pt_i;
}

uint32_t TetMesh::findEncroachingPoint_exact(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    static std::vector<uint64_t> enc_queue; // Static to avoid reallocation upon each call

    // Start collecting tetrahedra incident at the endpoints
    VT(ep0, enc_queue);

    for (uint64_t j : enc_queue) mark_Tet_1(j);

    uint32_t enc_pt_i = UINT32_MAX;

    marked_vertex[ep0] = marked_vertex[ep1] = 1;

    // Collect all encroaching points while expanding around insphere vertices
    for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
        const uint64_t tet = enc_queue[ti];
        const uint64_t tb = tet << 2;

        // Check each tet vertex for 'insphereness' and keep track of the one with largest sphere
        const uint32_t* tn = tet_node.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint32_t ui = tn[i];
            if (!marked_vertex[ui]) {
                if (genericPoint::dotProductSign3D(*vertices[ep0], *vertices[ep1], *vertices[ui]) <= 0) {
                    marked_vertex[ui] = 1;
                    if (enc_pt_i == UINT32_MAX || pointType::inGabrielSphere(*vertices[ui], *vertices[enc_pt_i], *vertices[ep0], *vertices[ep1]) < 0) {
                        enc_pt_i = ui;
                        tet_e = tb;
                    }
                }
                else marked_vertex[ui] = 2;
            }
        }

        const int nvmask[] = { (marked_vertex[tn[0]] == 1), (marked_vertex[tn[1]] == 1), (marked_vertex[tn[2]] == 1), (marked_vertex[tn[3]] == 1) };
        const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

        // Expand on adjacent tets if at least one common vertex is insphere
        const uint64_t* tg = tet_neigh.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint64_t nc = tg[i];
            const uint64_t n = nc >> 2;
            if (is_marked_Tet_1(n) == 2 || tet_node[nc] == INFINITE_VERTEX) continue;
            const int totmarked = totmarkeda - nvmask[i];
            if (totmarked) {
                mark_Tet_1(n);
                enc_queue.push_back(n);
            }
        }
    }

    // Clear all marks
    marked_vertex[ep0] = marked_vertex[ep1] = 0;
    for (uint64_t j : enc_queue) {
        unmark_Tet_1(j);
        j <<= 2;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j]] = 0;
    }
    enc_queue.clear();

    return enc_pt_i;
}


uint32_t TetMesh::findEncroachingPoint(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    uint32_t enc_pt_i = findEncroachingPoint_inexact(ep0, ep1, tet_e);
    if (enc_pt_i == UINT32_MAX) enc_pt_i = findEncroachingPoint_exact(ep0, ep1, tet_e);
    return enc_pt_i;
}
