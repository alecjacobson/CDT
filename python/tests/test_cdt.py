import numpy as np
import pytest

import cdt


def cube():
    V = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=np.float64)
    F = np.array(
        [[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
         [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
         [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], dtype=np.int64)
    return V, F


def tetrahedron_with_a_notch():
    """A shape whose faces cannot be recovered without Steiner points."""
    V, F = cube()
    # Push one corner inwards to make a reflex, non convex solid.
    V = V.copy()
    V[6] = [0.4, 0.4, 0.4]
    return V, F


def volume(V, T):
    a, b, c, d = (V[T[:, k]] for k in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def test_cube():
    V, F = cube()
    r = cdt.tetrahedralize(V, F)

    assert r.vertices.shape == (8, 3)
    assert r.tets.shape[1] == 4
    assert r.labels.shape == (r.tets.shape[0],)
    assert r.num_steiner_vertices == 0
    assert r.is_polyhedron
    assert r.face_recovery_ok
    # A cube's convex hull is the cube, so nothing lands outside.
    assert np.all(r.labels == cdt.INNER)
    assert np.isclose(volume(r.vertices, r.tets).sum(), 1.0)


def test_tets_are_positively_oriented():
    V, F = cube()
    r = cdt.tetrahedralize(V, F)
    assert np.all(volume(r.vertices, r.tets) > 0)


def test_vertex_map_recovers_the_input():
    V, F = cube()
    r = cdt.tetrahedralize(V, F)
    # Vertices are deduplicated and spatially reordered, so this is a real
    # permutation and not the identity.
    assert not np.array_equal(r.vertex_map, np.arange(len(V)))
    assert np.array_equal(r.vertices[r.vertex_map], V)


def test_duplicate_input_vertices_are_merged():
    V, F = cube()
    Vd = np.vstack([V, V[0]])
    r = cdt.tetrahedralize(Vd, F)
    assert r.vertices.shape == (8, 3)
    assert r.vertex_map[-1] == r.vertex_map[0]


def test_bounding_box_adds_eight_vertices_and_outer_tets():
    V, F = cube()
    plain = cdt.tetrahedralize(V, F)
    boxed = cdt.tetrahedralize(V, F, bounding_box=True)
    assert boxed.vertices.shape[0] == plain.vertices.shape[0] + 8
    assert np.any(boxed.labels == cdt.OUTER)
    inner = boxed.tets[boxed.labels == cdt.INNER]
    assert np.isclose(volume(boxed.vertices, inner).sum(), 1.0)


def test_inner_only_drops_the_outer_tets():
    V, F = cube()
    full = cdt.tetrahedralize(V, F, bounding_box=True)
    inner = cdt.tetrahedralize(V, F, bounding_box=True, inner_only=True)
    assert np.all(inner.labels == cdt.INNER)
    assert inner.tets.shape[0] == int((full.labels == cdt.INNER).sum())
    # Vertices are untouched by the filtering.
    assert inner.vertices.shape == full.vertices.shape


def test_segment_constraint_is_recovered():
    # A diagonal of the cube is not an edge of any input triangle, so it only
    # survives as a tet edge if it is passed as a constraint.
    V, F = cube()
    E = np.array([[0, 6]], dtype=np.int64)
    r = cdt.tetrahedralize(V, F, E=E)

    ends = r.vertex_map[E[0]]
    edges = np.vstack([r.tets[:, [i, j]] for i in range(4) for j in range(i + 1, 4)])
    edges = np.sort(edges, axis=1)
    assert np.any(np.all(edges == np.sort(ends), axis=1))


def test_out_of_range_index_raises():
    V, F = cube()
    F = F.copy()
    F[0, 0] = 99
    with pytest.raises(RuntimeError, match="out of range"):
        cdt.tetrahedralize(V, F)


def test_negative_index_raises():
    V, F = cube()
    F = F.copy()
    F[0, 0] = -1
    with pytest.raises(ValueError):
        cdt.tetrahedralize(V, F)


def test_too_few_vertices_raises():
    with pytest.raises(RuntimeError, match="at least 4"):
        cdt.tetrahedralize(np.zeros((3, 3)), np.array([[0, 1, 2]]))


def test_degenerate_triangle_raises_instead_of_exiting():
    # Exercises the nfg ip_error() path: upstream it calls exit(0), which
    # would take the interpreter down with it.
    V, F = cube()
    F = np.vstack([F, [[0, 1, 1]]])
    with pytest.raises(RuntimeError, match="degenerate"):
        cdt.tetrahedralize(V, F)


def test_repeated_calls_are_reproducible():
    # The shuffle in segment recovery uses a generator that is never reseeded
    # on its own, so without an explicit reset the second call in a process
    # continues the first call's sequence and lands its Steiner points
    # somewhere else.
    V, F = cube()
    # A cube needs no Steiner points; a shape that does exercises the shuffle.
    V, F = tetrahedron_with_a_notch()
    first = cdt.tetrahedralize(V, F, bounding_box=True)
    for _ in range(3):
        again = cdt.tetrahedralize(V, F, bounding_box=True)
        assert np.array_equal(again.vertices, first.vertices)
        assert np.array_equal(again.tets, first.tets)
        assert np.array_equal(again.labels, first.labels)


def test_seed_changes_nothing_about_validity():
    V, F = tetrahedron_with_a_notch()
    a = cdt.tetrahedralize(V, F, bounding_box=True, seed=1)
    b = cdt.tetrahedralize(V, F, bounding_box=True, seed=12345)
    # Both are valid tetrahedrizations of the same input...
    for r in (a, b):
        assert np.all(volume(r.vertices, r.tets) > 0)
        assert np.array_equal(r.vertices[r.vertex_map], V)
    # ...and the seed is actually plumbed through, so re-running it repeats.
    assert np.array_equal(
        cdt.tetrahedralize(V, F, bounding_box=True, seed=12345).tets, b.tets)


def test_accepts_lists_and_other_dtypes():
    V, F = cube()
    r = cdt.tetrahedralize(V.tolist(), F.astype(np.int32))
    assert r.tets.shape[0] > 0
