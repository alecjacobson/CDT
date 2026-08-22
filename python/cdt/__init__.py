"""Constrained Delaunay tetrahedrization made robust and practical.

Python bindings for https://github.com/MarcoAttene/CDT, the algorithm of
Diazzi, Panozzo, Vaxman and Attene, "Constrained Delaunay Tetrahedrization: A
robust and practical approach" (SIGGRAPH Asia 2023).
"""

from typing import NamedTuple, Optional

import numpy as np

from . import _cdt

#: Label of a tetrahedron outside the input polyhedron.
OUTER = _cdt.OUTER
#: Label of a tetrahedron inside the input polyhedron.
INNER = _cdt.INNER

__all__ = ["tetrahedralize", "Tetrahedralization", "OUTER", "INNER"]


class Tetrahedralization(NamedTuple):
    vertices: np.ndarray
    """``#V by 3`` positions. The first ``#V - num_steiner_vertices`` rows are
    the (deduplicated, reordered) input vertices."""

    tets: np.ndarray
    """``#T by 4`` indices into ``vertices``. Tets labelled ``INNER`` come
    first."""

    labels: np.ndarray
    """``#T`` array of ``INNER`` / ``OUTER``."""

    vertex_map: np.ndarray
    """``#Vin`` array; ``vertex_map[i]`` is the row of ``vertices`` holding the
    i-th input vertex. The input is deduplicated and spatially reordered, so
    this is not the identity."""

    num_steiner_vertices: int
    """Vertices the algorithm had to add for the CDT to exist."""

    is_polyhedron: bool
    """False if some input edge has an odd number of incident faces, in which
    case inside/outside is meaningless and every tet is labelled ``INNER``."""

    face_recovery_ok: bool
    """False if face recovery had to fall back on the slower method."""


def tetrahedralize(
    V,
    F,
    E: Optional[np.ndarray] = None,
    bounding_box: bool = False,
    verbose: bool = False,
    inner_only: bool = False,
    seed: int = 1,
) -> Tetrahedralization:
    """Constrained Delaunay tetrahedrization of a piecewise linear complex.

    The PLC is assumed valid; this is not verified. Steiner points are added
    as needed, so the output vertices are a superset of the input ones.

    Parameters
    ----------
    V : (#V, 3) array of float
        Vertex positions.
    F : (#F, 3) array of int
        Triangle constraints, indices into ``V``.
    E : (#E, 2) array of int, optional
        Segment constraints, indices into ``V``. These need not be edges of
        ``F``, but must not be degenerate or duplicate each other or the
        triangle edges. Note that a segment has no incident faces, so passing
        any makes the PLC non-polyhedral: ``is_polyhedron`` comes back False
        and every tet is labelled ``INNER``.
    bounding_box : bool
        Add eight vertices enclosing the input, so that the result covers a
        box rather than the convex hull of the input.
    verbose : bool
        Report progress on stdout.
    inner_only : bool
        Drop the tets outside the input polyhedron. Unreferenced vertices are
        kept, so ``vertices`` is unaffected.
    seed : int
        Segment recovery processes the missing edges in a shuffled order,
        which decides where the Steiner points land. Runs with the same input
        and seed produce the same output; different seeds produce different,
        equally valid tetrahedrizations. The default reproduces the command
        line tool.

    Returns
    -------
    Tetrahedralization

    Notes
    -----
    Steiner points are rounded to the nearest representable double on output,
    so an individual tet may come out degenerate or inverted even though the
    exact tetrahedrization is neither.
    """
    V = np.ascontiguousarray(V, dtype=np.float64)
    F = np.ascontiguousarray(F, dtype=np.int64)
    if E is not None:
        E = np.ascontiguousarray(E, dtype=np.int64)
        if E.size == 0:
            E = None
    return Tetrahedralization(
        *_cdt.tetrahedralize(V, F, E, bounding_box, verbose, inner_only, seed)
    )
