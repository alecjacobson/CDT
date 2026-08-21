This fork allows inputting a .obj file with

```
l [v1] [v2]
```

line-segment indices in addition to faces (basic `f [v1] [v2] [v3]`). There's no
attempt to check that line segments are valid and non-redundant. They shouldn't
be degenerate, copies of each other or copies of edges in the input triangles.

-------------------------------------------------------------

# CDT - Constrained Delaunay Tetrahedrization made robust and practical
This code implements an algorithm to calculate a Constrained Delaunay Tetrahedrization (CDT) of an input PLC represented by on OFF file.
Steiner points are possibly added to make the input admit a CDT.
Details of the algorithm are described in "**Constrained Delaunay Tetrahedrization: A robust and practical approach**" by L. Diazzi, D. Panozzo, A. Vaxman and <a href="http://saturno.ge.imati.cnr.it/ima/personal-old/attene/PersonalPage/attene.html">M. Attene</a> (ACM Trans Graphics Vol 42, N. 6, Procs of SIGGRAPH Asia 2023). 
You may download a copy here: http://arxiv.org/abs/2309.09805

<p align="center"><img src="teaser_img.png"></p>

## Usage
Clone this repository with:
```
git clone https://github.com/MarcoAttene/CDT
```

Once done, you may build the executable as follows:
```
cmake -B build -S .
```

This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a cdt.sln file.
On Linux/MacOS, this will produce a Makefile. 
Use it as usual to compile cdt. Alternatively, you can use the command line:
```
cmake --build build --config Release
```

When compiled, the code generates an executable called ``cdt``.
Launch it with no command line parameters to have a list of supported options.

Example:

```
cdt input_file.off
```
creates a file called ``input_file.off.tet`` representing the constrained tetrahedrization.


We tested our code on Linux (GCC-11), MacOS (GCC-11 and CLANG) and Windows (MSVC 2022 with both CL and CLANG).

## Python bindings

Built with [nanobind](https://github.com/wjakob/nanobind) and packaged with
scikit-build-core, so a plain `pip install` is enough:

```
pip install .
```

```python
import numpy as np
import cdt

V = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
              [0,0,1],[1,0,1],[1,1,1],[0,1,1]], dtype=float)
F = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
              [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])

r = cdt.tetrahedralize(V, F)
inside = r.tets[r.labels == cdt.INNER]     # #T by 4 into r.vertices
```

`tetrahedralize(V, F, E=None, bounding_box=False, verbose=False,
inner_only=False)` returns a named tuple of `vertices`, `tets`, `labels`,
`vertex_map`, `num_steiner_vertices`, `is_polyhedron` and
`face_recovery_ok`; see `help(cdt.tetrahedralize)` for the details. `E` is an
optional `#E by 2` list of segment constraints, the in-memory equivalent of
the `l` lines this fork reads from `.obj` files.

To build in place and run the tests instead of installing:

```
cmake -B build-python -S . -DCDT_BUILD_PYTHON=ON -DCDT_BUILD_EXECUTABLE=OFF
cmake --build build-python
(cd python && python -m pytest tests)
```

## MATLAB bindings

Built the way [gptoolbox](https://github.com/alecjacobson/gptoolbox) builds
its own mex files: the result lands in `mex/` and adding that directory to the
MATLAB path is all that is needed.

```
cmake -B build-matlab -S . -DCDT_BUILD_MATLAB=ON -DCDT_BUILD_EXECUTABLE=OFF
cmake --build build-matlab
```

Point CMake at a specific installation with `-DMatlab_ROOT_DIR=...` if it
picks the wrong one.

```matlab
addpath('/path/to/CDT/mex');
[V,F] = subdivided_sphere(2);
[TV,TT,TL] = cdt(V,F);
tetramesh(TT(TL==2,:),TV);
```

See `help cdt` for the full signature and `test_cdt` for the test suite.

## Bindings vs. the command line tool

Both bindings go through `cdt::tetrahedralize` in `src/cdt.h`, which returns
the tetrahedrization as arrays instead of writing a file. Two things differ
from what the `cdt` executable writes to a `.tet` file:

- Tets are reordered so that `det[b-a, c-a, d-a] > 0`, i.e. positively
  oriented, which is what the usual signed volume formulas expect.
- Both bindings report a map from input vertices to output vertices, since
  the algorithm deduplicates and spatially reorders the input.

Vertex positions are otherwise identical, including the rounding of Steiner
points to the nearest representable double.

### Reproducibility

Segment recovery processes the missing edges in a shuffled order, and that
order decides where the Steiner points land. The generator behind the shuffle
is never reseeded on its own, which is invisible to the executable (one run
per process) but not to the bindings: without a reset, the second call in a
session would continue the first call's sequence and produce a different,
equally valid tetrahedrization.

Both bindings therefore reset the generator on entry, so the same input gives
the same output every time, matching what the executable produces. Pass
`seed=` (python) or `'Seed'` (matlab) to explore other tetrahedrizations of
the same input.

## License
This program is distributed under the terms of either the GNU GPL or the GNU LGPL license.
The code can be compiled in two ways, depending on how CMake is invoked.
If you build using ``CMake -DLGPL=ON ..``, you may choose between GPL and LGPL at your option.
If you build using ``CMake -DLGPL=OFF ..`` or just ``CMake ..``, the code makes use of modified 
parts of a third-party code which requires you to accept the terms of the GPL license.
See ``src/delaunay.h`` for details.

In either case, the program is distributed in the hope that it will be      
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.