% CDT Constrained Delaunay tetrahedrization of a piecewise linear complex,
% adding Steiner points as needed.
%
% [TV,TT] = cdt(V,F)
% [TV,TT,TL,I,S] = cdt(V,F,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into rows of V
%   Optional:
%     'Edges'  followed by a #E by 2 list of segment constraints, indices into
%       rows of V. These need not be edges of F, but must not be degenerate or
%       duplicate each other or the edges of F.
%     'BoundingBox'  followed by whether to add eight vertices enclosing the
%       input, so the output covers a box rather than the convex hull of the
%       input {false}
%     'InnerOnly'  followed by whether to drop the tets outside the input
%       polyhedron. TV is unaffected, so it may then contain unreferenced
%       vertices {false}
%     'Verbose'  followed by whether to report progress {false}
% Outputs:
%   TV  #TV by 3 list of output vertex positions. The first
%     #TV-S.num_steiner_vertices rows are the input vertices.
%   TT  #TT by 4 list of tet indices into rows of TV, positively oriented so
%     that volume(TV,TT) is positive. Tets inside the input come first.
%   TL  #TT list of labels, 1 for a tet outside the input polyhedron and 2 for
%     a tet inside it
%   I  #V list of indices into rows of TV, so that TV(I,:) == V. The input is
%     deduplicated and spatially reordered, so this is not 1:#V.
%   S  struct with fields:
%     num_steiner_vertices  number of vertices the algorithm had to add
%     is_polyhedron  false if some input edge has an odd number of incident
%       faces, in which case inside/outside is meaningless and every tet is
%       labelled inside
%     face_recovery_ok  false if face recovery fell back on the slower method
%
% Validity of the input complex is assumed, not verified. Steiner points are
% rounded to the nearest representable double on output, so an individual tet
% may come out degenerate or inverted even though the exact tetrahedrization
% is neither.
%
% Example:
%   [V,F] = subdivided_sphere(2);
%   [TV,TT,TL] = cdt(V,F);
%   tetramesh(TT(TL==2,:),TV);
%
% Implements "Constrained Delaunay Tetrahedrization: A robust and practical
% approach", Diazzi, Panozzo, Vaxman and Attene, SIGGRAPH Asia 2023.
%
