function test_cdt()
  % TEST_CDT Exercise the cdt mex binding. Errors on the first failure.
  %
  % test_cdt()
  %

  [V,F] = cube();

  test_cube(V,F);
  test_orientation(V,F);
  test_vertex_map(V,F);
  test_duplicate_vertices(V,F);
  test_bounding_box(V,F);
  test_inner_only(V,F);
  test_edge_constraint(V,F);
  test_errors(V,F);
  test_integer_class_input(V,F);

  fprintf('test_cdt: all tests passed\n');
end

function [V,F] = cube()
  V = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
  F = [1 3 2;1 4 3;5 6 7;5 7 8;1 2 6;1 6 5;2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8];
end

function v = tet_volume(V,T)
  a = V(T(:,1),:); b = V(T(:,2),:); c = V(T(:,3),:); d = V(T(:,4),:);
  v = sum(cross(b-a,c-a,2).*(d-a),2)./6;
end

function test_cube(V,F)
  [TV,TT,TL,I,S] = cdt(V,F);
  assert(isequal(size(TV),[8 3]),'TV should have one row per input vertex');
  assert(size(TT,2) == 4,'TT should be #TT by 4');
  assert(isequal(size(TL),[size(TT,1) 1]),'TL should be #TT by 1');
  assert(isequal(size(I),[size(V,1) 1]),'I should be #V by 1');
  assert(S.num_steiner_vertices == 0,'a cube needs no Steiner points');
  assert(S.is_polyhedron && S.face_recovery_ok);
  % The convex hull of a cube is the cube, so nothing lands outside.
  assert(all(TL == 2),'every tet should be labelled inside');
  assert(abs(sum(tet_volume(TV,TT)) - 1) < 1e-12,'total volume should be 1');
end

function test_orientation(V,F)
  [TV,TT] = cdt(V,F);
  assert(all(tet_volume(TV,TT) > 0),'tets should be positively oriented');
end

function test_vertex_map(V,F)
  [TV,~,~,I] = cdt(V,F);
  % Vertices are deduplicated and spatially reordered, so this is a real
  % permutation and not 1:#V.
  assert(~isequal(I,(1:size(V,1))'),'I should not be the identity here');
  assert(isequal(TV(I,:),V),'TV(I,:) should recover the input vertices');
end

function test_duplicate_vertices(V,F)
  [TV,~,~,I] = cdt([V;V(1,:)],F);
  assert(size(TV,1) == 8,'duplicate input vertices should be merged');
  assert(I(end) == I(1),'both copies should map to the same output vertex');
end

function test_bounding_box(V,F)
  [TV0,~] = cdt(V,F);
  [TV,TT,TL] = cdt(V,F,'BoundingBox',true);
  assert(size(TV,1) == size(TV0,1)+8,'the box should add eight vertices');
  assert(any(TL == 1),'the box should put some tets outside');
  assert(abs(sum(tet_volume(TV,TT(TL==2,:))) - 1) < 1e-12);
end

function test_inner_only(V,F)
  [TV,~,TL] = cdt(V,F,'BoundingBox',true);
  [TVi,TTi,TLi] = cdt(V,F,'BoundingBox',true,'InnerOnly',true);
  assert(all(TLi == 2));
  assert(size(TTi,1) == sum(TL == 2),'InnerOnly should keep exactly the inside tets');
  assert(isequal(TVi,TV),'InnerOnly should leave the vertices alone');
end

function test_edge_constraint(V,F)
  % A diagonal of the cube is not an edge of any input triangle, so it only
  % survives as a tet edge if it is passed as a constraint.
  E = [1 7];
  [~,TT,~,I] = cdt(V,F,'Edges',E);
  ends = sort(reshape(I(E),1,2));
  edges = sort(reshape(TT(:,nchoosek(1:4,2)),[],2),2);
  assert(any(all(edges == ends,2)),'the constraint should appear as a tet edge');
end

function test_errors(V,F)
  bad = F; bad(1,1) = 99;
  assert_error(@() cdt(V,bad),'out of range index should error');
  bad = F; bad(1,1) = 0;
  assert_error(@() cdt(V,bad),'zero index should error');
  bad = F; bad(1,1) = 1.5;
  assert_error(@() cdt(V,bad),'non integer index should error');
  assert_error(@() cdt(V(1:3,:),[1 2 3]),'too few vertices should error');
  % Exercises the nfg ip_error() path: upstream it calls exit(0), which would
  % take matlab down with it.
  assert_error(@() cdt(V,[F;1 2 2]),'degenerate triangle should error');
  assert_error(@() cdt(V,F,'NoSuchOption',true),'unknown parameter should error');
end

function test_integer_class_input(V,F)
  [~,TT] = cdt(V,int32(F));
  assert(size(TT,1) > 0,'integer class index input should be accepted');
end

function assert_error(fh,msg)
  try
    fh();
  catch
    return;
  end
  error(msg);
end
