Point(1) = {-1, -1.4-0.05, 0, 1.0};
Point(2) = {1, -1.4-0.05, 0, 1.0};
Point(3) = {1, -1-0.05, 0, 1.0};
Point(4) = {0, -1-0.05, 0, 0.05};
Point(5) = {-1, -1-0.05, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line Loop(5) = {3, 4, 5, 1, 2};
Plane Surface(6) = {5};
Physical Line("FIXED") = {1};
Physical Line("CONTACT") = {3, 4};
Physical Line("LATS") = {5, 2};
Physical Surface("BLOCK") = {6};
lc = 0.1;
Mesh.CharacteristicLengthFactor = lc;

//// Unstructured Boundary Layer
//Field[1] = BoundaryLayer;
//Field[1].EdgesList = {3, 4};
//Field[1].hfar = 0.07;
//Field[1].hwall_n = 0.02;
//Field[1].hwall_t = 0.02;
//Background Field = 1;

