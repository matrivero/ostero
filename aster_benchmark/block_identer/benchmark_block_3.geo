Point(1) = {0, 3.5, 0, 1.0};
Point(2) = {1.6, 3.5, 0, 1.0};
Point(3) = {1.6, 4, 0, 1.0};
Point(4) = {0, 4, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
//lc = 0.35;
//Mesh.CharacteristicLengthFactor = lc;
Physical Line("FIXED") = {1};
Physical Line("CONTACTB") = {3};
Physical Line("LATSB") = {4, 2};
Physical Surface("BLOCK") = {6};

// Unstructured Boundary Layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {3};
//Field[1].hfar = 0.6;
//Field[1].hwall_n = 0.3;
//Field[1].hwall_t = 0.3;
Field[1].hfar = 0.07;
Field[1].hwall_n = 0.02;
Field[1].hwall_t = 0.02;
Background Field = 1;

