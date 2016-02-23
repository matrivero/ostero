Point(1) = {0.2, 4.5, 0, 1.0};
Point(2) = {1.4, 4.5, 0, 1.0};
Point(3) = {1.4, 5.0, 0, 1.0};
Point(4) = {0.2, 5.0, 0, 1.0};
Point(5) = {0.8, 4.95, 0, 1.0};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Circle(4) = {1, 5, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
lc = 0.35;
//Mesh.CharacteristicLengthFactor = lc;
Physical Line("PRESS") = {2};
Physical Line("CONTACT") = {4};
Physical Line("LATS") = {3, 1};
Physical Surface("IDENTER") = {6};

// Unstructured Boundary Layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {4};
//Field[1].hfar = 0.6;
//Field[1].hwall_n = 0.3;
//Field[1].hwall_t = 0.3;
Field[1].hfar = 0.07;
Field[1].hwall_n = 0.02;
Field[1].hwall_t = 0.02;
Background Field = 1;
