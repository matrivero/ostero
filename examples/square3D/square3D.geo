Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6};
}
Surface Loop(29) = {19, 6, 15, 28, 23, 27};
Volume(30) = {29};
Mesh.CharacteristicLengthFactor = 0.1;
Physical Surface("TOP") = {28};
Physical Surface("BOTTOM") = {6};
Physical Volume("BULK") = {1};
