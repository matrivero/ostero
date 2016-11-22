Point(1) = {0, 0, 0, 1.0};
Point(2) = {100, 0, 0, 1.0};
Point(3) = {120, 0, 0, 1.0};
Point(4) = {-120, 0, 0, 1.0};
Point(5) = {-100, 0, 0, 1.0};
Point(6) = {0, -100, 0, 1.0};
Point(7) = {0, -120, 0, 1.0};
Point(8) = {-100, 150, 0, 1.0};
Point(9) = {-120, 150, 0, 1.0};
Circle(1) = {2, 1, 6};
Circle(2) = {6, 1, 5};
Circle(3) = {3, 1, 7};
Circle(4) = {7, 1, 4};
Line(5) = {9, 4};
Line(6) = {9, 8};
Line(7) = {8, 5};
Line(8) = {2, 3};
Line Loop(9) = {5, -4, -3, -8, 1, 2, -7, -6};
Plane Surface(10) = {9};
Physical Line("FIXED") = {6};
Physical Line("PRESS") = {8};
Physical Surface("HOOK") = {10};
Mesh.CharacteristicLengthFactor = 3;
