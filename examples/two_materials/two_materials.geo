Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, 2, 0, 1.0};
Point(6) = {1, 2, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 5};
Line(6) = {5 ,6};
Line(7) = {6, 3};
//Transfinite Line{1,3} = 20;
//Transfinite Line{2,4} = 20;
Line Loop(5) = {1, 2, 3, 4};
Line Loop(6) = {3, 7, 6, 5};
Plane Surface(7) = {5};
Plane Surface(8) = {6};
//Transfinite Surface{7};
//Recombine Surface {7};

Mesh.CharacteristicLengthFactor = 0.07;

Physical Line("TOP") = {6};
Physical Line("BOTTOM") = {1};
Physical Surface("TRISBOTTOM") = {7};
Physical Surface("TRISTOP") = {8};
