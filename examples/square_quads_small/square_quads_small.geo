Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Line{1,3} = 4;
Transfinite Line{2,4} = 4;
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Transfinite Surface{6};
Recombine Surface {6};

Physical Line("TOP") = {3};
Physical Line("BOTTOM") = {1};
Physical Line("LATERALS") = {2, 4};
Physical Surface("BULK") = {6};
