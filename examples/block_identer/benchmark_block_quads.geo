Point(1) = {0, 3, 0, 1.0};
Point(2) = {1.6, 3, 0, 1.0};
Point(3) = {1.6, 4, 0, 1.0};
Point(4) = {0, 4, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Line{1,3} = 40;
Transfinite Line{2,4} = 40;
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Transfinite Surface{6};
Recombine Surface {6};

Physical Line("FIXED") = {1};
Physical Line("CONTACT") = {3};
Physical Line("LATS") = {2, 4};
Physical Surface("BLOCK") = {6};
