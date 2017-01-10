Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Transfinite Line {4, 1, 2, 3} = 2 Using Progression 1;
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};
Extrude {0, 0, 1} {
  Surface{6}; Layers{1}; Recombine;
}
Physical Surface("TOP") = {28};
Physical Surface("BOTTOM") = {6};
Physical Surface("LATERALS") = {23, 27, 15, 19};
Physical Volume("BULK") = {1};
