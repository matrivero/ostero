N = 2 ;

Point(1) = { 0.000 , 0.000,  0.000 ,1.0};
Point(2) = { 1.000 , 0.000,  0.000 ,1.0};
Point(3) = { 1.000 , 1.000,  0.000 ,1.0};
Point(4) = { 0.000 , 1.000,  0.000 ,1.0};
Point(5) = { 0.000 , 0.000,  1.000 ,1.0};
Point(6) = { 1.000 , 0.000,  1.000 ,1.0};
Point(7) = { 1.000 , 1.000,  1.000 ,1.0};
Point(8) = { 0.000 , 1.000,  1.000 ,1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = N Using Progression 1;

Physical Point("P000") = {1};
Physical Point("P100") = {2};
Physical Point("P110") = {3};
Physical Point("P010") = {4};
Physical Point("P001") = {5};
Physical Point("P101") = {6};
Physical Point("P111") = {7};
Physical Point("P011") = {8};

Line Loop(13) = {8, 5, 6, 7};
Plane Surface(14) = {13};
Line Loop(15) = {6, -11, -2, 10};
Plane Surface(16) = {15};
Line Loop(17) = {2, 3, 4, 1};
Plane Surface(18) = {17};
Line Loop(19) = {4, 9, -8, -12};
Plane Surface(20) = {19};
Line Loop(21) = {12, -7, -11, 3};
Plane Surface(22) = {21};
Line Loop(23) = {5, -10, -1, 9};
Plane Surface(24) = {23};

Transfinite Surface {14,16,18,20,22,24};

Recombine Surface {14, 16, 24, 18, 20, 22};
Surface Loop(25) = {14, 20, 18, 16, 22, 24};
Volume(26) = {25};
Recombine Volume {26};
Physical Volume("BULK") = {26};

