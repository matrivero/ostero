length  = 1;
Bdiameter= length/3;
Sdiameter= Bdiameter*0.6;
Slength =  length/3;
ND=10;
NL=10;

Point(1) = {0        , 0, 0, 1.0};
Point(2) = {Bdiameter, 0, 0, 1.0};
Point(3) = {Bdiameter, length, 0, 1.0};
Point(4) = {0        , length, 0, 1.0};

Point(5) = {Bdiameter/2-Sdiameter/2, length/2 -Slength/2, 0, 1.0};
Point(6) = {Bdiameter/2+Sdiameter/2, length/2 -Slength/2, 0, 1.0};
Point(7) = {Bdiameter/2+Sdiameter/2, length/2 +Slength/2, 0, 1.0};
Point(8) = {Bdiameter/2-Sdiameter/2, length/2 +Slength/2, 0, 1.0};

Point(9)  = {Bdiameter/2-Sdiameter/2, 0, 0, 1.0};
Point(10) = {Bdiameter/2+Sdiameter/2, 0, 0, 1.0};
Point(11) = {Bdiameter/2+Sdiameter/2, length, 0, 1.0};
Point(12) = {Bdiameter/2-Sdiameter/2, length, 0, 1.0};

Point(13) = {0, length/2 -Slength/2, 0, 1.0};
Point(14) = {Bdiameter, length/2 -Slength/2, 0, 1.0};
Point(15) = {Bdiameter, length/2 +Slength/2, 0, 1.0};
Point(16) = {0, length/2 +Slength/2, 0, 1.0};

Line(1) = {4, 16};
Line(2) = {16, 8};
Line(3) = {8, 5};
Line(4) = {5, 13};
Line(5) = {13, 1};
Line(6) = {1, 9};
Line(7) = {9, 10};
Line(8) = {10, 2};
Line(9) = {2, 14};
Line(10) = {14, 6};
Line(11) = {6, 7};
Line(12) = {7, 15};
Line(13) = {15, 3};
Line(14) = {3, 11};
Line(15) = {11, 12};
Line(16) = {12, 4};
Line(17) = {12, 8};
Line(18) = {11, 7};
Line(19) = {8, 5};
Line(20) = {5, 9};
Line(21) = {6, 10};
Line(22) = {8, 7};
Line(23) = {5, 6};

Transfinite Line {1, 17, 18, 13, 3, 11, 5, 20, 21, 9} = NL Using Progression 1;
Transfinite Line {16, 2, 4, 6, 8, 10, 12, 14} = ND/2 Using Progression 1;
Transfinite Line {7, 23, 22, 15} = ND Using Progression 1;
Line Loop(24) = {1, 2, -17, 16};
Plane Surface(25) = {24};
Line Loop(26) = {18, 12, 13, 14};
Plane Surface(27) = {26};
Line Loop(28) = {5, 6, -20, 4};
Plane Surface(29) = {28};
Line Loop(30) = {21, 8, 9, 10};
Plane Surface(31) = {30};
Line Loop(32) = {15, 17, 22, -18};
Plane Surface(33) = {32};
Line Loop(34) = {22, -11, -23, -3};
Plane Surface(35) = {34};
Line Loop(36) = {23, 21, -7, -20};
Plane Surface(37) = {36};

Transfinite Surface {25};
Transfinite Surface {33};
Transfinite Surface {27};
Transfinite Surface {35};
Transfinite Surface {29};
Transfinite Surface {37};
Transfinite Surface {31};

Recombine Surface {33, 27, 25, 35, 37, 31, 29};
Physical Surface("BULK") = {25, 33, 27, 35, 29, 37, 31};
Physical Line("TOP") = {16, 15, 14};
Physical Line("BOTTOM") = {6, 7, 8};
Physical Line("LATERALS") = {1, 2, 3, 4, 5, 9, 10, 11, 12, 13};
