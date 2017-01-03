length  = 1;
Slength =  length/2;
chan    =  length/10;
Bdiameter= length/5;
Sdiameter= Bdiameter*0.5;

ND        = 11;
NLMEDIUM  = 21;
NLBORDERS = 9;
NLCHAN    = 5;

Point(1) = {0        , 0, 0, 1.0};
Point(2) = {Bdiameter, 0, 0, 1.0};
Point(3) = {Bdiameter, length, 0, 1.0};
Point(4) = {0        , length, 0, 1.0};

Point(5) = {Bdiameter/2-Sdiameter/2, length/2 -Slength/2, 0, 1.0};
Point(6) = {Bdiameter/2+Sdiameter/2, length/2 -Slength/2, 0, 1.0};
Point(7) = {Bdiameter/2+Sdiameter/2, length/2 +Slength/2, 0, 1.0};
Point(8) = {Bdiameter/2-Sdiameter/2, length/2 +Slength/2, 0, 1.0};

Point(9)  = {0        , length/2 -Slength/2 - chan, 0, 1.0};
Point(10) = {Bdiameter, length/2 -Slength/2 - chan, 0, 1.0};
Point(11) = {Bdiameter, length/2 +Slength/2 + chan, 0, 1.0};
Point(12) = {0        , length/2 +Slength/2 + chan, 0, 1.0};


Line(1) = {4, 3};
Line(2) = {3, 11};
Line(3) = {11, 7};
Line(4) = {7, 6};
Line(5) = {6, 10};
Line(6) = {10, 2};
Line(7) = {1, 2};
Line(8) = {1, 9};
Line(9) = {9, 5};
Line(10) = {5, 8};
Line(11) = {8, 12};
Line(12) = {12, 4};
Line(13) = {8, 7};
Line(14) = {5, 6};
Line(15) = {9, 10};
Line(16) = {12, 11};

Transfinite Line {1, 16, 13, 14, 15, 7} = ND Using Progression 1;
Transfinite Line {12, 2, 8, 6} = NLBORDERS Using Progression 1;
Transfinite Line {10, 4} = NLMEDIUM Using Progression 1;
Transfinite Line {11, 3, 9, 5} = NLCHAN Using Progression 1;

Line Loop(17) = {12, 1, 2, -16};
Plane Surface(18) = {17};
Line Loop(19) = {11, 16, 3, -13};
Plane Surface(20) = {19};
Line Loop(21) = {10, 13, 4, -14};
Plane Surface(22) = {21};
Line Loop(23) = {9, 14, 5, -15};
Plane Surface(24) = {23};
Line Loop(25) = {8, 15, 6, -7};
Plane Surface(26) = {25};

Transfinite Surface {18};
Transfinite Surface {20};
Transfinite Surface {22};
Transfinite Surface {24};
Transfinite Surface {26};

Recombine Surface {18};
Recombine Surface {20};
Recombine Surface {22};
Recombine Surface {24};
Recombine Surface {26};

Physical Surface("BULK") = {18, 20, 22, 24, 26}; 
Physical Line("TOP") = {1};
Physical Line("BOTTOM") = {7};
Physical Line("LATERALS") =  {8, 9, 10, 11, 12, 2, 3, 4, 5, 6};
