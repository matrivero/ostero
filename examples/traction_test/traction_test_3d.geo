length  = 1;
Bdiameter= length/3;
Sdiameter= Bdiameter*0.6;
Slength =  length/3;
ND=10;
NLI=5;
NLB=5;
NAI=5;
NAE=5;

Point(1) = {0        , 0, 0, 1.0};
Point(2) = {Sdiameter, 0, 0, 1.0};
Point(3) = {Bdiameter, 0, 0, 1.0};
Point(4) = {Bdiameter, 0, length-2*Slength, 1.0};
Point(5) = {Sdiameter, 0, length-2*Slength, 1.0};
Point(6) = {Sdiameter, 0, length-  Slength, 1.0};
Point(7) = {Bdiameter, 0, length-  Slength, 1.0};
Point(8) = {Bdiameter, 0, length, 1.0};
Point(9) = {0, 0, length, 1.0};
Point(10) = {0, 0, length-2*Slength, 1.0};
Point(11) = {0, 0, length-  Slength, 1.0};
Point(12) = {Sdiameter, 0, length, 1.0};

Line(1) = {9, 11};
Line(2) = {11, 10};
Line(3) = {10, 1};
Line(4) = {1, 2};
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 5};
Line(8) = {5, 6};
Line(9) = {6, 7};
Line(10) = {7, 8};
Line(11) = {8, 12};
Line(12) = {12, 9};
Line(13) = {12, 6};
Line(14) = {5, 2};
Line(15) = {5, 10};
Line(16) = {6, 11};

Transfinite Line {1, 13, 10, 3, 14, 6} = NLB Using Progression 1;
Transfinite Line {4, 15, 16, 12} = NAI Using Progression 1;
Transfinite Line {2, 8} = NLI Using Progression 1;
Transfinite Line {5, 7, 9, 11} = NAE Using Progression 1;


Line Loop(17) = {1, -16, -13, 12};
Plane Surface(18) = {17};
Line Loop(19) = {11, 13, 9, 10};
Plane Surface(20) = {19};
Line Loop(21) = {16, 2, -15, 8};
Plane Surface(22) = {21};
Line Loop(23) = {14, -4, -3, -15};
Plane Surface(24) = {23};
Line Loop(25) = {14, 5, 6, 7};
Plane Surface(26) = {25};

Transfinite Surface {18};
Transfinite Surface {22};
Transfinite Surface {20};
Transfinite Surface {24};
Transfinite Surface {26};

Recombine Surface {18, 20, 22, 24, 26};

ext1[]=Extrude { {0,0,1}, {0,0,0}, Pi/2 } {
  Surface{18,22,20,24,26};Layers{3}; Recombine;
};
ext2[]=Extrude { {0,0,1}, {0,0,0}, Pi/2 } {
  Surface{43,82,60,121,99};Layers{3}; Recombine;
};
ext3[]=Extrude { {0,0,1}, {0,0,0}, Pi/2 } {
  Surface{216,177,138,160,199};Layers{3}; Recombine;
};
ext4[]=Extrude { {0,0,1}, {0,0,0}, Pi/2 } {
  Surface{311,233,250,267,289};Layers{3}; Recombine;
};



Physical Surface("TOP") = {147, 137, 276, 266, 380, 389, 69, 42};
Physical Surface("BOTTOM") = {112, 94, 324, 344, 228, 302, 190, 211};
Physical Surface("LATERALS") = {159, 81, 59, 176, 194, 116, 401, 364, 328, 288, 249, 306};

Physical Volume("BULK") = {6, 7, 9, 8, 10, 5, 4, 3, 2, 1, 11, 12, 15, 13, 14, 20, 19, 18, 17, 16};
