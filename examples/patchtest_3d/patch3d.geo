N = 2 ;

Point( 1) = { 0.249 , 0.342,  0.192 ,1.0};
Point( 2) = { 0.826 , 0.288,  0.288 ,1.0};
Point( 3) = { 0.850 , 0.649,  0.263 ,1.0};
Point( 4) = { 0.273 , 0.750,  0.230 ,1.0};
Point( 5) = { 0.320 , 0.186,  0.643 ,1.0};
Point( 6) = { 0.677 , 0.305,  0.683 ,1.0};
Point( 7) = { 0.788 , 0.693,  0.644 ,1.0};
Point( 8) = { 0.165 , 0.745,  0.702 ,1.0};
Point( 9) = { 0.000 , 0.000,  0.000 ,1.0};
Point(10) = { 1.000 , 0.000,  0.000 ,1.0};
Point(11) = { 1.000 , 1.000,  0.000 ,1.0};
Point(12) = { 0.000 , 1.000,  0.000 ,1.0};
Point(13) = { 0.000 , 0.000,  1.000 ,1.0};
Point(14) = { 1.000 , 0.000,  1.000 ,1.0};
Point(15) = { 1.000 , 1.000,  1.000 ,1.0};
Point(16) = { 0.000 , 1.000,  1.000 ,1.0};

Line(1) = {8, 4};
Line(2) = {4, 3};
Line(3) = {3, 7};
Line(4) = {7, 8};
Line(5) = {8, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {5, 1};
Line(9) = {1, 4};
Line(10) = {1, 2};
Line(11) = {2, 6};
Line(12) = {2, 3};
Line Loop(13) = {5, 6, 7, 4};
Plane Surface(14) = {13};
Line Loop(15) = {11, -6, 8, 10};
Plane Surface(16) = {15};
Line Loop(17) = {12, 3, -7, -11};
Plane Surface(18) = {17};
Line Loop(19) = {2, 3, 4, 1};
Plane Surface(20) = {19};
Line Loop(21) = {5, 8, 9, -1};
Plane Surface(22) = {21};
Line Loop(23) = {10, 12, -2, -9};
Plane Surface(24) = {23};
Transfinite Line {10, 11, 12, 2, 3, 7, 4, 1, 6, 8, 9, 5} = N Using Progression 1;

Transfinite Surface {14};
Recombine Surface {14};
Transfinite Surface {20};
Recombine Surface {20};
Transfinite Surface {18};
Recombine Surface {18};
Transfinite Surface {22};
Recombine Surface {22};
Transfinite Surface {16};
Recombine Surface {16};
Transfinite Surface {24};
Recombine Surface {24};
Surface Loop(25) = {22, 14, 16, 18, 24, 20};
Volume(26) = {25};
Transfinite Volume{26};
Recombine Volume {26};

Line(27) = {16, 12};
Line(28) = {12, 11};
Line(29) = {11, 15};
Line(30) = {15, 16};
Line(31) = {12, 4};
Line(32) = {11, 3};
Line(33) = {15, 7};
Line(34) = {8, 16};
Transfinite Line {27, 34, 31, 28, 32, 29, 33, 30} = N Using Progression 1;
Line Loop(35) = {27, 31, -1, 34};
Plane Surface(36) = {35};
Transfinite Surface {36};
Recombine Surface {36};
Line Loop(37) = {28, 32, -2, -31};
Plane Surface(38) = {37};
Transfinite Surface {38};
Recombine Surface {38};
Line Loop(39) = {29, 33, -3, -32};
Plane Surface(40) = {39};
Transfinite Surface {40};
Recombine Surface {40};
Line Loop(41) = {30, -34, -4, -33};
Plane Surface(42) = {41};
Transfinite Surface {42};
Recombine Surface {42};
Line Loop(43) = {27, 28, 29, 30};
Plane Surface(44) = {43};
Transfinite Surface {44};
Recombine Surface {44};
Surface Loop(45) = {44, 36, 38, 40, 42, 20};
Volume(46) = {45};
Transfinite Volume{46};
Recombine Volume {46};

Line(47) = {16, 13};
Line(48) = {13, 14};
Line(49) = {14, 15};
Line(50) = {14, 6};
Line(51) = {5, 13};
Transfinite Line {47, 51, 48, 50, 49} = N Using Progression 1;
Line Loop(52) = {48, 50, -6, 51};
Plane Surface(53) = {52};
Transfinite Surface {53};
Recombine Surface {53};
Line Loop(54) = {49, 33, -7, -50};
Plane Surface(55) = {54};
Transfinite Surface {55};
Recombine Surface {55};
Line Loop(56) = {5, 51, -47, -34};
Plane Surface(57) = {56};
Transfinite Surface {57};
Recombine Surface {57};
Line Loop(58) = {48, 49, 30, 47};
Plane Surface(59) = {58};
Transfinite Surface {59};
Recombine Surface {59};
Surface Loop(60) = {59, 53, 55, 57, 42, 14};
Volume(61) = {60};
Transfinite Volume{61};
Recombine Volume {61};

Line(62) = {14, 10};
Line(63) = {10, 11};
Line(64) = {2, 10};
Transfinite Line {62, 64, 63} = N Using Progression 1;
Line Loop(65) = {11, -50, 62, -64};
Plane Surface(66) = {65};
Transfinite Surface {66};
Recombine Surface {66};
Line Loop(67) = {63, 32, -12, 64};
Plane Surface(68) = {67};
Transfinite Surface {68};
Recombine Surface {68};
Line Loop(69) = {29, -49, 62, 63};
Plane Surface(70) = {69};
Transfinite Surface {70};
Recombine Surface {70};
Surface Loop(71) = {70, 66, 68, 55, 18, 40};
Volume(72) = {71};
Transfinite Volume{72};
Recombine Volume {72};

Line(73) = {10, 9};
Line(74) = {9, 1};
Line(75) = {12, 9};
Transfinite Line {73, 74, 75} = N Using Progression 1;
Line Loop(76) = {73, 74, 10, 64};
Plane Surface(77) = {76};
Transfinite Surface {77};
Recombine Surface {77};
Line Loop(78) = {74, 9, -31, 75};
Plane Surface(79) = {78};
Transfinite Surface {79};
Recombine Surface {79};
Line Loop(80) = {73, -75, 28, -63};
Plane Surface(81) = {80};
Transfinite Surface {81};
Recombine Surface {81};
Surface Loop(82) = {81, 77, 79, 24, 38, 68};
Volume(83) = {82};
Transfinite Volume{83};
Recombine Volume {83};

Line(84) = {9, 13};
Transfinite Line {84} = N Using Progression 1;
Line Loop(85) = {8, -74, 84, -51};
Plane Surface(86) = {85};
Transfinite Surface {86};
Recombine Surface {86};
Line Loop(87) = {84, 48, 62, 73};
Plane Surface(88) = {87};
Transfinite Surface {88};
Recombine Surface {88};
Surface Loop(89) = {88, 86, 66, 77, 16, 53};
Volume(90) = {89};
Transfinite Volume{90};
Recombine Volume {90};

Line Loop(91) = {75, 84, -47, 27};
Plane Surface(92) = {91};
Transfinite Surface {92};
Recombine Surface {92};
Surface Loop(93) = {22, 36, 79, 86, 92, 57};
Volume(94) = {93};
Transfinite Volume{94};
Recombine Volume {94};

Physical Point("P000") = {9};
Physical Point("P100") = {10};
Physical Point("P110") = {11};
Physical Point("P010") = {12};
Physical Point("P001") = {13};
Physical Point("P101") = {14};
Physical Point("P111") = {15};
Physical Point("P011") = {16};

//Physical Surface("TOP") = {44};
//Physical Surface("RIGHT") = {70};
//Physical Surface("FRONT") = {59};
//Physical Surface("LEFT") = {92};
//Physical Surface("BACK") = {81};
//Physical Surface("BOTTOM") = {88};

Physical Volume("BULK") = {94, 90, 61, 72, 26, 46, 83};









