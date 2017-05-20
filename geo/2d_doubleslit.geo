LX = 3;
LY = 3;
WX = 0.1; // witdh of the slit
WY = 0.1; // height of the slit 
D0 = 0.1; // distence between the two slits

Point(1) = {0,0,0};
Point(2) = {LX,0,0};
Point(3) = {LX+WX,0,0};
Point(4) = {LX+WX+D0,0,0};
Point(5) = {LX+2*WX+D0,0,0};
Point(6) = {2*LX+2*WX+D0,0,0};

Point(7) = {2*LX+2*WX+D0,LY,0};
Point(8) = {LX+2*WX+D0,LY,0};
Point(9) = {LX+WX+D0,LY,0};
Point(10) = {LX+WX,LY,0};
Point(11) = {LX,LY,0};
Point(12) = {0,LY,0};

Point(13) = {0,LY+WY,0};
Point(14) = {LX,LY+WY,0};
Point(15) = {LX+WX,LY+WY,0};
Point(16) = {LX+WX+D0,LY+WY,0};
Point(17) = {LX+2*WX+D0,LY+WY,0};
Point(18) = {2*LX+2*WX+D0,LY+WY,0};

Point(19) = {2*LX+2*WX+D0,2*LY+WY,0};
Point(20) = {LX+2*WX+D0,2*LY+WY,0};
Point(21) = {LX+WX+D0,2*LY+WY,0};
Point(22) = {LX+WX,2*LY+WY,0};
Point(23) = {LX,2*LY+WY,0};
Point(24) = {0,2*LY+WY,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {5,8};
Line(8) = {4,9};
Line(9) = {3,10};
Line(10) = {2,11};
Line(11) = {1,12};
Line(12) = {11,12};
Line(13) = {10,11};
Line(14) = {9,10};
Line(15) = {8,9};
Line(16) = {7,8};
Line(17) = {8,17};
Line(18) = {9,16};
Line(19) = {10,15};
Line(20) = {11,14};
Line(21) = {13,14};
Line(22) = {14,15};
Line(23) = {15,16};
Line(24) = {16,17};
Line(25) = {17,18};
Line(26) = {18,19};
Line(27) = {17,20};
Line(28) = {16,21};
Line(29) = {15,22};
Line(30) = {14,23};
Line(31) = {13,24};
Line(32) = {24,23};
Line(33) = {23,22};
Line(34) = {22,21};
Line(35) = {21,20};
Line(36) = {20,19};
//+
Line Loop(1) = {12, -11, 1, 10};
//+
Surface(1) = {1};
//+
Line Loop(2) = {10, -13, -9, -2};
//+
Surface(2) = {2};
//+
Line Loop(3) = {9, -14, -8, -3};
//+
Surface(3) = {3};
//+
Line Loop(4) = {8, -15, -7, -4};
//+
Surface(4) = {4};
//+
Line Loop(5) = {16, -7, 5, 6};
//+
Surface(5) = {5};
//+
Line Loop(6) = {20, 22, -19, 13};
//+
Surface(6) = {6};
//+
Line Loop(7) = {18, 24, -17, 15};
//+
Surface(7) = {7};
//+
Line Loop(8) = {31, 32, -30, -21};
//+
Surface(8) = {8};
//+
Line Loop(9) = {30, 33, -29, -22};
//+
Surface(9) = {9};
//+
Line Loop(10) = {34, -28, -23, 29};
//+
Surface(10) = {10};
//+
Line Loop(11) = {28, 35, -27, -24};
//+
Surface(11) = {11};
//+
Line Loop(12) = {25, 26, -36, -27};
//+
Surface(12) = {12};
//+
Coherence;
//+
Physical Surface(1) = {8, 9, 10, 11, 12, 6, 7, 1, 2, 3, 4, 5};
//+
Physical Line(2) = {32, 33, 34, 35, 36, 26, 25, 17, 16, 6, 5, 4, 3, 2, 1, 11, 12, 20, 21, 31, 19, 23, 14, 18};


//+
Transfinite Line {32, 21, 12, 1, 5, 16, 25, 36} = 50 Using Progression 1;
//+
Transfinite Line {31, 30, 29, 28, 27, 26, 11, 10, 9, 8, 7, 6} = 50 Using Progression 1;
//+
Transfinite Line {33, 22, 13, 15, 24, 35, 2, 4} = 5 Using Progression 1;
//+
Transfinite Line {3, 14, 23, 34} = 5 Using Progression 1;
//+
Transfinite Line {20, 19, 18, 17} = 5 Using Progression 1;
//+
Transfinite Surface {1:12};
//+
Recombine Surface {1:12};
