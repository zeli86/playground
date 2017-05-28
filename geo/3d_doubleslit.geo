LX = 2;
LY = 2;
LZ = 2;
WX = 2; // witdh of the slit
WY = 2; // height of the slit 
WZ = 2; 
D0 = 2; // distence between the two slits

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

Line(newl) = {1,2};
Line(newl) = {2,3};
Line(newl) = {3,4};
Line(newl) = {4,5};
Line(newl) = {5,6};
Line(newl) = {6,7};
Line(newl) = {5,8};
Line(newl) = {4,9};
Line(newl) = {3,10};
Line(newl) = {2,11};
Line(newl) = {1,12};
Line(newl) = {11,12};
Line(newl) = {10,11};
Line(newl) = {9,10};
Line(newl) = {8,9};
Line(newl) = {7,8};
Line(newl) = {8,17};
Line(newl) = {9,16};
Line(newl) = {10,15};
Line(newl) = {11,14};
Line(newl) = {13,14};
Line(newl) = {14,15};
Line(newl) = {15,16};
Line(newl) = {16,17};
Line(newl) = {17,18};
Line(newl) = {18,19};
Line(newl) = {17,20};
Line(newl) = {16,21};
Line(newl) = {15,22};
Line(newl) = {14,23};
Line(newl) = {13,24};
Line(newl) = {24,23};
Line(newl) = {23,22};
Line(newl) = {22,21};
Line(newl) = {21,20};
Line(newl) = {20,19};
Line(newl) = {12,13};
Line(newl) = {7,18};

///

Point(newp) = {0,0,LZ};
Point(newp) = {LX,0,LZ};
Point(newp) = {LX+WX,0,LZ};
Point(newp) = {LX+WX+D0,0,LZ};
Point(newp) = {LX+2*WX+D0,0,LZ};
Point(newp) = {2*LX+2*WX+D0,0,LZ};

Point(newp) = {2*LX+2*WX+D0,LY,LZ};
Point(newp) = {LX+2*WX+D0,LY,LZ};
Point(newp) = {LX+WX+D0,LY,LZ};
Point(newp) = {LX+WX,LY,LZ};
Point(newp) = {LX,LY,LZ};
Point(newp) = {0,LY,LZ};

Point(newp) = {0,LY+WY,LZ};
Point(newp) = {LX,LY+WY,LZ};
Point(newp) = {LX+WX,LY+WY,LZ};
Point(newp) = {LX+WX+D0,LY+WY,LZ};
Point(newp) = {LX+2*WX+D0,LY+WY,LZ};
Point(newp) = {2*LX+2*WX+D0,LY+WY,LZ};

Point(newp) = {2*LX+2*WX+D0,2*LY+WY,LZ};
Point(newp) = {LX+2*WX+D0,2*LY+WY,LZ};
Point(newp) = {LX+WX+D0,2*LY+WY,LZ};
Point(newp) = {LX+WX,2*LY+WY,LZ};
Point(newp) = {LX,2*LY+WY,LZ};
Point(newp) = {0,2*LY+WY,LZ};

Line(newl) = {1+24,2+24};
Line(newl) = {2+24,3+24};
Line(newl) = {3+24,4+24};
Line(newl) = {4+24,5+24};
Line(newl) = {5+24,6+24};
Line(newl) = {6+24,7+24};
Line(newl) = {5+24,8+24};
Line(newl) = {4+24,9+24};
Line(newl) = {3+24,10+24};
Line(newl) = {2+24,11+24};
Line(newl) = {1+24,12+24};
Line(newl) = {11+24,12+24};
Line(newl) = {10+24,11+24};
Line(newl) = {9+24,10+24};
Line(newl) = {8+24,9+24};
Line(newl) = {7+24,8+24};
Line(newl) = {8+24,17+24};
Line(newl) = {9+24,16+24};
Line(newl) = {10+24,15+24};
Line(newl) = {11+24,14+24};
Line(newl) = {13+24,14+24};
Line(newl) = {14+24,15+24};
Line(newl) = {15+24,16+24};
Line(newl) = {16+24,17+24};
Line(newl) = {17+24,18+24};
Line(newl) = {18+24,19+24};
Line(newl) = {17+24,20+24};
Line(newl) = {16+24,21+24};
Line(newl) = {15+24,22+24};
Line(newl) = {14+24,23+24};
Line(newl) = {13+24,24+24};
Line(newl) = {24+24,23+24};
Line(newl) = {23+24,22+24};
Line(newl) = {22+24,21+24};
Line(newl) = {21+24,20+24};
Line(newl) = {20+24,19+24};
Line(newl) = {12+24,13+24};
Line(newl) = {7+24,18+24};

Line(newl) = {1,1+24};
Line(newl) = {2,2+24};
Line(newl) = {3,3+24};
Line(newl) = {4,4+24};
Line(newl) = {5,5+24};
Line(newl) = {6,6+24};
Line(newl) = {7,7+24};
Line(newl) = {8,8+24};
Line(newl) = {9,9+24};
Line(newl) = {10,10+24};
Line(newl) = {11,11+24};
Line(newl) = {12,12+24};
Line(newl) = {13,13+24};
Line(newl) = {14,14+24};
Line(newl) = {15,15+24};
Line(newl) = {16,16+24};
Line(newl) = {17,17+24};
Line(newl) = {18,18+24};
Line(newl) = {19,19+24};
Line(newl) = {20,20+24};
Line(newl) = {21,21+24};
Line(newl) = {22,22+24};
Line(newl) = {23,23+24};
Line(newl) = {24,24+24};

///

Point(newp) = {0,0,LZ+WZ};
Point(newp) = {LX,0,LZ+WZ};
Point(newp) = {LX+WX,0,LZ+WZ};
Point(newp) = {LX+WX+D0,0,LZ+WZ};
Point(newp) = {LX+2*WX+D0,0,LZ+WZ};
Point(newp) = {2*LX+2*WX+D0,0,LZ+WZ};

Point(newp) = {2*LX+2*WX+D0,LY,LZ+WZ};
Point(newp) = {LX+2*WX+D0,LY,LZ+WZ};
Point(newp) = {LX+WX+D0,LY,LZ+WZ};
Point(newp) = {LX+WX,LY,LZ+WZ};
Point(newp) = {LX,LY,LZ+WZ};
Point(newp) = {0,LY,LZ+WZ};

Point(newp) = {0,LY+WY,LZ+WZ};
Point(newp) = {LX,LY+WY,LZ+WZ};
Point(newp) = {LX+WX,LY+WY,LZ+WZ};
Point(newp) = {LX+WX+D0,LY+WY,LZ+WZ};
Point(newp) = {LX+2*WX+D0,LY+WY,LZ+WZ};
Point(newp) = {2*LX+2*WX+D0,LY+WY,LZ+WZ};

Point(newp) = {2*LX+2*WX+D0,2*LY+WY,LZ+WZ};
Point(newp) = {LX+2*WX+D0,2*LY+WY,LZ+WZ};
Point(newp) = {LX+WX+D0,2*LY+WY,LZ+WZ};
Point(newp) = {LX+WX,2*LY+WY,LZ+WZ};
Point(newp) = {LX,2*LY+WY,LZ+WZ};
Point(newp) = {0,2*LY+WY,LZ+WZ};

Line(newl) = {1+48,2+48};
Line(newl) = {2+48,3+48};
Line(newl) = {3+48,4+48};
Line(newl) = {4+48,5+48};
Line(newl) = {5+48,6+48};
Line(newl) = {6+48,7+48};
Line(newl) = {5+48,8+48};
Line(newl) = {4+48,9+48};
Line(newl) = {3+48,10+48};
Line(newl) = {2+48,11+48};
Line(newl) = {1+48,12+48};
Line(newl) = {11+48,12+48};
Line(newl) = {10+48,11+48};
Line(newl) = {9+48,10+48};
Line(newl) = {8+48,9+48};
Line(newl) = {7+48,8+48};
Line(newl) = {8+48,17+48};
Line(newl) = {9+48,16+48};
Line(newl) = {10+48,15+48};
Line(newl) = {11+48,14+48};
Line(newl) = {13+48,14+48};
Line(newl) = {14+48,15+48};
Line(newl) = {15+48,16+48};
Line(newl) = {16+48,17+48};
Line(newl) = {17+48,18+48};
Line(newl) = {18+48,19+48};
Line(newl) = {17+48,20+48};
Line(newl) = {16+48,21+48};
Line(newl) = {15+48,22+48};
Line(newl) = {14+48,23+48};
Line(newl) = {13+48,24+48};
Line(newl) = {24+48,23+48};
Line(newl) = {23+48,22+48};
Line(newl) = {22+48,21+48};
Line(newl) = {21+48,20+48};
Line(newl) = {20+48,19+48};
Line(newl) = {12+48,13+48};
Line(newl) = {7+48,18+48};

Line(newl) = {1+24,1+48};
Line(newl) = {2+24,2+48};
Line(newl) = {3+24,3+48};
Line(newl) = {4+24,4+48};
Line(newl) = {5+24,5+48};
Line(newl) = {6+24,6+48};
Line(newl) = {7+24,7+48};
Line(newl) = {8+24,8+48};
Line(newl) = {9+24,9+48};
Line(newl) = {10+24,10+48};
Line(newl) = {11+24,11+48};
Line(newl) = {12+24,12+48};
Line(newl) = {13+24,13+48};
Line(newl) = {14+24,14+48};
Line(newl) = {15+24,15+48};
Line(newl) = {16+24,16+48};
Line(newl) = {17+24,17+48};
Line(newl) = {18+24,18+48};
Line(newl) = {19+24,19+48};
Line(newl) = {20+24,20+48};
Line(newl) = {21+24,21+48};
Line(newl) = {22+24,22+48};
Line(newl) = {23+24,23+48};
Line(newl) = {24+24,24+48};

///

Point(newp) = {0,0,2*LZ+WZ};
Point(newp) = {LX,0,2*LZ+WZ};
Point(newp) = {LX+WX,0,2*LZ+WZ};
Point(newp) = {LX+WX+D0,0,2*LZ+WZ};
Point(newp) = {LX+2*WX+D0,0,2*LZ+WZ};
Point(newp) = {2*LX+2*WX+D0,0,2*LZ+WZ};

Point(newp) = {2*LX+2*WX+D0,LY,2*LZ+WZ};
Point(newp) = {LX+2*WX+D0,LY,2*LZ+WZ};
Point(newp) = {LX+WX+D0,LY,2*LZ+WZ};
Point(newp) = {LX+WX,LY,2*LZ+WZ};
Point(newp) = {LX,LY,2*LZ+WZ};
Point(newp) = {0,LY,2*LZ+WZ};

Point(newp) = {0,LY+WY,2*LZ+WZ};
Point(newp) = {LX,LY+WY,2*LZ+WZ};
Point(newp) = {LX+WX,LY+WY,2*LZ+WZ};
Point(newp) = {LX+WX+D0,LY+WY,2*LZ+WZ};
Point(newp) = {LX+2*WX+D0,LY+WY,2*LZ+WZ};
Point(newp) = {2*LX+2*WX+D0,LY+WY,2*LZ+WZ};

Point(newp) = {2*LX+2*WX+D0,2*LY+WY,2*LZ+WZ};
Point(newp) = {LX+2*WX+D0,2*LY+WY,2*LZ+WZ};
Point(newp) = {LX+WX+D0,2*LY+WY,2*LZ+WZ};
Point(newp) = {LX+WX,2*LY+WY,2*LZ+WZ};
Point(newp) = {LX,2*LY+WY,2*LZ+WZ};
Point(newp) = {0,2*LY+WY,2*LZ+WZ};

Line(newl) = {1+72,2+72};
Line(newl) = {2+72,3+72};
Line(newl) = {3+72,4+72};
Line(newl) = {4+72,5+72};
Line(newl) = {5+72,6+72};
Line(newl) = {6+72,7+72};
Line(newl) = {5+72,8+72};
Line(newl) = {4+72,9+72};
Line(newl) = {3+72,10+72};
Line(newl) = {2+72,11+72};
Line(newl) = {1+72,12+72};
Line(newl) = {11+72,12+72};
Line(newl) = {10+72,11+72};
Line(newl) = {9+72,10+72};
Line(newl) = {8+72,9+72};
Line(newl) = {7+72,8+72};
Line(newl) = {8+72,17+72};
Line(newl) = {9+72,16+72};
Line(newl) = {10+72,15+72};
Line(newl) = {11+72,14+72};
Line(newl) = {13+72,14+72};
Line(newl) = {14+72,15+72};
Line(newl) = {15+72,16+72};
Line(newl) = {16+72,17+72};
Line(newl) = {17+72,18+72};
Line(newl) = {18+72,19+72};
Line(newl) = {17+72,20+72};
Line(newl) = {16+72,21+72};
Line(newl) = {15+72,22+72};
Line(newl) = {14+72,23+72};
Line(newl) = {13+72,24+72};
Line(newl) = {24+72,23+72};
Line(newl) = {23+72,22+72};
Line(newl) = {22+72,21+72};
Line(newl) = {21+72,20+72};
Line(newl) = {20+72,19+72};
Line(newl) = {12+72,13+72};
Line(newl) = {7+72,18+72};

Line(newl) = {1+48,1+72};
Line(newl) = {2+48,2+72};
Line(newl) = {3+48,3+72};
Line(newl) = {4+48,4+72};
Line(newl) = {5+48,5+72};
Line(newl) = {6+48,6+72};
Line(newl) = {7+48,7+72};
Line(newl) = {8+48,8+72};
Line(newl) = {9+48,9+72};
Line(newl) = {10+48,10+72};
Line(newl) = {11+48,11+72};
Line(newl) = {12+48,12+72};
Line(newl) = {13+48,13+72};
Line(newl) = {14+48,14+72};
Line(newl) = {15+48,15+72};
Line(newl) = {16+48,16+72};
Line(newl) = {17+48,17+72};
Line(newl) = {18+48,18+72};
Line(newl) = {19+48,19+72};
Line(newl) = {20+48,20+72};
Line(newl) = {21+48,21+72};
Line(newl) = {22+48,22+72};
Line(newl) = {23+48,23+72};
Line(newl) = {24+48,24+72};

/// 

//+
Line Loop(1) = {1, 78, -39, -77};
//+
Surface(1) = {1};
//+
Line Loop(2) = {79, -40, -78, 2};
//+
Surface(2) = {2};
//+
Line Loop(3) = {41, -80, -3, 79};
//+
Surface(3) = {3};
//+
Line Loop(4) = {4, 81, -42, -80};
//+
Surface(4) = {4};
//+
Line Loop(5) = {82, -43, -81, 5};
//+
Surface(5) = {5};
//+
Line Loop(6) = {139, 101, -140, -39};
//+
Surface(6) = {6};
//+
Line Loop(7) = {202, -163, -201, 101};
//+
Surface(7) = {7};
//+
Line Loop(8) = {40, 141, -102, -140};
//+
Surface(8) = {8};
//+
Line Loop(9) = {141, 103, -142, -41};
//+
Surface(9) = {9};
//+
Line Loop(10) = {42, 143, -104, -142};
//+
Surface(10) = {10};
//+
Line Loop(11) = {144, -105, -143, 43};
//+
Surface(11) = {11};
//+
Line Loop(12) = {105, 206, -167, -205};
//+
Surface(12) = {12};
//+
Line Loop(13) = {166, -205, -104, 204};
//+
Surface(13) = {13};
//+
Line Loop(14) = {103, 204, -165, -203};
//+
Surface(14) = {14};
//+
Line Loop(15) = {164, -203, -102, 202};
//+
Surface(15) = {15};
//+
Line Loop(16) = {12, 88, -50, -87};
//+
Surface(16) = {16};
//+
Line Loop(17) = {51, -87, -13, 86};
//+
Surface(17) = {17};
//+
Line Loop(18) = {14, 86, -52, -85};
//+
Surface(18) = {18};
//+
Line Loop(19) = {15, 85, -53, -84};
//+
Surface(19) = {19};
//+
Line Loop(20) = {16, 84, -54, -83};
//+
Surface(20) = {20};
//+
Line Loop(21) = {53, 147, -115, -146};
//+
Surface(21) = {21};
//+
Line Loop(22) = {54, 146, -116, -145};
//+
Surface(22) = {22};
//+
Line Loop(23) = {147, 114, -148, -52};
//+
Surface(23) = {23};
//+
Line Loop(24) = {51, 149, -113, -148};
//+
Surface(24) = {24};
//+
Line Loop(25) = {50, 150, -112, -149};
//+
Surface(25) = {25};
//+
Line Loop(26) = {212, -174, -211, 112};
//+
Surface(26) = {26};
//+
Line Loop(27) = {113, 211, -175, -210};
//+
Surface(27) = {27};
//+
Line Loop(28) = {114, 210, -176, -209};
//+
Surface(28) = {28};
//+
Line Loop(29) = {177, -209, -115, 208};
//+
Surface(29) = {29};
//+
Line Loop(30) = {178, -208, -116, 207};
//+
Surface(30) = {30};
//+
Line Loop(31) = {94, -63, -93, 25};
//+
Surface(31) = {31};
//+
Line Loop(32) = {62, -93, -24, 92};
//+
Surface(32) = {32};
//+
Line Loop(33) = {23, 92, -61, -91};
//+
Surface(33) = {33};
//+
Line Loop(34) = {22, 91, -60, -90};
//+
Surface(34) = {34};
//+
Line Loop(35) = {21, 90, -59, -89};
//+
Surface(35) = {35};
//+
Line Loop(36) = {152, -121, -151, 59};
//+
Surface(36) = {36};
//+
Line Loop(37) = {60, 153, -122, -152};
//+
Surface(37) = {37};
//+
Line Loop(38) = {61, 154, -123, -153};
//+
Surface(38) = {38};
//+
Line Loop(39) = {124, -155, -62, 154};
//+
Surface(39) = {39};
//+
Line Loop(40) = {125, -156, -63, 155};
//+
Surface(40) = {40};
//+
Line Loop(41) = {218, -187, -217, 125};
//+
Surface(41) = {41};
//+
Line Loop(42) = {124, 217, -186, -216};
//+
Surface(42) = {42};
//+
Line Loop(43) = {123, 216, -185, -215};
//+
Surface(43) = {43};
//+
Line Loop(44) = {122, 215, -184, -214};
//+
Surface(44) = {44};
//+
Line Loop(45) = {121, 214, -183, -213};
//+
Surface(45) = {45};
//+
Line Loop(46) = {36, 95, -74, -96};
//+
Surface(46) = {46};
//+
Line Loop(47) = {96, -73, -97, 35};
//+
Surface(47) = {47};
//+
Line Loop(48) = {34, 97, -72, -98};
//+
Surface(48) = {48};
//+
Line Loop(49) = {33, 98, -71, -99};
//+
Surface(49) = {49};
//+
Line Loop(50) = {32, 99, -70, -100};
//+
Surface(50) = {50};
//+
Line Loop(51) = {161, -132, -162, 70};
//+
Surface(51) = {51};
//+
Line Loop(52) = {71, 160, -133, -161};
//+
Surface(52) = {52};
//+
Line Loop(53) = {134, -159, -72, 160};
//+
Surface(53) = {53};
//+
Line Loop(54) = {135, -158, -73, 159};
//+
Surface(54) = {54};
//+
Line Loop(55) = {157, -136, -158, 74};
//+
Surface(55) = {55};
//+
Line Loop(56) = {136, 219, -198, -220};
//+
Surface(56) = {56};
//+
Line Loop(57) = {135, 220, -197, -221};
//+
Surface(57) = {57};
//+
Line Loop(58) = {196, -221, -134, 222};
//+
Surface(58) = {58};
//+
Line Loop(59) = {195, -222, -133, 223};
//+
Surface(59) = {59};
//+
Line Loop(60) = {132, 223, -194, -224};
//+
Surface(60) = {60};
//+
Line Loop(61) = {219, -188, -218, 126};
//+
Surface(61) = {61};
//+
Line Loop(62) = {126, -157, -64, 156};
//+
Surface(62) = {62};
//+
Line Loop(63) = {26, 95, -64, -94};
//+
Surface(63) = {63};
//+
Line Loop(64) = {38, 94, -76, -83};
//+
Surface(64) = {64};
//+
Line Loop(65) = {76, 156, -138, -145};
//+
Surface(65) = {65};
//+
Line Loop(66) = {200, -218, -138, 207};
//+
Surface(66) = {66};
//+
Line Loop(67) = {207, -168, -206, 106};
//+
Surface(67) = {67};
//+
Line Loop(68) = {106, -145, -44, 144};
//+
Surface(68) = {68};
//+
Line Loop(69) = {44, -83, -6, 82};
//+
Surface(69) = {69};
//+
Line Loop(70) = {189, -220, -127, 217};
//+
Surface(70) = {70};
//+
Line Loop(71) = {127, -158, -65, 155};
//+
Surface(71) = {71};
//+
Line Loop(72) = {65, -96, -27, 93};
//+
Surface(72) = {72};
//+
Line Loop(73) = {179, -217, -117, 208};
//+
Surface(73) = {73};
//+
Line Loop(74) = {55, 155, -117, -146};
//+
Surface(74) = {74};
//+
Line Loop(75) = {55, -93, -17, 84};
//+
Surface(75) = {75};
//+
Line Loop(76) = {169, -208, -107, 205};
//+
Surface(76) = {76};
//+
Line Loop(77) = {107, -146, -45, 143};
//+
Surface(77) = {77};
//+
Line Loop(78) = {45, -84, -7, 81};
//+
Surface(78) = {78};
//+
Line Loop(79) = {190, -221, -128, 216};
//+
Surface(79) = {79};
//+
Line Loop(80) = {128, -159, -66, 154};
//+
Surface(80) = {80};
//+
Line Loop(81) = {66, -97, -28, 92};
//+
Surface(81) = {81};
//+
Line Loop(82) = {180, -216, -118, 209};
//+
Surface(82) = {82};
//+
Line Loop(83) = {118, -154, -56, 147};
//+
Surface(83) = {83};
//+
Line Loop(84) = {56, -92, -18, 85};
//+
Surface(84) = {84};
//+
Line Loop(85) = {170, -209, -108, 204};
//+
Surface(85) = {85};
//+
Line Loop(86) = {108, -147, -46, 142};
//+
Surface(86) = {86};
//+
Line Loop(87) = {85, -46, -80, 8};
//+
Surface(87) = {87};
//+
Line Loop(88) = {191, -222, -129, 215};
//+
Surface(88) = {88};
//+
Line Loop(89) = {129, -160, -67, 153};
//+
Surface(89) = {89};
//+
Line Loop(90) = {67, -98, -29, 91};
//+
Surface(90) = {90};
//+
Line Loop(91) = {181, -215, -119, 210};
//+
Surface(91) = {91};
//+
Line Loop(92) = {119, -153, -57, 148};
//+
Surface(92) = {92};
//+
Line Loop(93) = {57, -91, -19, 86};
//+
Surface(93) = {93};
//+
Line Loop(94) = {171, -210, -109, 203};
//+
Surface(94) = {94};
//+
Line Loop(95) = {109, -148, -47, 141};
//+
Surface(95) = {95};
//+
Line Loop(96) = {9, 86, -47, -79};
//+
Surface(96) = {96};
//+
Line Loop(97) = {192, -223, -130, 214};
//+
Surface(97) = {97};
//+
Line Loop(98) = {130, -161, -68, 152};
//+
Surface(98) = {98};
//+
Line Loop(99) = {68, -99, -30, 90};
//+
Surface(99) = {99};
//+
Line Loop(100) = {182, -214, -120, 211};
//+
Surface(100) = {100};
//+
Line Loop(101) = {120, -152, -58, 149};
//+
Surface(101) = {101};
//+
Line Loop(102) = {58, -90, -20, 87};
//+
Surface(102) = {102};
//+
Line Loop(103) = {211, -172, -202, 110};
//+
Surface(103) = {103};
//+
Line Loop(104) = {110, -149, -48, 140};
//+
Surface(104) = {104};
//+
Line Loop(105) = {48, -87, -10, 78};
//+
Surface(105) = {105};
//+
Line Loop(106) = {193, -224, -131, 213};
//+
Surface(106) = {106};
//+
Line Loop(107) = {131, -162, -69, 151};
//+
Surface(107) = {107};
//+
Line Loop(108) = {69, -100, -31, 89};
//+
Surface(108) = {108};
//+
Line Loop(109) = {199, -213, -137, 212};
//+
Surface(109) = {109};
//+
Line Loop(110) = {137, -151, -75, 150};
//+
Surface(110) = {110};
//+
Line Loop(111) = {75, -89, -37, 88};
//+
Surface(111) = {111};
//+
Line Loop(112) = {173, -212, -111, 201};
//+
Surface(112) = {112};
//+
Line Loop(113) = {111, -150, -49, 139};
//+
Surface(113) = {113};
//+
Line Loop(114) = {49, -88, -11, 77};
//+
Surface(114) = {114};
//+
Line Loop(115) = {31, 32, -30, -21};
//+
Surface(115) = {115};
//+
Line Loop(116) = {30, 33, -29, -22};
//+
Surface(116) = {116};
//+
Line Loop(117) = {29, 34, -28, -23};
//+
Surface(117) = {117};
//+
Line Loop(118) = {28, 35, -27, -24};
//+
Surface(118) = {118};
//+
Line Loop(119) = {36, -26, -25, 27};
//+
Surface(119) = {119};
//+
Line Loop(120) = {37, 21, -20, 12};
//+
Surface(120) = {120};
//+
Line Loop(121) = {22, -19, 13, 20};
//+
Surface(121) = {121};
//+
Line Loop(122) = {19, 23, -18, 14};
//+
Surface(122) = {122};
//+
Line Loop(123) = {18, 24, -17, 15};
//+
Surface(123) = {123};
//+
Line Loop(124) = {25, -38, 16, 17};
//+
Surface(124) = {124};
//+
Line Loop(125) = {11, -12, -10, -1};
//+
Surface(125) = {125};
//+
Line Loop(126) = {9, 13, -10, 2};
//+
Surface(126) = {126};
//+
Line Loop(127) = {9, -14, -8, -3};
//+
Surface(127) = {127};
//+
Line Loop(128) = {7, 15, -8, 4};
//+
Surface(128) = {128};
//+
Line Loop(129) = {7, -16, -6, -5};
//+
Surface(129) = {129};
//+
Line Loop(130) = {69, 70, -68, -59};
//+
Surface(130) = {130};
//+
Line Loop(131) = {68, 71, -67, -60};
//+
Surface(131) = {131};
//+
Line Loop(132) = {67, 72, -66, -61};
//+
Surface(132) = {132};
//+
Line Loop(133) = {66, 73, -65, -62};
//+
Surface(133) = {133};
//+
Line Loop(134) = {65, 74, -64, -63};
//+
Surface(134) = {134};
//+
Line Loop(135) = {75, 59, -58, 50};
//+
Surface(135) = {135};
//+
Line Loop(136) = {58, 60, -57, 51};
//+
Surface(136) = {136};
//+
Line Loop(137) = {57, 61, -56, 52};
//+
Surface(137) = {137};
//+
Line Loop(138) = {56, 62, -55, 53};
//+
Surface(138) = {138};
//+
Line Loop(139) = {55, 63, -76, 54};
//+
Surface(139) = {139};
//+
Line Loop(140) = {49, -50, -48, -39};
//+
Surface(140) = {140};
//+
Line Loop(141) = {48, -51, -47, -40};
//+
Surface(141) = {141};
//+
Line Loop(142) = {47, -52, -46, -41};
//+
Surface(142) = {142};
//+
Line Loop(143) = {46, -53, -45, -42};
//+
Surface(143) = {143};
//+
Line Loop(144) = {45, -54, -44, -43};
//+
Surface(144) = {144};
//+
Line Loop(145) = {131, 132, -130, -121};
//+
Surface(145) = {145};
//+
Line Loop(146) = {130, 133, -129, -122};
//+
Surface(146) = {146};
//+
Line Loop(147) = {129, 134, -128, -123};
//+
Surface(147) = {147};
//+
Line Loop(148) = {128, 135, -127, -124};
//+
Surface(148) = {148};
//+
Line Loop(149) = {127, 136, -126, -125};
//+
Surface(149) = {149};
//+
Line Loop(150) = {137, 121, -120, 112};
//+
Surface(150) = {150};
//+
Line Loop(151) = {120, 122, -119, 113};
//+
Surface(151) = {151};
//+
Line Loop(152) = {119, 123, -118, 114};
//+
Surface(152) = {152};
//+
Line Loop(153) = {118, 124, -117, 115};
//+
Surface(153) = {153};
//+
Line Loop(154) = {117, 125, -138, 116};
//+
Surface(154) = {154};
//+
Line Loop(155) = {111, -112, -110, -101};
//+
Surface(155) = {155};
//+
Line Loop(156) = {110, -113, -109, -102};
//+
Surface(156) = {156};
//+
Line Loop(157) = {109, -114, -108, -103};
//+
Surface(157) = {157};
//+
Line Loop(158) = {108, -115, -107, -104};
//+
Surface(158) = {158};
//+
Line Loop(159) = {107, -116, -106, -105};
//+
Surface(159) = {159};
//+
Line Loop(160) = {188, -198, -189, 187};
//+
Surface(160) = {160};
//+
Line Loop(161) = {189, -197, -190, 186};
//+
Surface(161) = {161};
//+
Line Loop(162) = {190, -196, -191, 185};
//+
Surface(162) = {162};
//+
Line Loop(163) = {195, -191, -184, 192};
//+
Surface(163) = {163};
//+
Line Loop(164) = {194, -192, -183, 193};
//+
Surface(164) = {164};
//+
Line Loop(165) = {200, -187, -179, -178};
//+
Surface(165) = {165};
//+
Line Loop(166) = {179, -186, -180, -177};
//+
Surface(166) = {166};
//+
Line Loop(167) = {180, -185, -181, -176};
//+
Surface(167) = {167};
//+
Line Loop(168) = {181, -184, -182, -175};
//+
Surface(168) = {168};
//+
Line Loop(169) = {182, -183, -199, -174};
//+
Surface(169) = {169};
//+
Line Loop(170) = {168, 178, -169, 167};
//+
Surface(170) = {170};
//+
Line Loop(171) = {169, 177, -170, 166};
//+
Surface(171) = {171};
//+
Line Loop(172) = {176, -171, 165, 170};
//+
Surface(172) = {172};
//+
Line Loop(173) = {175, -172, 164, 171};
//+
Surface(173) = {173};
//+
Line Loop(174) = {174, -173, 163, 172};
//+
Surface(174) = {174};

//Surface Loop() = {,,,,,};

Surface Loop(1) = {1,16,114,105,125,140};
Surface Loop(2) = {6,25,113,104,140,155};
Surface Loop(3) = {7,26,112,103,155,174};
Surface Loop(4) = {17,105,96,2,126,141};
Surface Loop(5) = {24,104,95,8,141,156};
Surface Loop(6) = {27,103,94,15,156,173};
Surface Loop(7) = {18,3,96,87,127,142};
Surface Loop(8) = {23,9,95,86,142,157};
Surface Loop(9) = {28,14,94,85,157,172};
Surface Loop(10) = {19,4,87,78,128,143};
Surface Loop(11) = {21,10,86,77,143,158};
Surface Loop(12) = {29,13,85,76,158,171};
Surface Loop(13) = {20,5,78,69,129,144};
Surface Loop(14) = {22,11,77,68,144,159};
Surface Loop(15) = {30,12,76,67,159,170};
Surface Loop(16) = {35,16,111,102,120,135};
Surface Loop(17) = {36,25,110,101,135,150};
Surface Loop(18) = {45,26,109,100,150,169};
Surface Loop(19) = {34,17,102,93,121,136};
Surface Loop(20) = {37,24,101,92,136,151};
Surface Loop(21) = {44,27,100,91,151,168};
Surface Loop(22) = {33,18,93,84,122,137};
Surface Loop(23) = {38,23,92,83,137,152};
Surface Loop(24) = {43,28,91,82,152,167};
Surface Loop(25) = {32,19,84,75,123,138};
Surface Loop(26) = {39,21,83,74,138,153};
Surface Loop(27) = {42,29,82,73,153,166};
Surface Loop(28) = {31,20,75,64,124,139};
Surface Loop(29) = {40,22,74,65,139,154};
Surface Loop(30) = {41,30,73,66,154,165};

Surface Loop(31) = {50,35,108,99,115,130};
Surface Loop(32) = {51,36,107,98,130,145};
Surface Loop(33) = {60,45,106,97,145,164};

Surface Loop(34) = {49,34,99,90,116,131};
Surface Loop(35) = {52,37,98,89,131,146};
Surface Loop(36) = {59,44,97,88,146,163};

Surface Loop(37) = {48,33,90,81,117,132};
Surface Loop(38) = {53,38,89,80,132,147};
Surface Loop(39) = {58,43,88,79,147,162};

Surface Loop(40) = {47,32,81,72,118,133};
Surface Loop(41) = {54,39,80,71,133,148};
Surface Loop(42) = {57,42,79,70,148,161};

Surface Loop(43) = {46,31,72,63,119,134};
Surface Loop(44) = {55,40,71,62,134,149};
Surface Loop(45) = {56,41,70,61,149,160};


//Volume() = {};
Volume(1) = {1};
//Volume(2) = {2};
Volume(3) = {3};
Volume(4) = {4};
//Volume(5) = {5};
Volume(6) = {6};
Volume(7) = {7};
//Volume(8) = {8};
Volume(9) = {9};
Volume(10) = {10};
//Volume(11) = {11};
Volume(12) = {12};
Volume(13) = {13};
//Volume(14) = {14};
Volume(15) = {15};
Volume(16) = {16};
//Volume(17) = {17};
Volume(18) = {18};
Volume(19) = {19};
Volume(20) = {20};
Volume(21) = {21};
Volume(22) = {22};
//Volume(23) = {23};
Volume(24) = {24};
Volume(25) = {25};
Volume(26) = {26};
Volume(27) = {27};
Volume(28) = {28};
//Volume(29) = {29};
Volume(30) = {30};
Volume(31) = {31};
//Volume(32) = {32};
Volume(33) = {33};
Volume(34) = {34};
//Volume(35) = {35};
Volume(36) = {36};
Volume(37) = {37};
//Volume(38) = {38};
Volume(39) = {39};
Volume(40) = {40};
//Volume(41) = {41};
Volume(42) = {42};
Volume(43) = {43};
//Volume(44) = {44};
Volume(45) = {45};

//Physical Surface(1) = {8, 9, 10, 11, 12, 6, 7, 1, 2, 3, 4, 5};
//Physical Line(2) = {32, 33, 34, 35, 36, 26, 25, 17, 16, 6, 5, 4, 3, 2, 1, 11, 12, 20, 21, 31, 19, 23, 14, 18};
//+
Physical Volume(1) = {3, 6, 9, 12, 15, 30, 27, 24, 21, 18, 16, 19, 22, 25, 28, 26, 20, 1, 4, 7, 10, 13, 31, 34, 37, 40, 43, 45, 42, 39, 36, 33};
//+
Physical Surface(2) = {7,15,14,13,12,60,59,58,57,57,106,109,112,61,66,67,108,111,114,50:46,68,64,69,1:5,164:170,115:129,130,145,135,150,140,155,132,147,137,152,142,157,134,149,139,154,144,159,131,146,141,156,133,148,143,158};

Mesh.CharacteristicLengthMin = 5;
Mesh.CharacteristicLengthMax = 5;

//+
//Recombine Surface {1:12};

//+
Transfinite Line {1, 10, 11, 12, 88, 77, 49, 50, 48, 39, 87, 78} = 10 Using Progression 1;
Transfinite Surface {1,16,125,140,114,105};
Recombine Surface {1,16,125,140,114,105};
Transfinite Volume {1};

//+
Transfinite Line {101, 112, 111, 201, 212, 173, 174, 202, 172, 211, 110, 163} = 10 Using Progression 1;
Transfinite Surface {7,26,174,155,112,103};
Recombine Surface {7,26,174,155,112,103};
Transfinite Volume {3};
//+
Transfinite Line {178, 169, 167, 208, 205, 107, 116, 105, 206, 207, 168, 106} = 10 Using Progression 1;
Transfinite Surface {12,30,67,76,170,159};
Recombine Surface {12,30,67,76,170,159};
Transfinite Volume {15};
//+
Transfinite Line {54, 43, 45, 44, 83, 82, 6, 16, 5, 7, 84, 81} = 10 Using Progression 1;
Transfinite Surface {5,20,129,144,78,69};
Recombine Surface {5,20,129,144,78,69};
Transfinite Volume {13};
