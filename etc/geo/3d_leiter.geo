l = 1;
r1 = 3;
r2 = 1;
r3 = 0.5;
n = 20;
n2 = 20;
progr = 1;

// exterior cube
Point(1) = {-r1,0,-r1};
Point(2) = {0,r1,-r1};
Point(3) = {r1,0,-r1};
Point(4) = {0,-r1,-r1};
Point(5) = {0,r2,-r1};
Point(6) = {0,-r2,-r1};
Point(7) = {r2,0,-r1};
Point(8) = {-r2,0,-r1};
Point(9) = {0,0,-r1};
Point(10) = {0,r3,-r1};
Point(11) = {0,-r3,-r1};
Point(12) = {r3,0,-r1};
Point(13) = {-r3,0,-r1};

Line(1) = {1,8};
Line(2) = {2,5};
Line(3) = {3,7};
Line(4) = {4,6};
 
Circle(5) = {5,9,7};
Circle(6) = {7,9,6};
Circle(7) = {6,9,8};
Circle(8) = {8,9,5};
Circle(9) = {1,9,2};
Circle(10) = {2,9,3};
Circle(11) = {3,9,4};
Circle(12) = {4,9,1};

Line(13) = {11,12};
Line(14) = {12,10};
Line(15) = {10,13};
Line(16) = {13,11};
Line(17) = {7,12};
Line(18) = {5,10};
Line(19) = {8,13};
Line(20) = {6,11};

Line Loop(2) = {9,2,-8,-1};
Line Loop(3) = {10,3,-5,-2};
Line Loop(4) = {11,4,-6,-3};
Line Loop(5) = {12,1,-7,-4};
Line Loop(6) = {13,14,15,16}; // inneres Quadrat
Line Loop(7) = {13,-17,6,20};
Line Loop(8) = {14,-18,5,17};
Line Loop(9) = {15,-19,8,18};
Line Loop(10) = {16,-20,7,19};

Plane Surface(2) = {2};
Extrude Surface {2, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(3) = {3};
Extrude Surface {3, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(4) = {4};
Extrude Surface {4, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(5) = {5};
Extrude Surface {5, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(6) = {6};
Extrude Surface {6, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(7) = {7};
Extrude Surface {7, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(8) = {8};
Extrude Surface {8, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(9) = {9};
Extrude Surface {9, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

Plane Surface(10) = {10};
Extrude Surface {10, {0.0,0.0,2*r1}};
Delete { Volume{1}; }

 
// // connection volumes
Surface Loop(1) = {6,130,117,121,129,125}; Volume(1) = {1};
Surface Loop(2) = {10,218,129,187,103,151}; Volume(2) = {2};
Surface Loop(3) = {7,152,117,151,81,143}; Volume(3) = {3};
Surface Loop(4) = {8,174,121,165,59,143}; Volume(4) = {4};
Surface Loop(5) = {9,196,187,125,165,37}; Volume(5) = {5};
Surface Loop(6) = {5,108,95,103,77,41}; Volume(6) = {6};
Surface Loop(7) = {4,86,81,77,73,55}; Volume(7) = {7};
Surface Loop(8) = {3,64,59,51,55,33}; Volume(8) = {8};
Surface Loop(9) = {33,41,2,42,37,29}; Volume(9) = {9};

// define transfinite mesh
Transfinite Line {5,6,7,8,13,14,15,16,9,10,11,12,88,90,113,110,68,66,111,46,44,112,24,22} = 15 Using Progression progr; // polar richtung
Transfinite Line {1,2,3,4,17,18,19,20,155,23,177,25,135,67,133,45} = 15 Using Progression progr; // radial richtung
Transfinite Line {72,76,27,36,115,124,116,120,54,50,32,28} = 15 Using Progression progr; // z richtung

//Transfinite Surface {2,3,4,5,6,7,8,9,10,29,41,29,41,37,187,95,103,129,125,165,33,121,59,51,117,151,77,143,81,55,73,42,108,218,152,86,130,196,174,64}={};
//Recombine Surface {2,3,4,5,6,7,8,9,10,29,41,29,41,37,187,95,103,129,125,165,33,121,59,51,117,151,77,143,81,55,73,42,108,218,152,86,130,196,174,64};
//Transfinite Volume {1,2,3,4,5,6,7,8,9}={};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Surface(1) = {29,51,73,95}; // Mantelfläche
Physical Surface(2) = {6,7,8,9,10}; // kleiner Kreis 
Physical Surface(3) = {130,174,196,152,218}; //kleiner Kreis
Physical Surface(4) = {2,3,4,5,42,108,86,64}; 
//Physical Volume(1) = {1,2,3,4,5,6,7,8,9};
Physical Volume(1) = {1,2,3,4,5};
Physical Volume(2) = {6,7,8,9};
