lh = 3;
n = 20;
progr = 1;

// exterior cube
Point(1) = {-lh,-lh,-lh};
Point(2) = {lh,-lh,-lh};
Point(3) = {lh,lh,-lh};
Point(4) = {-lh,lh,-lh};
Point(5) = {-lh,-lh,lh};
Point(6) = {lh,-lh,lh};
Point(7) = {lh,lh,lh};
Point(8) = {-lh,lh,lh};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {1,10,-5,-9};
Line Loop(4) = {2,11,-6,-10};
Line Loop(5) = {3,12,-7,-11};
Line Loop(6) = {4,9,-8,-12}; 

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};

// connection volumes
Surface Loop(1) = {1,2,3,4,5,6}; Volume(1) = {1};

// define transfinite mesh
//Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = n Using Progression progr; // polar richtung
//Transfinite Surface {1,2,3,4,5,6}={};
//Recombine Surface {1,2,3,4,5,6};
//Transfinite Volume {1}={};

Physical Surface(1) = {1,2,3,4,5,6}; 
Physical Volume(1) = {1};

