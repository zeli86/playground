cl1 = 1;
a=10;
b=5;
Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};
Point(3) = {0, b,  0};
Point(4) = {-a, 0, 0};
Point(5) = {0, -b, 0};

// the first cutout:
Ellipse(1) = {2, 1, 5, 3};
Ellipse(2) = {3, 1, 5, 4};
Line(3) = {4,2};

Line Loop(1) = {1,2,3};

Plane Surface(1) = {1};

vol01[] = Extrude {{1, 0, 0}, {1, 0, 0}, Pi/2} {Surface{1};};
  
// these define the boundary indicators in deal.II:
//Physical Line(5) = {1, 2, 3, 4};
// you need the physical surface, because that is what deal.II reads in
//Plane Surface(1) = {1};
//Physical Surface(1) = {1};

//Volume(1) = {vol01};

/*
Transfinite Line {3} = 11 Using Progression 1;
Transfinite Line {9} = 11 Using Progression 1;
Transfinite Line {1,2,5,6} = 6 Using Progression 1;
Transfinite Volume{1} = {};
*/
/*
Transfinite Surface {1} = {};
Transfinite Surface {2} = {};
Transfinite Surface {17} = {};
Transfinite Surface {29} = {};
Transfinite Surface {12} = {};
Transfinite Surface {27} = {};
Transfinite Surface {15} = {};
Transfinite Volume{2} = {};
Recombine Surface {1,2,12,27,17,29,15};
*/

// some parameters for the meshing:
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.5;
Mesh.SubdivisionAlgorithm = 2;
Mesh.Smoothing = 10;
