cl1 = 1;
a=10;
b=5;
Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};
Point(3) = {0, b,  0};
Point(4) = {-a, 0, 0};
Point(5) = {0, -b, 0};

Ellipse(1) = {2, 1, 5, 3};
Ellipse(2) = {3, 1, 5, 4};
Ellipse(3) = {5, 1, 2, 4};
Ellipse(4) = {2, 1, 4, 5};
Line(5) = {1,2};
Line(6) = {1,3};
Line(7) = {1,4};
Line(8) = {1,5};

Line Loop(1) = {1, 5, -6};
Line Loop(2) = {2, -7, 6};
Line Loop(3) = {3, -7, 8};
Line Loop(4) = {4, -8, 5};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

//Transfinite Line {1,2,3,4} = 11 Using Progression 1;
//Transfinite Line {5,6,7,8} = 11 Using Progression 1;
//Transfinite Surface {1} = {1,2,3};
//Transfinite Surface {2} = {1,3,4};
//Transfinite Surface {3} = {1,4,5};
//Transfinite Surface {4} = {1,2,5};

Recombine Surface{1};

// these define the boundary indicators in deal.II:
Physical Line(1) = {1, 2, 3, 4};
// you need the physical surface, because that is what deal.II reads in
Physical Surface(1) = {1,2,3,4};

// some parameters for the meshing:
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.5;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;


Coherence;
