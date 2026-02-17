SetFactory("OpenCASCADE");

length = 40.0;
nelements = 41;

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {length, 0, 0, 1.0};

// Line
Line(1) = {1, 2};

// Optional: control number of elements
Transfinite Line {1} = nelements;

// Mesh only 1D
Mesh 1;
