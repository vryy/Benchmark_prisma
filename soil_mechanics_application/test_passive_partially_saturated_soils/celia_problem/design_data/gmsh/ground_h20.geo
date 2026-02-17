SetFactory("OpenCASCADE");

width = 1.0;
height = 40.0;
nelements = 40;

// Base surface
Rectangle(1) = {0, 0, 0, width, width};

// 1 element in x and y
Transfinite Curve{1,2,3,4} = 2;
Transfinite Surface{1};
Recombine Surface{1};

// Extrude → creates a VOLUME with hexes
vol[] = Extrude {0, 0, height} {
  Surface{1};
  Layers{nelements};
  Recombine;
};

// Make volume transfinite + recombined
Transfinite Volume{vol[1]};
Recombine Volume{vol[1]};

// HEX20
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh.HighOrderOptimize = 1;

// Force 3D meshing
Mesh 3;
