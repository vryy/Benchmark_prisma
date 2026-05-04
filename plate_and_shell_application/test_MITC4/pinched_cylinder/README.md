The geometry and simulation data is taken from the paper:

[1]   Kiendl et al, Isogeometric shell analysis with Kirchhoff–Love elements, 2009

The analytical solution provided by Kiendl for nu=0.3 is 1.8248e-5. He obtained 1.8264e-5 in the analysis.

Another more accurate analytical solution is from:

[2]   Yamamoto et al, A quadrilateral shell element incorporating thickness-stretch for nearly incompressible hyperelasticity, 2019

where it is 1.826797e-5 for nu=0.3 (so Kiendl is correct than he's used to imagine:))

As of 2022, the solution obtained by MITC4 is 1.87181654621e-05 using the 160x160 mesh. From [2], the MITC4 results converge quite closed to the reference solution. Therefore, it is suggested that problem still existed in the implementation.
