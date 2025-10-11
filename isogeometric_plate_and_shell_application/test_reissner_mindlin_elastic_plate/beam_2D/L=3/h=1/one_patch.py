import sys
import os

import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *

E = 2.1e11
nu = 0.3
L = 3.0
h = 1.0
f = -1e6
plinear_solver = MKLPardisoSolver()

def test():
    order = 3
    nsampling = [40, 10]
    model1 = simulation_include.Model(E, nu, L, h, f, order, nsampling, plinear_solver)
    model = model1.Run(output=False)
    # for node in model.model_part.Nodes:
    #     print(node.Id)
    #     print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
    tol = 1e-11
    ref_ux = 0.000118934642103
    ref_uy = -0.000593273173087
    test_node = model.model_part.Nodes[559]
    ux = test_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
    uy = test_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
    assert(abs(ux - ref_ux) / abs(ux + ref_ux) < tol)
    assert(abs(uy - ref_uy) / abs(uy + ref_uy) < tol)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        order = 3
        nsampling = [40, 10]
        model1 = simulation_include.Model(E, nu, L, h, f, order, nsampling, plinear_solver)
        model1.Run(output=True)
