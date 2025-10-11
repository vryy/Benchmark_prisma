import sys
import os

import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *

E = 2.1e11
nu = 0.3
L = 10.0
D = 10.0
h = 1.0e-2
f = -1e1
plinear_solver = MKLPardisoSolver()

def test():
    order = 3
    nsampling = [40, 40]
    model1 = simulation_include.Model(E, nu, L, D, h, f, order, nsampling, plinear_solver)
    model = model1.Run(output=False)
    # for node in model.model_part.Nodes:
    #     print(node.Id)
    #     print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
    tol = 1e-8
    ref_w = -0.66159505807
    test_node = model.model_part.Nodes[1849]
    w = test_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
    # print(w)
    test = abs(w - ref_w) / abs(w + ref_w)
    # print(test)
    assert(test < tol)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        order = 3
        nsampling = [40, 40]
        model1 = simulation_include.Model(E, nu, L, D, h, f, order, nsampling, plinear_solver)
        model1.Run()
