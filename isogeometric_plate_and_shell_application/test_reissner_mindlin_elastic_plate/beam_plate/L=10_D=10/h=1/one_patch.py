import sys
import os

import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

E = 2.1e11
nu = 0.3
L = 10.0
D = 10.0
h = 1.0
f = -1e6
if KratosMKLSolversApplication.Has("MKLPardisoSolver"):
    plinear_solver = MKLPardisoSolver()
else:
    plinear_solver = SuperLUSolver()

def test():
    order = 3
    nsampling = [40, 40]
    model1 = simulation_include.Model(E, nu, L, D, h, f, order, nsampling, plinear_solver)
    model = model1.Run(output=False)
    # for node in model.model_part.Nodes:
    #     print(node.Id)
    #     print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
    tol = 1e-8
    ref_w = -0.0674644630407
    test_node = model.model_part.Nodes[1849]
    w = test_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
    # print(w)
    test = abs(w - ref_w) / abs(w + ref_w)
    # print(test)
    assert(test < tol)
    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        order = 3
        nsampling = [40, 40]
        model1 = simulation_include.Model(E, nu, L, D, h, f, order, nsampling, plinear_solver)
        model1.Run()
