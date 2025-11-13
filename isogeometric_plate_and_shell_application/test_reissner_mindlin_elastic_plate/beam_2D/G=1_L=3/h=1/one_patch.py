import sys
import os

import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

nu = 0.3
E = 2.0*(1+nu)
L = 3.0
h = 1.0
f = 1.0/h # fbar = 1
if KratosMKLSolversApplication.Has("MKLPardisoSolver"):
    plinear_solver = MKLPardisoSolver()
else:
    plinear_solver = SuperLUSolver()

def test():
    order = 2
    nsampling = [8, 2]
    model1 = simulation_include.Model(E, nu, L, h, f, order, nsampling, plinear_solver)
    model = model1.Run(output=False)
    # for node in model.model_part.Nodes:
    #     print(node.Id)
    #     print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
    tol = 1e-11
    ref_ux = -9.552463329985748
    ref_uy = 47.485131757202225
    test_node = model.model_part.Nodes[40]
    ux = test_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
    uy = test_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
    # print("ux: ", ux)
    # print("uy: ", uy)
    assert(abs(ux - ref_ux) / abs(ux + ref_ux) < tol)
    assert(abs(uy - ref_uy) / abs(uy + ref_uy) < tol)
    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        order = 3
        nsampling = [40, 10]
        model1 = simulation_include.Model(E, nu, L, h, f, order, nsampling, plinear_solver)
        model1.Run(output=True)
