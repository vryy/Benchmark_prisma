##################################################################
import sys
import os
import math
##################################################################
import two_triangles_include
from two_triangles_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = two_triangles_include.Model('two_triangles',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    #model.WriteOutput(time)

    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            node.SetSolutionStepValue(FORCE_X, 1.0e3)
            node.SetSolutionStepValue(FORCE_Y, 0.0)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT_X))
    #     print(node.GetSolutionStepValue(DISPLACEMENT_Y))

    return model

def test():
    model = main(logging=False, output=False)

    work_done = 0.0
    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            work_done += 0.5 * 1.0e3 * node.GetSolutionStepValue(DISPLACEMENT_X)
    print("work_done: %.10e" % (work_done))
    strain_energy = model.solver.solver.GetStrainEnergy()
    print("strain energy = %.10e" % (strain_energy))

    ######### pytesting results #########
    ref_energy = 6.6666666667e-02
    assert(abs(strain_energy - ref_energy) < 1e-12)
    assert(abs(work_done - ref_energy) < 1e-12)
    #####################################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
