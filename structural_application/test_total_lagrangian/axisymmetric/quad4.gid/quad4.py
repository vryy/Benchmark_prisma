##################################################################
import sys
import os
import math
##################################################################
import quad4_include
from quad4_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(logging=True, output=True):
    model = quad4_include.Model('quad4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)

    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            node.SetSolutionStepValue(FACE_LOAD_X, 1.0e3)
            node.SetSolutionStepValue(FACE_LOAD_Y, 0.0)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(output=False, logging=False)

    ######### pytesting results #########
    tol = 1e-6
    ref_disp_x = 0.000530094901587
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            disp_x = node.GetSolutionStepValue(DISPLACEMENT_X)
            assert(abs(disp_x - ref_disp_x) < 1e-12)
    #####################################
    print("Test passed")

def tag():
    tags = ""
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output = False)
