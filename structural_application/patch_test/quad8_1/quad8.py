##################################################################
import sys
import os
import math
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import quad8_include
from quad8_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = quad8_include.Model('quad8',current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)

    model.model_part.Nodes[2].SetSolutionStepValue(FORCE_X, 1.0/3)
    model.model_part.Nodes[3].SetSolutionStepValue(FORCE_X, 1.0/3)
    model.model_part.Nodes[6].SetSolutionStepValue(FORCE_X, 4.0/3)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    ######### pytesting results #########
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            disp_x = node.GetSolutionStepValue(DISPLACEMENT_X)
            assert(abs(disp_x - 1.0) < 1.0e-12)
    #####################################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=False, output=False)
