##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import quad9_include
from quad9_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = quad9_include.Model('quad9',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)

    model.model_part.Nodes[9].SetSolutionStepValue(FACE_LOAD_X, 2.0)
    model.model_part.Nodes[6].SetSolutionStepValue(FACE_LOAD_X, 2.0)
    model.model_part.Nodes[8].SetSolutionStepValue(FACE_LOAD_X, 2.0)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False,output=False)

    ######### pytesting results #########
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            disp_x = node.GetSolutionStepValue(DISPLACEMENT_X)
            # print(disp_x)
            assert(abs(disp_x - 1.0) < 1.0e-12)
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
        main(logging=True,output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
