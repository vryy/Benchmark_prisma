##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import quad4_2_include
from quad4_2_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = quad4_2_include.Model('quad4_2',current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Y0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    time = 1.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Y, -0.1)
    model.SolveModel(time)

    react_y = 0.0
    for node in prescribed_nodes:
        ry = node.GetSolutionStepValue(REACTION_Y)
        print("Reaction at node %d: %.16e" % (node.Id, ry))
        react_y += ry
    print("Sum of reaction forces: %.16e" % react_y)

    ######### pytesting results #########
    ref_reac_y = -2.3076923076923119e+10
    test = abs(react_y - ref_reac_y) / abs(ref_reac_y)
    assert(test < 1e-12)
    #####################################
    print("Test passed")

def test():
    main(output=False, logging=False)

def tag():
    tags = ""
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=True, logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
