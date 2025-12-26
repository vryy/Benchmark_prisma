##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import square1_include as simulation_include
from square1_include import *
model_name_ = 'square1'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    # ============================================ #
    # |       USER CALCULATION SCRIPT            | #
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv #

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol):
            node.Fix(PRESSURE)
            node.SetSolutionStepValue(PRESSURE, 1.0)
        if (abs(node.X0 - 1.0) < tol) or (abs(node.Y0) < tol) or (abs(node.Y0 - 1.0) < tol):
            node.Fix(PRESSURE)
            node.SetSolutionStepValue(PRESSURE, 0.0)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 0.5) < tol) and (abs(node.Y0 - 0.5) < tol):
            mon_node = node

    pres = mon_node.GetSolutionStepValue(PRESSURE)
    # print("%.16e" % pres)
    ref_pres = 2.3254786506765721e-01
    assert(abs(pres - ref_pres) < 1e-10)
    print("Test passed")

def tag():
    return "helmholtz"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
