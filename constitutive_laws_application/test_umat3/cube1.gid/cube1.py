##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import cube1_include
from cube1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = cube1_include.Model('cube1',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    time = 2.0
    model.model_part.Properties[1].SetValue(DENSITY,         1.0e4 )
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(output=False, logging=False)

    mon_node = model.model_part.Nodes[4]
    uz = mon_node.GetSolutionStepValue(DISPLACEMENT_Z)
    ref_uz = -2.3357142857142856e-03
    print("uz: %.16e" % (uz))
    print("diff uz: %.16e" % (uz-ref_uz))
    assert(abs(uz - ref_uz) < 1e-10)

    print("Test passed")

def tag():
    return "umat"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
