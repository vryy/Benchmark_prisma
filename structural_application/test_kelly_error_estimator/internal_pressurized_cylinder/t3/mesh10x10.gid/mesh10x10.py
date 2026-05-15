##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import mesh10x10_include
from mesh10x10_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = mesh10x10_include.Model('mesh10x10',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    P = 1.0e2
    load = Vector(3)
    for node in model.model_part.Nodes:
        x = node.X0
        y = node.Y0
        r = math.sqrt(x*x + y*y)
        t = math.acos(x/r)
        load[0] = P * math.cos(t)
        load[1] = P * math.sin(t)
        load[2] = 0.0
        node.SetSolutionStepValue(FACE_LOAD, load)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(output=False, logging=False)

    stress_util = RecoverStressUtility()
    ee = stress_util.ComputeKellyErrorEstimation(model.model_part, DISPLACEMENT)

    print("ee: %.16e" % (ee))
    ee_ref = 1.1756147477114069e-03
    assert(abs(ee - ee_ref) < 1e-10)

    print("Test passed")

def tag():
    return "kelly"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
