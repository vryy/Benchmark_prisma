##################################################################
import sys
import os
import time
##################################################################
##################################################################
model_path = os.getcwd() + "/../design_data/slope_stability_hyplas.gid"
sys.path.append(model_path)
import slope_stability_hyplas_include
from slope_stability_hyplas_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

sys.path.append(os.getcwd() + "/../")
import simulator

def main(output=True, logging=True, case_number=7):
    model = slope_stability_hyplas_include.Model('slope_stability_hyplas',model_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    simulator.Run(model, output=output, case_number=case_number, coarsening=True)

    return model

def test():
    model = main(output=False, logging=False, case_number=7)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 35.0) < tol and abs(node.Y0 - 40.) < tol:
            pointA = node

    ux = pointA.GetSolutionStepValue(DISPLACEMENT_X)
    uy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    ux_ref = 4.4788988252431611e-01
    uy_ref = -5.6328076138229455e-01
    print("ux: %.16e, uy: %.16e" % (ux, uy))
    assert(abs(ux - ux_ref) < 1e-10)
    assert(abs(uy - uy_ref) < 1e-10)
    print("Test passed")

def tag():
    return "p4est"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True)
