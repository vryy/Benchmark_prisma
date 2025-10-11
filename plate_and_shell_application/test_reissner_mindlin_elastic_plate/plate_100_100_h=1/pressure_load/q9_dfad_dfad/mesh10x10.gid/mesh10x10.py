##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Fr 21. Jun 11:32:20 CEST 2024
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import mesh10x10_include
from mesh10x10_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = mesh10x10_include.Model('mesh10x10',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1e-6
    w = 100.0
    h = 100.0
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) or (abs(node.X0 - w) < tol) or (abs(node.Y0) < tol) or (abs(node.Y0 - h) < tol):
            node.Fix(DISPLACEMENT_Z)
            node.Fix(ROTATION_X)
            node.Fix(ROTATION_Y)

    ## DEBUGGING
    # for node in model.model_part.Nodes:
    #     node.Fix(ROTATION_X)
    #     node.Fix(ROTATION_Y)
    ## END OF DEBUGGING

    ## pressure load
    p = 1e2
    for element in model.model_part.Elements:
        element.SetValue(POSITIVE_FACE_PRESSURE, 0.5*p)
        element.SetValue(NEGATIVE_FACE_PRESSURE, -0.5*p)

    ## analysis
    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    return model

def test():
    model = main(logging=False, output=False)

    ##########pytesting results
    tol = 1e-6
    w = 100.0
    h = 100.0
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 0.5*w) < tol) and (abs(node.Y0 - 0.5*h) < tol):
            test_node = node
    disp_z = test_node.GetSolutionStepValue(DISPLACEMENT_Z)
    print("%.15e" % disp_z)
    ref_disp_z = -6.713149224870795e+02
    error = abs(disp_z - ref_disp_z) / abs(ref_disp_z)
    assert(error < 1e-10)
    print("Test passed")
    ###########################

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
