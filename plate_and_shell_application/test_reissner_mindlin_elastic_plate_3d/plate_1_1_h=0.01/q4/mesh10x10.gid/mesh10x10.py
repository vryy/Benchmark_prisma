##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mo 28. Aug 14:29:34 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh10x10.gid')
import mesh10x10_include
from mesh10x10_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
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

    # debugging
    # for node in model.model_part.Nodes:
    #     node.Fix(DISPLACEMENT_X)
    #     node.Fix(DISPLACEMENT_Y)
        # node.Fix(ROTATION_Z)
    # end debugging

    ## DEBUGGING
    # for node in model.model_part.Nodes:
    #     node.Fix(ROTATION_X)
    #     node.Fix(ROTATION_Y)
    ## END OF DEBUGGING

    ## body force
    body_force = Vector(3)
    body_force[0] = 0
    body_force[1] = 0
    body_force[2] = -1.e2
    model.model_part.Properties[1].SetValue(BODY_FORCE, body_force)

    ## analysis
    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))
    #     print(node.GetSolutionStepValue(ROTATION))

    return model

def test():
    model = main(logging=False, output=False)

    ##### pytesting results

    psi_x = model.model_part.Nodes[4].GetSolutionStepValue(ROTATION_X)
    psi_y = model.model_part.Nodes[4].GetSolutionStepValue(ROTATION_Y)
    disp_z = model.model_part.Nodes[4].GetSolutionStepValue(DISPLACEMENT_Z)
    print("psi_x: %.10e, psi_y: %.10e, disp_z: %.10e" % (psi_x, psi_y, disp_z))
    ref_psi_x = 1.1496044229e-03
    ref_psi_y = -1.1496044229e-03
    ref_disp_z = -7.0049739520e-03
    assert(abs(disp_z - ref_disp_z) / abs(ref_disp_z) < 1e-10)
    assert(abs(psi_x - ref_psi_x) / abs(ref_psi_x) < 1e-10)
    assert(abs(psi_y - ref_psi_y) / abs(ref_psi_y) < 1e-10)

    ##### end pytesting results

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
