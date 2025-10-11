##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Sa 14. Mar 00:15:32 CET 2020
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import cube1_include
from cube1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = cube1_include.Model('cube1',current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(COMPLEX_DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(COMPLEX_DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(COMPLEX_DISPLACEMENT_Z)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(COMPLEX_DISPLACEMENT_X)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    time = 1.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(COMPLEX_DISPLACEMENT_X, 0.1)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    react_p = 0.0
    for node in prescribed_nodes:
        # print("Node %d: react_x = %f, disp_x = %f" % (node.Id, node.GetSolutionStepValue(COMPLEX_REACTION_X), node.GetSolutionStepValue(COMPLEX_DISPLACEMENT_X)))
        # print(node.Id)
        # print(node.GetSolutionStepValue(COMPLEX_DISPLACEMENT))
        # print(node.GetSolutionStepValue(COMPLEX_REACTION))
        react_p += node.GetSolutionStepValue(COMPLEX_REACTION_X)
    work_done = 0.5 * react_p * node.GetSolutionStepValue(COMPLEX_DISPLACEMENT_X)
    # print("work done = %.10e" % (work_done))
    print("work_done: ", work_done)
    # strain_energy = model.solver.solver.GetStrainEnergy()
    # # print("strain energy = %.10e" % (strain_energy))
    # print("strain energy: ", strain_energy)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    ######### pytesting results #########
    # assert(abs(strain_energy - 10500.0) / strain_energy < 1e-12)
    assert(abs((work_done - 10500.0) / work_done) < 1e-12)
    #####################################

def test():
    main(logging=False, output=False)
    print("Test passed")

def tag():
    return "kinematic_linear,linear_elastic,gcomplex"

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
