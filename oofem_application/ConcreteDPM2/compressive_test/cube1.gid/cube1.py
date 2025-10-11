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
# This example tests the tensile part of C-DPM2. The stress will
# reach and bounded by ft
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
    model = cube1_include.Model('cube1',current_dir_,current_dir_,logging)
    model.InitializeModel()

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    if logging:
        ifile = open("reaction.log", "w")
        ifile.write("%-*s%s\n" % (20, "disp", "x-reaction"))

    disp = 0.0
    delta_disp = -3e-6
    delta_time = 3e6*delta_disp

    for i in range(0, 300): #1):
        print("#######################################")
        print(f"##########LOADING STEP {i+1} STARTED##############")
        print("#######################################")

        time += delta_time
        disp += delta_disp
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            react_p = 0.0
            for node in prescribed_nodes:
                react_p += node.GetSolutionStepValue(REACTION_X)
            ifile.write("%-*.10e%.10e\n" % (20, disp, react_p))
            ifile.flush()

    if logging:
        ifile.close()

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            prescribed_nodes.append(node)

    react_p = 0.0
    for node in prescribed_nodes:
        react_p += node.GetSolutionStepValue(REACTION_X)
    print("%.16e, %.3e %%" % (react_p, abs(react_p + 3e6) / 3e6 * 100))

    ref_reac = -2.9770397441691551e+06
    assert(abs(react_p - ref_reac) / abs(ref_reac) < 1e-8)

    print("Test passed")

def tag():
    return "oofem,concrete,cdpm2"

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
