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
## This file is generated on Di 29. Aug 11:14:37 CEST 2023
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
model = quad4_2_include.Model('quad4_2',current_dir_,current_dir_)
model.InitializeModel()
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main():
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

def test():
    main()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
