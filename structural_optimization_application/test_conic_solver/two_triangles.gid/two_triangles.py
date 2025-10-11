##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
#####     2023 by Hoang-Giang Bui                            #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Di 19. Sep 13:01:59 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import two_triangles_include
from two_triangles_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################
# Reference example:
#   /home/hbui/workspace/python/cvxpy/fela_socp.py
##################################################################

def main(output=True, logging=True):
    model = two_triangles_include.Model('two_triangles',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1e-06
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)

    time = 0.0
    model.SolveModel(time)
    # print(model.solver.solver.A)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(output=False, logging=False)

    lambda_ref = 2.6842866097108981e+02
    print("%.16e" % (model.solver.conic_solver.gamma))
    assert(abs(model.solver.conic_solver.gamma - lambda_ref) < 1e-10)

    print("Test passed")

def tag():
    return "FELA"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
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
