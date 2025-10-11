##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
import math
##################################################################
##################################################################
import square10x10_include
from square10x10_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################
def main(output=True, logging=True):
    model = square10x10_include.Model('square10x10',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    temp = 10.0
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, temp)
    #    if abs(node.Y0 - 0.0) < tol:
    #        node.Fix(TEMPERATURE)
    #        node.SetSolutionStepValue(TEMPERATURE, 0)
    #    if abs(node.Y0 - 1.0) < tol:
    #        node.Fix(TEMPERATURE)
    #        node.SetSolutionStepValue(TEMPERATURE, temp)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if logging:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.5) < tol:
            temp = node.GetSolutionStepValue(TEMPERATURE)
            # print(temp)
            assert(abs(temp - 5.0) < 1e-6)

    print("Test passed")

def tag():
    return "thermal,amg"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
