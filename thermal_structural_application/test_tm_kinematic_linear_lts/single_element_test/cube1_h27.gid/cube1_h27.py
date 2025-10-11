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
## This file is generated on Sa 21. Mar 11:46:17 CET 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./cube1_h27.gid')
import cube1_h27_include
from cube1_h27_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=False, logging=False):

    model = cube1_h27_include.Model('cube1_h27',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## set reference temperature
    ref_temp = 0.0 # 300.0 # 22 oC
    values = [ref_temp]*27

    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(REFERENCE_TEMPERATURE, values, model.model_part.ProcessInfo)

    ## boundary condition
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(TEMPERATURE)
            prescribed_nodes.append(node)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
        node.SetSolutionStepValue(TEMPERATURE, ref_temp)

    ## analysis step 0

    time = 0.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    ## analysis step 1

    time = 1.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(TEMPERATURE, ref_temp + 50.0)

    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    ## reporting
    if logging:
        react_p = 0.0
        for node in prescribed_nodes:
            print("Node %d: react_x = %f, disp_x = %f" % (node.Id, node.GetSolutionStepValue(REACTION_X), node.GetSolutionStepValue(DISPLACEMENT_X)))
            react_p += node.GetSolutionStepValue(REACTION_X)
        work_done = 0.5 * react_p * node.GetSolutionStepValue(DISPLACEMENT_X)
        print("work done = %f" % (work_done))
        strain_energy = model.solver.solver.GetStrainEnergy()
        print("strain energy = %f" % (strain_energy))

    return model

def test():
    model = main(logging=False, output=False)

    ### pytesting results
    tol = 1e-06
    assert(abs(model.model_part.Nodes[11].GetSolutionStepValue(DISPLACEMENT_X) - 0.000650) < 1e-10)
    assert(abs(model.model_part.Nodes[11].GetSolutionStepValue(DISPLACEMENT_Y) - 0.000650) < 1e-10)
    assert(abs(model.model_part.Nodes[11].GetSolutionStepValue(DISPLACEMENT_Z) - 0.000650) < 1e-10)
    print("Test passed")

def tag():
    return "thermal-mechanical"

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
