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
## This file is generated on Mi 12. Jan 11:18:00 CET 2022 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh1.gid')
import mesh1_include
from mesh1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = mesh1_include.Model('mesh1',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    # boundary condition
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Y0-250.0) < tol:
            node1 = node
        if abs(node.X0-250.0) < tol and abs(node.Y0) < tol:
            node2 = node
        if abs(node.X0) < tol and abs(node.Y0-500.0) < tol:
            node3 = node
        if abs(node.X0-500.0) < tol and abs(node.Y0) < tol:
            node4 = node

    node1.Fix(DISPLACEMENT_X)
    node1.Fix(DISPLACEMENT_Y)
    node2.Fix(DISPLACEMENT_X)
    node2.Fix(DISPLACEMENT_Y)

    # colinear constraint on edges
    prop = model.model_part.Properties[2]
    last_cond_id = len(model.model_part.Conditions)
    for node_id in model.layer_nodes_sets['left']:
        if not (node_id == node1.Id) and not (node_id == node3.Id):
            node = model.model_part.Nodes[node_id]
            cond = CollinearConstraint2D(last_cond_id + 1, node1, node, node3, prop)
            model.model_part.AddCondition(cond)
            last_cond_id = last_cond_id + 1

    for node_id in model.layer_nodes_sets['bottom']:
        if not (node_id == node2.Id) and not (node_id == node4.Id):
            node = model.model_part.Nodes[node_id]
            cond = CollinearConstraint2D(last_cond_id + 1, node2, node, node4, prop)
            model.model_part.AddCondition(cond)
            last_cond_id = last_cond_id + 1

    # sys.exit(0)
    # time stepping

    time = 0.0
    model.SolveModel(time)

    disp = 0.0
    delta_disp = 1.0e-2
    delta_time = delta_disp

    if logging:
        ifile = open('reaction_result.txt', "w")
        ifile.write('u\tf\n')
        ifile.write('0\t0\n')

    for n in range(0, 250):
        disp = disp + delta_disp
        node1.SetSolutionStepValue(DISPLACEMENT_Y, disp)
        node2.SetSolutionStepValue(DISPLACEMENT_X, disp)

        time = time + delta_time
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            ifile.write(str(disp) + "\t" + str(node1.GetSolutionStepValue(REACTION_Y)) + "\n")
            ifile.flush()

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Y0-250.0) < tol:
            node1 = node

    ###### pytesting results
    reac = node1.GetSolutionStepValue(REACTION_Y)
    ref_reac = 13.43852898311006
    print("%.15e" % reac)
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    print("Test passed")
    ########################

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
