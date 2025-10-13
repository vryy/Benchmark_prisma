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
## This file is generated on Sa 26. Mar 13:54:04 CET 2022 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh_1.gid')
import mesh_1_include
from mesh_1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, output_more=False, logging=True, to_continue=True):
    model = mesh_1_include.Model('mesh_1',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 200.0) < tol and abs(node.Y0 + 10.0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 420.0) < tol and abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)

    load_node_1 = model.model_part.Nodes[model.layer_nodes_sets['load_1'][0]]
    load_node_2 = model.model_part.Nodes[model.layer_nodes_sets['load_2'][0]]
    measured_nodes = [load_node_1, load_node_2]

    ## create point load
    model.model_part.CreateNewCondition("PointForce2D", 1, [load_node_1.Id], model.model_part.Properties[1])
    model.model_part.CreateNewCondition("PointForce2D", 2, [load_node_2.Id], model.model_part.Properties[1])
    # print(model.model_part)

    ## callback function to set external load
    P=-1e3
    def SetLoad(model_part, lmbda):
        load_node_1.SetSolutionStepValue(FORCE_X, 0.0)
        load_node_1.SetSolutionStepValue(FORCE_Y, lmbda*P*10.0/11)
        load_node_2.SetSolutionStepValue(FORCE_X, 0.0)
        load_node_2.SetSolutionStepValue(FORCE_Y, lmbda*P/11)
    model.solver.solver.set_load_callback = SetLoad

    ## initialization
    time = 0.0
    model.SolveModel(time)

    if logging:
        ifile = open("Disp_Reaction.txt", 'w')

        ifile.write('node\tdx\tdy\tdz\trx\try\trz\n')
        for node in measured_nodes:
            dx = node.GetSolutionStepValue(DISPLACEMENT_X)
            dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
            dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
            rx = node.GetSolutionStepValue(REACTION_X)
            ry = node.GetSolutionStepValue(REACTION_Y)
            rz = node.GetSolutionStepValue(REACTION_Z)
            ifile.write(str(node.Id) + '\t' + str(dx) + '\t' + str(dy) + '\t' + str(dz) + '\t' + str(rx) + '\t' + str(ry) + '\t' + str(rz) + '\n')
        ifile.write('\n')
        ifile.flush()

    ## loading
    step_list = []
    step_list.append([3, 5.0e-1])   # 1.5
    step_list.append([10, 1.0e-1])  # 2.5
    step_list.append([10, 5.0e-2])  # 3.0
    step_list.append([10, 5.0e-2])  # 3.5
    step_list.append([10, 5.0e-2])  # 4.0
    step_list.append([20, 1.0e-1])  # 6.0
    step_list.append([10, 2.0e-1])  # 8.0
    step_list.append([10, 2.0e-1])  # 10.0
    step_list.append([10, 5.0e-1])  # 15.0
    step_list.append([10, 5.0e-1])  # 20.0
    step_list.append([10, 5.0e-1])  # 25.0
    step_list.append([10, 1.0])  # 35.0
    step_list.append([10, 1.0])  # 45.0
    step_list.append([3, 1.0])   # 48.0
    step = 0
    cnt_list = 0
    for tmp in step_list:
        nsteps = tmp[0]
        model.solver.solver.arc_length_control_process.GetConstraint().SetRadius(tmp[1])
        delta_time = tmp[1]
        for i in range(0, nsteps):
            time += delta_time
            step += 1

            print("Loading step " + str(step) + ", time " + str(time) + " started")
            model.SolveModel(time)
            if output_more:
                model.WriteOutput(time)
            print("Loading step " + str(step) + ", time " + str(time) + " completed")

            if logging:
                ifile.write('node\tdx\tdy\tdz\trx\try\trz\n')
                for node in measured_nodes:
                    dx = node.GetSolutionStepValue(DISPLACEMENT_X)
                    dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
                    dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
                    rx = node.GetSolutionStepValue(REACTION_X)
                    ry = node.GetSolutionStepValue(REACTION_Y)
                    rz = node.GetSolutionStepValue(REACTION_Z)
                    ifile.write(str(node.Id) + '\t' + str(dx) + '\t' + str(dy) + '\t' + str(dz) + '\t' + str(rx) + '\t' + str(ry) + '\t' + str(rz) + '\n')
                ifile.write('\n')
                ifile.flush()

        if output:
            model.WriteOutput(time)

        cnt_list += 1

        if to_continue:
            if cnt_list == len(step_list):
                cont = input("Do you want to continue (y/n)?")
                if cont == "y":
                    tmp = input("Please give the values to continue [nsteps, radius]:")
                    step_list.append(tmp)
                else:
                    print("Analysis is completed")

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False, output_more=False, to_continue=False)

    ###### pytesting results
    load_node_1 = model.model_part.Nodes[model.layer_nodes_sets['load_1'][0]]
    load_node_2 = model.model_part.Nodes[model.layer_nodes_sets['load_2'][0]]
    reac_1 = load_node_1.GetSolutionStepValue(REACTION_Y)
    reac_2 = load_node_2.GetSolutionStepValue(REACTION_Y)
    ref_reac_1 = -2375.5126235417656
    ref_reac_2 = -237.5512623559085
    print("reac_1: %.15e, reac_2: %.15e" % (reac_1, reac_2))
    assert(abs(reac_1 - ref_reac_1) / abs(ref_reac_1) < 1e-10)
    assert(abs(reac_2 - ref_reac_2) / abs(ref_reac_2) < 1e-10)
    print("Test passed")
    ########################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True, output_more=False, to_continue=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
