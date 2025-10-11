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
## This file is generated on Do 30. Jun 20:03:30 CEST 2022 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./quad9.gid')
import quad9_include
from quad9_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################
def main(output=True, logging=True):
    model = quad9_include.Model('quad9',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    if logging:
        ifile = open("Disp_Reaction.txt", 'w')

    ## boundary condition
    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol):
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    if logging:
        ifile.write('node\tdx\tdy\tdz\trx\try\trz\n')
        for node in prescribed_nodes:
            dx = node.GetSolutionStepValue(DISPLACEMENT_X)
            dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
            dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
            rx = node.GetSolutionStepValue(REACTION_X)
            ry = node.GetSolutionStepValue(REACTION_Y)
            rz = node.GetSolutionStepValue(REACTION_Z)
            ifile.write(str(node.Id) + '\t' + str(dx) + '\t' + str(dy) + '\t' + str(dz) + '\t' + str(rx) + '\t' + str(ry) + '\t' + str(rz) + '\n')
        ifile.write('\n')
        ifile.flush()

    step_list = []
    step_list.append([8, 1.0e-4])
    step_list.append([100, 1.0e-5])
    step_list.append([100, 1.0e-5])
    step_list.append([100, 1.0e-5])

    step = 0
    for tmp in step_list:
        nsteps = tmp[0]
        delta_disp = tmp[1]
        delta_time = delta_disp

        for i in range(0, nsteps):

            for node in prescribed_nodes:
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)

            time += delta_time

            print("Load step " + str(step+1) + ", time " + str(time) + " starts")

            model.SolveModel(time)
            if output:
                model.WriteOutput(time)

            print("Load step " + str(step+1) + ", time " + str(time) + " completed")

            if logging:
                ifile.write('node\tdx\tdy\tdz\trx\try\trz\n')
                for node in prescribed_nodes:
                    dx = node.GetSolutionStepValue(DISPLACEMENT_X)
                    dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
                    dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
                    rx = node.GetSolutionStepValue(REACTION_X)
                    ry = node.GetSolutionStepValue(REACTION_Y)
                    rz = node.GetSolutionStepValue(REACTION_Z)
                    ifile.write(str(node.Id) + '\t' + str(dx) + '\t' + str(dy) + '\t' + str(dz) + '\t' + str(rx) + '\t' + str(ry) + '\t' + str(rz) + '\n')
                ifile.write('\n')
                ifile.flush()

            step += 1

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    ### pytesting results
    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            prescribed_nodes.append(node)

    reac_force_x = 0.0
    for node in prescribed_nodes:
        reac_force_x += node.GetSolutionStepValue(REACTION_X)
    ref_reac = 2.137171991876853e+01
    print("%.15e" % reac_force_x)
    assert(abs(reac_force_x - ref_reac) / ref_reac < 1e-10)

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
