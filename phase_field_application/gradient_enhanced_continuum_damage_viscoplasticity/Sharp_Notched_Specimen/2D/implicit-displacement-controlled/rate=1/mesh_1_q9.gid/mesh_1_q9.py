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
## This file is generated on Do 29. Sep 19:04:00 CEST 2022 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh_1_q9.gid')
import mesh_1_q9_include
from mesh_1_q9_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True, nsteps=1000, dstep=4e-7):

    model = mesh_1_q9_include.Model('mesh_1_q9',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    tol = 1.0e-6
    prescribed_nodes = []

    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) and (abs(node.Y0) < tol):
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Y0 - 2.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    if logging:
        ifile = open("Disp_Reaction.txt", 'w')
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

    # rate = 1.0 # [mm/s]
    # rate = 1.0e3 # [mm/s]
    # rate = 1.0e6 # [mm/s]
    #rate = 1.0e-3 # [mm/s]
    rate = 1.0 # [mm/s]

    step_list = []
    step_list.append([nsteps, dstep])

    step = 0
    for tmp in step_list:
        nsteps = tmp[0]
        delta_disp = tmp[1]
        delta_time = tmp[1] / rate

        for i in range(0, nsteps):

            for node in prescribed_nodes:
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, delta_disp)

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

    if logging and (not output):
        model.WriteOutput(time) # just write the last result

    if logging:
        ifile.close()

    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False, nsteps=10, dstep=4e-5)

    ### pytesting results
    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 2.0) < tol:
            prescribed_nodes.append(node)

    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ref_reac = 1.158853480303876e-01
    print("%.15e" % reac_force_y)
    assert(abs(reac_force_y - ref_reac) / ref_reac < 1e-10)

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
