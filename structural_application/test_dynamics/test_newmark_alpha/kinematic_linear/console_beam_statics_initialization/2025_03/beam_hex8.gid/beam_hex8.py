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
## This file is generated on Mi 8. Feb 18:05:36 CET 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./beam_hex8.gid')
import beam_hex8_include
from beam_hex8_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################
def WriteNode(ifile, step, time, node):
    # dx = node.GetSolutionStepValue(DISPLACEMENT_X)
    # dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
    # dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
    # vx = node.GetSolutionStepValue(DISPLACEMENT_DT_X)
    # vy = node.GetSolutionStepValue(DISPLACEMENT_DT_Y)
    # vz = node.GetSolutionStepValue(DISPLACEMENT_DT_Z)
    # ax = node.GetSolutionStepValue(ACCELERATION_X)
    # ay = node.GetSolutionStepValue(ACCELERATION_Y)
    # az = node.GetSolutionStepValue(ACCELERATION_Z)
    dx = node.GetSolutionStepValue(DISPLACEMENT_EINS_X)
    dy = node.GetSolutionStepValue(DISPLACEMENT_EINS_Y)
    dz = node.GetSolutionStepValue(DISPLACEMENT_EINS_Z)
    vx = node.GetSolutionStepValue(DISPLACEMENT_EINS_DT_X)
    vy = node.GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y)
    vz = node.GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z)
    ax = node.GetSolutionStepValue(ACCELERATION_EINS_X)
    ay = node.GetSolutionStepValue(ACCELERATION_EINS_Y)
    az = node.GetSolutionStepValue(ACCELERATION_EINS_Z)
    ifile.write(str(step) + "\t" + str(time) + "\t" + str(dx) + "\t" + str(dy) + "\t" + str(dz) + "\t" + str(vx) + "\t" + str(vy) + "\t" + str(vz) + "\t" + str(ax) + "\t" + str(ay) + "\t" + str(az) + "\n")
##################################################################

def main(output=True, logging=True):
    model = beam_hex8_include.Model('beam_hex8',os.getcwd()+"/",os.getcwd()+"/", analysis_type=0, logging=logging)
    model.InitializeModel()

    ## perform static analysis
    tol = 1.0e-6

    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)
        if (abs(node.X0 - 10.0) < tol) and (abs(node.Y0 - 1.0) < tol):
            prescribed_nodes.append(node)

    for node in prescribed_nodes:
        node.Fix(DISPLACEMENT_Y)

    time = 0.0
    model.SolveModel(time)

    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.1)

    time = 1.0
    model.SolveModel(time)
    # model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     dx = node.GetSolutionStepValue(DISPLACEMENT_X)
    #     dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
    #     dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
    #     print(dy)
    # sys.exit(0)

    ####################################################

    ## perform dynamics analysis
    dynmodel = beam_hex8_include.Model('beam_hex8',os.getcwd()+"/",os.getcwd()+"/", analysis_type=2, logging=logging)
    dynmodel.InitializeModel()
    dynmodel.SetModelPart(model.model_part)

    model.model_part.Properties[1].SetValue(RAYLEIGH_DAMPING_ALPHA, 0.0 )
    model.model_part.Properties[1].SetValue(RAYLEIGH_DAMPING_BETA, 0.0 )

    for node in dynmodel.model_part.Nodes:
        node.Free(DISPLACEMENT_X)
        node.Free(DISPLACEMENT_Y)
        node.Free(DISPLACEMENT_Z)
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)

    for node in prescribed_nodes:
        dnode = dynmodel.model_part.Nodes[node.Id]
        # print("fix node: ", dnode.Id)
        dnode.Fix(DISPLACEMENT_Y)

    for node in dynmodel.model_part.Nodes:
        dx = node.GetSolutionStepValue(DISPLACEMENT_X)
        dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
        dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
        # print(dy)
        node.SetSolutionStepValue(DISPLACEMENT_NULL_X, dx)
        node.SetSolutionStepValue(DISPLACEMENT_EINS_X, dx)
        node.SetSolutionStepValue(DISPLACEMENT_NULL_Y, dy)
        node.SetSolutionStepValue(DISPLACEMENT_EINS_Y, dy)
        node.SetSolutionStepValue(DISPLACEMENT_NULL_Z, dz)
        node.SetSolutionStepValue(DISPLACEMENT_EINS_Z, dz)
        # node.SetSolutionStepValue(DISPLACEMENT_NULL_DT_X, 0.0)
        # node.SetSolutionStepValue(DISPLACEMENT_NULL_DT_Y, 0.0)
        # node.SetSolutionStepValue(DISPLACEMENT_NULL_DT_Z, 0.0)
        # node.SetSolutionStepValue(ACCELERATION_NULL_X, 0.0)
        # node.SetSolutionStepValue(ACCELERATION_NULL_Y, 0.0)
        # node.SetSolutionStepValue(ACCELERATION_NULL_Z, 0.0)

    if logging:
        tip_node = dynmodel.model_part.Nodes[42]
        ifile = open("displacement.txt", "w")
        ifile.write("# step\ttime\tdx\tdy\tdz\tvx\tvy\tvz\tax\tay\taz\n")

    step = 0
    time = 0.0

    if logging:
        WriteNode(ifile, step, time, tip_node)

    dynmodel.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    delta_time = 0.1
    step += 1
    time += delta_time
    dynmodel.SolveModel(time)
    if output:
        dynmodel.WriteOutput(time)
    if logging:
        WriteNode(ifile, step, time - delta_time, tip_node)

    dynmodel.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    for node in prescribed_nodes:
        dnode = dynmodel.model_part.Nodes[node.Id]
        # print("fix node: ", dnode.Id)
        dnode.Free(DISPLACEMENT_Y)

    nsteps = 1000
    # nsteps = 20
    for i in range(0, nsteps):
        step += 1
        time += delta_time
        dynmodel.SolveModel(time)
        if output:
            dynmodel.WriteOutput(time)
        if logging:
            WriteNode(ifile, step, time - delta_time, tip_node)

    return dynmodel

def test():
    dynmodel = main(logging=False, output=False)

    tip_node = dynmodel.model_part.Nodes[42]
    ref_values = [0.0038273908555991813, -0.018772871678550584, 0.0031342780018649716]

    assert(abs(tip_node.GetSolutionStepValue(DISPLACEMENT_EINS_Y) - ref_values[0]) < 1e-10)
    assert(abs(tip_node.GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y) - ref_values[1]) < 1e-10)
    assert(abs(tip_node.GetSolutionStepValue(ACCELERATION_EINS_Y) - ref_values[2]) < 1e-10)

    print("Test passed")

def tag():
    return "dynamics"

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
