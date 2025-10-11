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
## This file is generated on Mo 28. Feb 10:33:42 CET 2022
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./square1_q9.gid')
import square1_q9_include
from square1_q9_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, step, disp, elem, reaction_nodes, process_info):
    damage = elem.CalculateOnIntegrationPoints(DAMAGE, process_info)
    stress = elem.CalculateOnIntegrationPoints(STRESSES, process_info)
    reac_x = 0.0
    for node in reaction_nodes:
        reac_x += node.GetSolutionStepValue(REACTION_X)
    ifile.write("%-*d%-*.10e%-*.10e%-*.10e%.10e\n" % (10, step, 20, disp, 20, damage[0][0], 20, reac_x, stress[0][0]))
    ifile.flush()

def main(logging=True, output=True, output_last=False, delta_disp = 1e-6, nsteps_phase = [150, 50, 150]):
    model = square1_q9_include.Model('square1_q9',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1.0e-6
    prescribed_nodes = []
    reaction_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            reaction_nodes.append(node)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("results.txt", "w")
        ifile.write("%-*s%-*s%-*s%-*s%s\n" % (10, "step", 20, "disp", 20, "damage", 20, "reaction", "stress_xx"))

    time = 0
    model.SolveModel(time)

    if logging:
        ifile.write("%-*d%-*.10e%-*.10e%-*.10e%.10e\n" % (10, 0, 20, 0.0, 20, 0.0, 20, 0.0, 0.0))

    ## loading
    delta_time = 1.0
    disp = 0.0
    for step in range(0, nsteps_phase[0]):
        print("Load step " + str(step + 1) + " starts")
        disp += delta_disp
        time += delta_time
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        print("Load step " + str(step + 1) + " completed")

        if logging:
            WriteLog(ifile, step, disp, model.model_part.Elements[1], reaction_nodes, model.model_part.ProcessInfo)

    ## unloading
    delta_disp = -delta_disp
    for step in range(0, nsteps_phase[1]):
        print("Unload step " + str(step + 1) + " starts")
        disp += delta_disp
        time += delta_time
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        print("Unload step " + str(step + 1) + " completed")

        if logging:
            WriteLog(ifile, step, disp, model.model_part.Elements[1], reaction_nodes, model.model_part.ProcessInfo)

    ## reloading
    delta_disp = -delta_disp
    for step in range(0, nsteps_phase[2]):
        print("Reload step " + str(step + 1) + " starts")
        disp += delta_disp
        time += delta_time
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        print("Reload step " + str(step + 1) + " completed")

        if logging:
            WriteLog(ifile, step, disp, model.model_part.Elements[1], reaction_nodes, model.model_part.ProcessInfo)

    if logging:
        ifile.close()

    if output_last:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1.0e-6
    reaction_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            reaction_nodes.append(node)

    ###### pytesting results
    damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)[0][0]
    reac = 0.0
    for node in reaction_nodes:
        reac += node.GetSolutionStepValue(REACTION_X)
    ref_damage = 9.846895143068409e-01
    ref_reac = -8.590778640180035e-01
    print("reac: %.15e, damage: %.15e" % (reac, damage))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    assert(abs(damage - ref_damage) / abs(ref_damage) < 1e-10)
    print("Test passed")
    ########################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, output_last=True, nsteps_phase = [150, 50, 150])

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
