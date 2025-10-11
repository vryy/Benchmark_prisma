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

def main(logging=True, output=True, nsteps_phase = [150, 50, 150]):
    model = square1_q9_include.Model('square1_q9',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("results_damage.txt", "w")
        ifile.write("step\tdisp\tdamage\n")

        ifile2 = open("results_reaction.txt", "w")
        ifile2.write("step\tdisp\treaction\teqv_strain\n")

    time = 0
    model.SolveModel(time)

    if logging:
        ifile.write("0\t0.0\t0.0\n")
        ifile2.write("0\t0.0\t0.0\t0.0\n")

    ## set initial damage
    dv = [1e-2] * 9
    model.model_part.Elements[1].SetValuesOnIntegrationPoints(DAMAGE, dv, model.model_part.ProcessInfo)

    ## furthermore, we enforce initial elastic state
    for node in model.model_part.Nodes:
        node.Fix(EQUIVALENT_STRAIN)

    ## load the system elastically
    print("Elastic load step started")
    disp = 0.0
    initial_disp = 1e-5
    delta_disp = initial_disp
    delta_time = 10.0

    disp += delta_disp
    time += delta_time
    step = 0
    for node in prescribed_nodes:
        # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
        node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)
    print("Elastic load step completed")

    if logging:
        damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
        ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
        ifile.flush()

        ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(EQUIVALENT_STRAIN)) + "\n")

    # sys.exit(0)

    ## enable the damage again and set the initial equivalent strain
    ## here the damage must not be propagated
    ## in other words, the equivalent strain must stay zero
    print("Setting step started")
    for node in model.model_part.Nodes:
        node.Free(EQUIVALENT_STRAIN)

    es = model.model_part.Elements[1].GetValuesOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, model.model_part.ProcessInfo)
    new_es = [v[0] for v in es]
    model.model_part.Elements[1].SetValuesOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, new_es, model.model_part.ProcessInfo)

    step = 1
    delta_time = 1.0
    time += delta_time
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)
    print("Setting step completed")

    if logging:
        damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
        ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
        ifile.flush()

        ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(EQUIVALENT_STRAIN)) + "\n")

    # sys.exit(0)

    ## loading
    delta_disp = 1e-6
    delta_time = 1.0
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
            damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
            ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
            ifile.flush()

            ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\n")
            ifile2.flush()

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
            damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
            ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
            ifile.flush()

            ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(EQUIVALENT_STRAIN)) + "\n")
            ifile2.flush()

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
            damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
            ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
            ifile.flush()

            ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(EQUIVALENT_STRAIN)) + "\n")
            ifile2.flush()

    if logging:
        ifile.close()
        ifile2.close()

    return model

def test():
    model = main(logging=False, output=False)

    ###### pytesting results
    damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)[0][0]
    reac = model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)
    ref_damage = 0.8722151287422735
    ref_reac = -0.0686029471679605
    print("reac: %.15e, damage: %.15e" % (reac, damage))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    assert(abs(damage - ref_damage) / abs(ref_damage) < 1e-10)
    print("Test passed")
    ########################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, nsteps_phase = [150, 50, 150])

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
