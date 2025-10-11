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
## This file is generated on Mi 9. Mar 11:25:54 CET 2022 
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
def SetDisplacement(model_part, lmbda):
    tol = 1e-6
    for node in model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            node.SetSolutionStepValue(PRESCRIBED_DISPLACEMENT_X, 1.0e-1)
            node.SetSolutionStepValue(DISPLACEMENT_X, lmbda*1.0e-1)

def main(output=True, logging=True):

    model = square1_q9_include.Model('square1_q9',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1.0e-6
    load_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)

    model.solver.solver.set_displacement_callback = SetDisplacement
    model.solver.solver.time_scheme = ArcLengthDisplacementControlEnergyReleaseResidualBasedIncrementalUpdateStaticDeactivationScheme(ResidualBasedIncrementalUpdateStaticDeactivationScheme())

    if logging:
        ifile = open("results_damage.txt", "w")
        ifile.write("step\tdisp\tdamage\n")

        ifile2 = open("results_reaction.txt", "w")
        ifile2.write("step\tdisp\treaction\n")

    time = 0
    model.SolveModel(time)
    # sys.exit(0)

    if logging:
        ifile.write("0\t0.0\t0.0\n")
        ifile2.write("0\t0.0\t0.0\n")

    ## loading
    delta_time = 1.0
    step = 0

    step_list = []
    step_list.append([1, 1.0e-4])
    step_list.append([5, 1.0e-5])
    for tmp in step_list:
        nsteps = tmp[0]
        model.solver.solver.arc_length_control_process.GetConstraint().SetRadius(tmp[1])
        for i in range(0, nsteps):
            print("Load step " + str(step + 1) + " starts")
            time += delta_time
            model.SolveModel(time)
            if output:
                model.WriteOutput(time)
            print("Load step " + str(step + 1) + " completed")

            if logging:
                disp = model.model_part.Nodes[9].GetSolutionStepValue(DISPLACEMENT_X)
                damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
                ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
                ifile.flush()

                ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\n")
                ifile2.flush()

            step += 1

    # ## loading using energy release control
    # new_arc_length_control_process = ArcLengthControlProcess(ArcLengthDisplacementControlEnergyReleaseConstraint(1.0e-1))
    # model.solver.solver.SetArcLengthControlProcess(new_arc_length_control_process)

    # step_list = []
    # step_list.append([1, 1.0e-4])
    # for tmp in step_list:
    #     nsteps = tmp[0]
    #     model.solver.solver.arc_length_control_process.GetConstraint().SetRadius(tmp[1])
    #     for i in range(0, nsteps):
    #         print("Load step " + str(step + 1) + " starts")
    #         time += delta_time
    #         model.SolveModel(time)
    #         # model.WriteOutput(time)
    #         print("Load step " + str(step + 1) + " completed")

    #         disp = model.model_part.Nodes[9].GetSolutionStepValue(DISPLACEMENT_X)
    #         damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
    #         ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
    #         ifile.flush()

    #         ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\n")
    #         ifile2.flush()

    #         step += 1

    # ## unloading
    # lambda_current = model.solver.solver.arc_length_control_process.GetLambda()
    # print("Lambda at the end of loading: " + str(lambda_current))
    # delta_lambda = -lambda_current / 100.0
    # model.solver.solver.DeactivateArcLengthControl() # deactivate arc-length, load is user-controlled
    # for i in range(0, 50):
    #     lambda_current += delta_lambda
    #     SetLoad(model.model_part, lambda_current)
    #     model.solver.solver.arc_length_control_process.Update(delta_lambda)

    #     print("Unload step " + str(step + 1) + " starts")
    #     time += delta_time
    #     model.SolveModel(time)
    #     # model.WriteOutput(time)
    #     print("Unload step " + str(step + 1) + " completed")

    #     disp = model.model_part.Nodes[9].GetSolutionStepValue(DISPLACEMENT_X)
    #     damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
    #     ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
    #     ifile.flush()

    #     ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\n")
    #     ifile2.flush()

    #     step += 1

    # lambda_current = model.solver.solver.arc_length_control_process.GetLambda()
    # print("Lambda at the end of unloading: " + str(lambda_current))

    # ## reloading
    # step_list = []
    # step_list.append([100, 4.0e-6])
    # model.solver.solver.ActivateArcLengthControl()
    # first_reload = True
    # for tmp in step_list:
    #     nsteps = tmp[0]
    #     model.solver.solver.arc_length_control_process.SetRadius(tmp[1])
    #     for i in range(0, nsteps):
    #         if first_reload:
    #             model.solver.solver.arc_length_control_process.SetForcedForward(True)
    #             first_reload = False
    #         else:
    #             model.solver.solver.arc_length_control_process.SetForcedForward(False)

    #         print("Reload step " + str(step + 1) + " starts")
    #         time += delta_time
    #         model.SolveModel(time)
    #         # model.WriteOutput(time)
    #         print("Reload step " + str(step + 1) + " completed")

    #         disp = model.model_part.Nodes[9].GetSolutionStepValue(DISPLACEMENT_X)
    #         damage = model.model_part.Elements[1].CalculateOnIntegrationPoints(DAMAGE, model.model_part.ProcessInfo)
    #         ifile.write(str(step+1) + "\t" + str(disp) + "\t" + str(damage[0][0]) + "\n")
    #         ifile.flush()

    #         ifile2.write(str(step+1) + "\t" + str(disp) + "\t" + str(model.model_part.Nodes[1].GetSolutionStepValue(REACTION_X)) + "\n")
    #         ifile2.flush()

    #         step += 1

    # lambda_current = model.solver.solver.arc_length_control_process.GetLambda()
    # print("Lambda at the end of reloading: " + str(lambda_current))

    if logging:
        ifile.close()
        ifile2.close()

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
    ref_reac = 5.119083456470422e+01
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
