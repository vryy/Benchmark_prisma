##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
#####     2023 by Hoang-Giang Bui                            #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mi 30. Aug 14:44:00 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./r2gtn_elpl.gid')
import r2gtn_elpl_include
from r2gtn_elpl_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, time, disp, prescribed_nodes, elem, process_info):
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    d = elem.CalculateOnIntegrationPoints(DAMAGE, process_info)
    p = elem.CalculateOnIntegrationPoints(PLASTICITY, process_info)
    ld = elem.CalculateOnIntegrationPoints(LOCAL_DAMAGE, process_info)
    lp = elem.CalculateOnIntegrationPoints(PLASTICITY_INDICATOR, process_info)
    ifile.write("%-*.6e%-*.10e%-*.10e%-*.10e%-*.10e%-*.10e%.10e\n" % (20, time, 20, disp, 20, reac_force_y, 20, ld[0][0], 20, lp[0][0], 20, d[0][0], p[0][0]))
    ifile.flush()

def main(output=True, logging=True):
    model = r2gtn_elpl_include.Model('r2gtn_elpl',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Y0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)
        # node.Fix(DAMAGE)
        # node.Fix(PLASTICITY)

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("%-*s%-*s%-*s%-*s%-*s%-*s%s\n" % (20, "time(us)", 20, "displacement", 20, "reaction", 20, "loc-d", 20, "loc-p", 20, "nonloc-d", "nonloc-p"))

    time = 0.0
    model.SolveModel(time)
    if logging:
        WriteLog(ifile, time*1e6, 0.0, prescribed_nodes, model.model_part.Elements[1], model.model_part.ProcessInfo)
    if output:
        model.WriteOutput(time)
    # sys.exit(0)

    total_disp = 2.0
    # total_disp = 1e-2

    delta_time = 1.0e-5 # 10 us
    delta_disp = 1e-2
    # delta_disp = 1.
    disp = 0.
    nsteps = int(total_disp / delta_disp)
    for i in range(0, nsteps):
        disp += delta_disp
        time += delta_time
        print("Load step disp %e starts" % (disp))
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
            # node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, delta_disp)
        model.SolveModel(time)
        if output:
            model.WriteOutput(disp)
        if logging:
            WriteLog(ifile, time*1e6, disp, prescribed_nodes, model.model_part.Elements[1], model.model_part.ProcessInfo)

    if (not output) and logging:
        model.WriteOutput(disp)

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    ### pytesting results
    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 1.0) < tol:
            prescribed_nodes.append(node)

    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    # ref_reac = 1.320684247667211e+02  # results obtained by Hereon laptop
    # ref_reac = 1.320680528025937e+02  # results obtained by AMD laptop
    ref_reac = 1.320679080759643e+02    # results obtained by Intel i9-13950HX laptop
    print("reac_force_y = %.15e, ref_reac = %.15e, diff = %e" % (reac_force_y, ref_reac, reac_force_y - ref_reac))
    assert(abs(reac_force_y - ref_reac) / ref_reac < 1e-10)
    print("Test passed")

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
