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

def WriteLog(ifile, time, disp, detF, prescribed_nodes):
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ifile.write("%.6e\t%.10e\t%.10e\t%10e\n" % (time, disp, reac_force_y, detF))
    ifile.flush()

def main(output=True, logging=True):
    model = quad9_include.Model('quad9',os.getcwd()+"/",os.getcwd()+"/",logging)
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

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("time(us)\t\tdisplacement\t\t\t\t\treaction\t\tdetF\n")

    time = 0.0
    model.SolveModel(time)
    if logging:
        WriteLog(ifile, time*1e6, 0.0, 1.0, prescribed_nodes)
    if output:
        model.WriteOutput(time)

    # total_disp = 0.75
    total_disp = 2.0
    # total_disp = 1.0

    delta_disp = 1e-2
    delta_time = delta_disp
    # delta_disp = 1.
    disp = 0.
    nsteps = int(total_disp / delta_disp)
    for i in range(0, nsteps):
        disp += delta_disp
        time += delta_time
        print("********************Solve step %d, displacement %f" % (i, disp))
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, delta_disp)
        model.SolveModel(time)
        if output:
            model.WriteOutput(time*1e6)
        if logging:
            detF = model.model_part.Elements[1].CalculateOnIntegrationPoints(CURRENT_DEFORMATION_GRADIENT_DETERMINANT, model.model_part.ProcessInfo)
            WriteLog(ifile, time*1e6, disp, detF[0][0], prescribed_nodes)

    if logging:
        ifile.close()

    return model

def test():
    model = main(output=False, logging=False)

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 1.0) < tol:
            prescribed_nodes.append(node)

    ######### pytesting results #########
    # ref_reac = 2.8766917611e+04
    ref_reac = 2.876691764946716e+04
    reac = 0.0
    for node in prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)
    error = abs(reac - ref_reac) / abs(ref_reac)
    print("reac: %.15e, ref_reac: %.15e, error: %.10e" % (reac, ref_reac, error))
    assert(error < 1e-7)
    #####################################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=False, logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
