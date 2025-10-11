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
## This file is generated on Mi 27. Sep 12:16:41 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./solid_necking_bar_q4.gid')
import solid_necking_bar_q4_include
from solid_necking_bar_q4_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, prescribed_nodes, time, disp):
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ifile.write("%.6e\t%.10e\t%.10e\n" % (time, disp, reac_force_y))
    ifile.flush()

def main(total_disp=4.0,delta_disp=1e-2,output=True,logging=True,output_final=False):

    model = solid_necking_bar_q4_include.Model('solid_necking_bar_q4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    # ============================================ #
    # |       USER CALCULATION SCRIPT            | #
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv #

    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Y0 - 15.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("time(us)\t\tdisplacement\t\t\t\t\treaction\n")

    time = 0.0
    model.SolveModel(time)
    if logging:
        WriteLog(ifile, prescribed_nodes, time*1e6, 0.0)
    if output:
        model.WriteOutput(time)

    # total_disp = 4.0

    # delta_time = 1.0e-5 # 10 us
    # delta_disp = 5e-3
    # delta_disp = 1e-2
    # delta_disp = 1.
    disp = 0.
    nsteps = int(total_disp / delta_disp)
    delta_time = delta_disp
    for i in range(0, nsteps):
        disp += delta_disp
        time += delta_time
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, delta_disp)
        model.SolveModel(time)
        if logging:
            WriteLog(ifile, prescribed_nodes, time*1e3, disp)
        if output:
            model.WriteOutput(time*1e6)

    if output_final:
        model.WriteOutput(time*1e6)

    if logging:
        ifile.close()

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(total_disp=1.0,delta_disp=1e-2,output=False,logging=False)

    ######pytesting######
    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 15.0) < tol:
            prescribed_nodes.append(node)

    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ref_reac = 3.7729225631e+03
    diff = abs(reac_force_y - ref_reac) / abs(ref_reac)
    print("reac_force_y: %e, ref_reac: %e, diff: %e" % (reac_force_y, ref_reac, diff))
    assert(diff < 1e-10)
    #####################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=False,logging=True,output_final=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
