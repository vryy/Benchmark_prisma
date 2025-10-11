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
## This file is generated on Di 19. Sep 13:01:59 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import quad4_include
from quad4_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, time, disp, prescribed_nodes):
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ifile.write("%.6e\t%.10e\t%.10e\n" % (time, disp, reac_force_y))
    ifile.flush()

def main(logging=True,output=True):
    model = quad4_include.Model('quad4',current_dir_,current_dir_,logging)
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
        ifile.write("time(us)\t\tdisplacement\t\t\t\t\treaction\n")

    time = 0.0
    model.SolveModel(time)
    if logging:
        WriteLog(ifile, time*1e6, 0.0, prescribed_nodes)
    if output:
        model.WriteOutput(time)

    total_disp = 2.0

    delta_time = 1.0e-5 # 10 us
    delta_disp = 1e-2
    # delta_disp = 1.
    disp = 0.
    nsteps = int(total_disp / delta_disp)
    for i in range(0, nsteps):
        disp += delta_disp
        time += delta_time
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, delta_disp)
        model.SolveModel(time)
        if logging:
            WriteLog(ifile, time*1e6, disp, prescribed_nodes)
        if output:
            model.WriteOutput(time*1e6)

    if logging:
        ifile.close()

    ######### pytesting results #########
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ref_reac = 2.8767014628e+04 # this result differ to Abaqus but is the same as with Hyplas
    test = abs(reac_force_y - ref_reac) / (abs(reac_force_y) + abs(ref_reac))
    assert(test < 1e-11)
    #####################################

def test():
    main(logging=False, output = False)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output = False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
