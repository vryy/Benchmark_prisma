##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mi 4. Aug 11:38:04 CEST 2021 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh_hyplas_q4.gid')
import mesh_hyplas_q4_include
from mesh_hyplas_q4_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, disp, nodes):
    reac = 0.0
    for node in nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)
    ifile.write("%.10e\t%.10e\n" % (disp, reac))
    ifile.flush()

def main(logging=True,output=True):
    model = mesh_hyplas_q4_include.Model('mesh_hyplas_q4',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    ## boundary condition
    ymin = 0.0
    ymax = 2.666700E+01
    xmin = 0.0
    tol = 1.0e-6

    prescribed_nodes = []

    for node in model.model_part.Nodes:
        if abs(node.X0 - xmin) < tol:
            node.Fix(DISPLACEMENT_X)

        if abs(node.Y0 - ymin) < tol:
            node.Fix(DISPLACEMENT_Y)

        if abs(node.Y0 - ymax) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    # print("prescribed_nodes:")
    # for node in prescribed_nodes:
    #     print(node.Id)
    # sys.exit(0)

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("disp\treaction\n")

    ## load increment - displacement control

    time = 0.0
    disp = 0.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)
    if logging:
        WriteLog(ifile, disp, prescribed_nodes)

    delta_disp_list = []
    # for i in range(0, 10):
    #     delta_disp_list.append(0.5)

    for i in range(0, 50):
        delta_disp_list.append(0.1)

    for du in delta_disp_list:
        disp += du
        for node in prescribed_nodes:
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, du)
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

        delta_time = du
        time = time + delta_time
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            WriteLog(ifile, disp, prescribed_nodes)

    if logging:
        ifile.close()
    print("Analysis completed")

    ######### pytesting results #########
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    print(reac_force_y)
    ref_reac = 4.34311102815 # this is the value with one thread. The obtained value will change with number of threads.
    test = abs(reac_force_y - ref_reac) / (abs(reac_force_y) + abs(ref_reac))
    print(test)
    assert(test < 1e-7)
    #####################################

def test():
    main(False, False)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(True, True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
