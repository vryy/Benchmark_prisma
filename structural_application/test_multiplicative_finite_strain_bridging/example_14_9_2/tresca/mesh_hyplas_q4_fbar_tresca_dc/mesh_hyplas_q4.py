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
    ifile.write("%.10e\t%.15e\n" % (disp, reac))
    ifile.flush()

def main(output=True, logging=True):
    model = mesh_hyplas_q4_include.Model('mesh_hyplas_q4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## boundary condition
    ymin = 0.0
    ymax = 2.667000E+01
    xmin = 0.0
    xmax = 6.413
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
    # model.WriteOutput(time)
    if logging:
        WriteLog(ifile, disp, prescribed_nodes)

    delta_disp_list = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.125, 0.125, 0.25, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.125, 0.125, 0.25]

    print("*********LOADING STARTED**********")

    for du in delta_disp_list:
        disp += du
        print("*********LOAD STEP " + str(disp) + " STARTED")
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

        print("Displacement")
        for node in model.model_part.Nodes:
            print("%d  %.16e   %.16e" % (node.Id, node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y)))

    if logging:
        ifile.close()
    print("Analysis completed")

    return model

def test():
    model = main(output=False, logging=False)

    tol = 1.0e-6
    ymax = 2.667000E+01
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - ymax) < tol:
            prescribed_nodes.append(node)

    ######### pytesting results #########
    ref_reac = 2.2276955901865147e+01
    reac = 0.0
    for node in prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)
    print("reac: %.16e" % (reac))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    #####################################
    print("Test passed")

if __name__ == "__main__":
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
