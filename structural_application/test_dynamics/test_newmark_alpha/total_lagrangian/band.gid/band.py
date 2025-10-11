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
## This file is generated on Do 15. Apr 19:55:16 CEST 2021
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./band.gid')
import band_include
from band_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True, delta_time = 1e-3, max_time = 5.0):
    tol = 1.0e-6

    # create and solve a static model
    model_static = band_include.Model('band',os.getcwd()+"/",os.getcwd()+"/",0,0.0,0.0,0.0, logging=logging)
    model_static.InitializeModel()
    print("Static model is initialized")

    for node in model_static.model_part.Nodes:
        if abs(node.Y0 - 0.21) < tol:
            node.SetSolutionStepValue(FACE_LOAD_Y, 1.0e1)
        if abs(node.X0 - 0.24899) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

    model_static.SolveModel(0.0)
    for node in model_static.model_part.Nodes:
        print("node " + str(node.Id) + " DISPLACEMENT: " + str(node.GetSolutionStepValue(DISPLACEMENT)))

    # solve the dynamic model

    model = band_include.Model('band',os.getcwd()+"/",os.getcwd()+"/",2,1.0,0.0,0.0, logging=logging)
    model.InitializeModel()

    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.24899) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

        snode = model_static.model_part.Nodes[node.Id]
        node.SetSolutionStepValue(DISPLACEMENT_NULL_X, snode.GetSolutionStepValue(DISPLACEMENT_X))
        node.SetSolutionStepValue(DISPLACEMENT_NULL_Y, snode.GetSolutionStepValue(DISPLACEMENT_Y))

        if abs(node.X0 - 0.6) < tol and abs(node.Y0 - 0.21) < tol:
            probe_node = node

    if logging:
        ifile = open("tip_deflection.txt", "w")
        ifile.write("t\tdy\n")

    time = 0.0
    nsteps = int(max_time/delta_time)

    for step in range(0, nsteps):
        time = time + delta_time
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            ifile.write(str(time)+"\t"+str(probe_node.GetSolutionStepValue(DISPLACEMENT_Y))+"\n")

    if logging:
        ifile.close()

    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False, max_time=1.0)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.6) < tol and abs(node.Y0 - 0.21) < tol:
            probe_node = node

    ref_disp_y = -3.0600858324625782e-03
    assert(abs(probe_node.GetSolutionStepValue(DISPLACEMENT_Y) - ref_disp_y) < 1e-10)

    print("Test passed")

def tag():
    return "dynamics,neo-hookean"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
