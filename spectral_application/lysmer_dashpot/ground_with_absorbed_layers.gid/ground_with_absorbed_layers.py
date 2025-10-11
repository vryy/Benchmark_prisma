##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025-2026 by Hoang-Giang Bui (UoB)            #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Fri May 16 12:13:19 PM BST 2025
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import ground_with_absorbed_layers_include as simulation_include
from ground_with_absorbed_layers_include import *
model_name_ = 'ground_with_absorbed_layers'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True, nsteps=500, delta_time=0.025):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging)
    model.InitializeModel()

    ## add a point load to generate Ricker wave
    ## This simulation is adopted from Chau et al, Hybrid asynchronous isogeometric Perfectly Matched Layer for transient elastodynamics, Computers and Geotechnics, 2023.
    tol = 1e-06
    load_node = None
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) and (abs(node.Y0) < tol):
            load_node = node
            break

    if load_node == None:
        raise Exception("Can't find load node")

    last_cond_id = model.model_part.GetLastConditionId()
    model.model_part.CreateNewCondition("PointForce3D", last_cond_id+1, [load_node.Id], model.model_part.Properties[1])

    ## boundary condition
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol):
            node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Z)

    ## time stepping
    time = 0.0

    A = 1e6
    ts = 3.0
    tp = 3.0

    disp_bounds = [-1e99, 1e99]

    if logging:
        pointC = None
        pointD = None

        for node in model.model_part.Nodes:
            if (abs(node.X0-100.0) < tol) and (abs(node.Y0-0.0) < tol):
                pointC = node
            if (abs(node.X0-100.0) < tol) and (abs(node.Y0+250.0) < tol):
                pointD = node

        if pointC == None:
            raise Exception("Point C cannot be found")

        if pointD == None:
            raise Exception("Point D cannot be found")

        ifile = open("monitoring.log", "w")
        ifile.write("%-*s%-*s%-*s%-*s%s\n" % (20, "time", 20, "C-ux", 20, "C-uy", 20, "D-ux", "D-uy"))

    for step in range(0, nsteps):

        time += delta_time

        aux = math.pi*(time-ts)/tp
        force = A * (2*aux**2 - 1.0) * math.exp(-aux**2)

        load_node.SetSolutionStepValue(FORCE_Y, force)

        model.SolveModel(time)

        for node in model.model_part.Nodes:
            disp = node.GetSolutionStepValue(DISPLACEMENT)
            if (disp[1] > disp_bounds[0]):
                disp_bounds[0] = disp[1]
            if (disp[1] < disp_bounds[1]):
                disp_bounds[1] = disp[1]

        if output:
            model.WriteOutput(time)

        if logging:
            Cux = pointC.GetSolutionStepValue(DISPLACEMENT_X)
            Cuy = pointC.GetSolutionStepValue(DISPLACEMENT_Y)
            Dux = pointD.GetSolutionStepValue(DISPLACEMENT_X)
            Duy = pointD.GetSolutionStepValue(DISPLACEMENT_Y)
            ifile.write("%-*.10e%-*.10e%-*.10e%-*.10e%.10e\n" % (20, time, 20, Cux, 20, Cuy, 20, Dux, Duy))

    if logging:
        ifile.close()

    print("disp_bounds:", disp_bounds)

    return model

def test():
    model = main(logging=False, output=False, nsteps=400)

    pointC = None

    tol = 1e-06
    for node in model.model_part.Nodes:
        if (abs(node.X0-100.0) < tol) and (abs(node.Y0-0.0) < tol):
            pointC = node

    ux = pointC.GetSolutionStepValue(DISPLACEMENT_X)
    uy = pointC.GetSolutionStepValue(DISPLACEMENT_Y)
    # print("%.16e" % ux)
    # print("%.16e" % uy)
    ref_ux = 1.3773257410783072e-03
    ref_uy = 1.7926657836923360e-03
    assert(abs(ux - ref_ux) < 1e-10)
    assert(abs(uy - ref_uy) < 1e-10)
    print("Test passed")

def tag():
    return "Lysmer,ABC"

def print_tag():
    print("Tag(s): " + tag())

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
