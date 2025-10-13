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
## This file is generated on Mi 17. Jun 17:23:39 CEST 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./two_beams.gid')
import two_beams_include
from two_beams_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

B1 = 1.0
B2 = 0.1
H = 0.7
E = 2.0e7
A = 1.0
rho = 1.0
L1 = math.sqrt(H**2 + B1**2)
alpha = 100.0
beta = 1803.0 #100.0

def main(output=True, logging=True):
    model = two_beams_include.Model('two_beams',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    model.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model.model_part.CreateNewNode(2, -B2, H, 0.0)
    model.model_part.CreateNewNode(3, B1, H, 0.0)
    model.AddDofs(model.model_part)

    prop1 = model.model_part.Properties[1]
    prop1.SetValue(YOUNG_MODULUS, E)
    prop1.SetValue(AREA, A)
    prop1.SetValue(DENSITY, rho)
    prop1.SetValue(RAYLEIGH_DAMPING_ALPHA, 0.0)
    prop1.SetValue(RAYLEIGH_DAMPING_BETA, 0.0)

    prop2 = model.model_part.Properties[2]
    prop2.SetValue(YOUNG_MODULUS, alpha*E)
    prop2.SetValue(AREA, A)
    prop2.SetValue(DENSITY, beta*rho)
    prop2.SetValue(RAYLEIGH_DAMPING_ALPHA, 0.0)
    prop2.SetValue(RAYLEIGH_DAMPING_BETA, 0.0)

    model.model_part.CreateNewElement("TrussElement3D2N", 1, [1, 3], prop1)
    model.model_part.CreateNewElement("TrussElement3D2N", 2, [1, 2], prop2)

    print(model.model_part)

    model.model_part.Nodes[1].Fix(DISPLACEMENT_Y)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_X)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_Y)
    model.model_part.Nodes[3].Fix(DISPLACEMENT_X)

    model.model_part.Nodes[1].Fix(DISPLACEMENT_Z)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_Z)
    model.model_part.Nodes[3].Fix(DISPLACEMENT_Z)

    model.model_part.Elements[1].Initialize(model.model_part.ProcessInfo)
    model.model_part.Elements[2].Initialize(model.model_part.ProcessInfo)

    ux0 = 0.0183311538269 # obtain from static analysis
    uy0 = 0.160902425421

    vx0 = 0.0
    vy0 = 0.0

    ax0 = 0.0
    ay0 = -1.783983985694714e6 # u0dd = -inv(M)*K*u0

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_NULL_X, ux0)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_NULL_Y, uy0)

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_NULL_DT_X, vx0)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_NULL_DT_Y, vy0)

    model.model_part.Nodes[1].SetSolutionStepValue(ACCELERATION_NULL_X, ax0)
    model.model_part.Nodes[3].SetSolutionStepValue(ACCELERATION_NULL_Y, ay0)

    if logging:
        ifile = open("output.txt", "w")
        ifile.write("time\tux/b2\tuy/h\tvx/b2\tvy/h\tax/b2\tay/h\n")

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    time = 1.0e-4
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    if logging:
        ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
        uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
        vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_DT_X)
        vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_DT_Y)
        ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_X)
        ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_Y)
        ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(vx/B2) + '\t' + str(vy/H) + '\t' + str(ax/B2) + '\t' + str(ay/H) + '\n')

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    delta_time = 1.0e-4
    nsteps = 1000
    for i in range(0, nsteps):
        time = time + delta_time

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
            uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
            vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_DT_X)
            vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_DT_Y)
            ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_X)
            ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_Y)
            ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(vx/B2) + '\t' + str(vy/H) + '\t' + str(ax/B2) + '\t' + str(ay/H) + '\n')

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
    uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
    vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_DT_X)
    vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_DT_Y)
    ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_X)
    ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_Y)

    ref_values = [0.05588869130767097, -0.06443579099413001, 70.25053257948413, -595.1953617115364, -17039.123482617855, 1004344.1300272287]
    assert(abs(ux/B2 - ref_values[0]) < 1e-10)
    assert(abs(uy/H - ref_values[1]) < 1e-10)
    assert(abs(vx/B2 - ref_values[2]) < 1e-10)
    assert(abs(vy/H - ref_values[3]) < 1e-10)
    assert(abs(ax/B2 - ref_values[4]) < 1e-10)
    assert(abs(ay/H - ref_values[5]) < 1e-10)
    print("Test passed")

def tag():
    return "dynamics"

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
