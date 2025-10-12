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
import two_beams_include
from two_beams_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True, max_time=0.1, delta_time=1e-4):
    model = two_beams_include.Model('two_beams',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    B1 = 1.0
    B2 = 0.1
    H = 0.7
    E = 2.0e7
    A = 1.0
    rho = 1.0
    L1 = math.sqrt(H**2 + B1**2)
    alpha = 100.0
    # beta = 1803.0 #100.0
    beta = 100.0 # 1803.0

    model.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model.model_part.CreateNewNode(2, -B2, H, 0.0)
    model.model_part.CreateNewNode(3, B1, H, 0.0)
    model.AddDofs(model.model_part)

    prop1 = model.model_part.Properties[1]
    prop1.SetValue(YOUNG_MODULUS, E)
    prop1.SetValue(AREA, A)
    prop1.SetValue(DENSITY, rho)
    prop1.SetValue(RAYLEIGH_DAMPING_ALPHA, 100.0)
    prop1.SetValue(RAYLEIGH_DAMPING_BETA, 0.0)

    prop2 = model.model_part.Properties[2]
    prop2.SetValue(YOUNG_MODULUS, alpha*E)
    prop2.SetValue(AREA, A)
    prop2.SetValue(DENSITY, beta*rho)
    prop2.SetValue(RAYLEIGH_DAMPING_ALPHA, 100.0)
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

    time = 0.0

    ux0 = 0.0183311538269 # obtain from static analysis
    uy0 = 0.160902425421

    vx0 = 0.0
    vy0 = 0.0

    ax0 = 0.0
    ay0 = -1.783983985694714e6 # u0dd = -inv(M)*K*u0

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_NULL_X, ux0 - vx0*delta_time + ax0*delta_time*delta_time/2)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_NULL_Y, uy0 - vy0*delta_time + ay0*delta_time*delta_time/2)

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_NULL_DT_X, 0.0)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_NULL_DT_Y, 0.0)

    model.model_part.Nodes[1].SetSolutionStepValue(ACCELERATION_NULL_X, 0.0)
    model.model_part.Nodes[3].SetSolutionStepValue(ACCELERATION_NULL_Y, 0.0)

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_X, ux0 - vx0*delta_time + ax0*delta_time*delta_time/2)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_Y, uy0 - vy0*delta_time + ay0*delta_time*delta_time/2)

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_DT_X, 0.0)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_DT_Y, 0.0)

    model.model_part.Nodes[1].SetSolutionStepValue(ACCELERATION_X, 0.0)
    model.model_part.Nodes[3].SetSolutionStepValue(ACCELERATION_Y, 0.0)

    if logging:
        ifile = open("output.txt", "w")
        ifile.write("time\tux/b2\tuy/h\tvx_n/b2\tvy_n/h\tax_n/b2\tay_n/h\n")

        ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_X)
        uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_Y)
        vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
        vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
        ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_NULL_X)
        ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_NULL_Y)
        ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(vx/B2) + '\t' + str(vy/H) + '\t' + str(ax/B2) + '\t' + str(ay/H) + '\n')

    time += delta_time

    model.model_part.CloneTimeStep(time)

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_NULL_X, ux0)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_NULL_Y, uy0)

    model.model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_NULL_DT_X, 0.0)
    model.model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_NULL_DT_Y, 0.0)

    model.model_part.Nodes[1].SetSolutionStepValue(ACCELERATION_NULL_X, 0.0)
    model.model_part.Nodes[3].SetSolutionStepValue(ACCELERATION_NULL_Y, 0.0)

    if logging:
        ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_X)
        uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_Y)
        vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
        vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
        ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_NULL_X)
        ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_NULL_Y)
        ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(vx/B2) + '\t' + str(vy/H) + '\t' + str(ax/B2) + '\t' + str(ay/H) + '\n')


    # model.model_part.CloneTimeStep(time + delta_time)

    # print("displacement of node 3")
    # print(model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_Y))
    # print(model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_Y, 1))
    # print(model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_Y, 2))
    # sys.exit(0)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    time += delta_time
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    if logging:
        # ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
        # uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
        # vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_DT_X)
        # vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_DT_Y)
        # ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_X)
        # ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_Y)
        ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_EINS_X)
        uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_EINS_Y)
        vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
        vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
        ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_NULL_X)
        ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_NULL_Y)
        ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(vx/B2) + '\t' + str(vy/H) + '\t' + str(ax/B2) + '\t' + str(ay/H) + '\n')

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    nsteps = int(max_time/delta_time)
    for i in range(0, nsteps):
        time = time + delta_time

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            # ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
            # uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
            # vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_DT_X)
            # vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_DT_Y)
            # ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_X)
            # ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_Y)
            ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_EINS_X)
            uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_EINS_Y)
            vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
            vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
            ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_NULL_X)
            ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_NULL_Y)
            ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(vx/B2) + '\t' + str(vy/H) + '\t' + str(ax/B2) + '\t' + str(ay/H) + '\n')
            ifile.flush()

    if logging:
        ifile.close()

    return model

def test():
    model = main(output=False, logging=False, max_time=0.01)

    ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_EINS_X)
    uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_EINS_Y)
    vx = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
    vy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
    ax = model.model_part.Nodes[1].GetSolutionStepValue(ACCELERATION_NULL_X)
    ay = model.model_part.Nodes[3].GetSolutionStepValue(ACCELERATION_NULL_Y)

    print("ux: %.16e" % ux)
    print("uy: %.16e" % uy)
    print("vx: %.16e" % vx)
    print("vy: %.16e" % vy)
    print("ax: %.16e" % ax)
    print("ay: %.16e" % ay)
    ref_ux = -1.3561024000023994e-02
    ref_uy = 5.4664759447915524e-02
    ref_vx = -1.4814689811476180e+01
    ref_vy = 7.0643628877126204e+01
    ref_ax = 5.3046306217050267e+04
    ref_ay = -9.3271732478672848e+05
    assert(abs(ux - ref_ux) < 1e-10)
    assert(abs(uy - ref_uy) < 1e-10)
    assert(abs(vx - ref_vx) / abs(ref_vx) < 1e-10)
    assert(abs(vy - ref_vy) / abs(ref_vy) < 1e-10)
    assert(abs(ax - ref_ax) / abs(ref_ax) < 1e-10)
    assert(abs(ay - ref_ay) / abs(ref_ay) < 1e-10)

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
