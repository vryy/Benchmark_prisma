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
# calculate insitu-stress for geology_virgin.gid
model = two_beams_include.Model('two_beams',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()
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
beta = 100.0
u0 = -1.0

def SetDisplacement(model_part, factor):
    model_part.Nodes[3].SetSolutionStepValue(PRESCRIBED_DISPLACEMENT_X, 0.0)
    model_part.Nodes[3].SetSolutionStepValue(PRESCRIBED_DISPLACEMENT_Y, u0)
    model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_X, 0.0)
    model_part.Nodes[3].SetSolutionStepValue(DISPLACEMENT_Y, factor*u0)

def main(output=True, logging=True):
    model.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model.model_part.CreateNewNode(2, -B2, H, 0.0)
    model.model_part.CreateNewNode(3, B1, H, 0.0)
    model.AddDofs(model.model_part)

    prop1 = model.model_part.Properties[1]
    prop1.SetValue(YOUNG_MODULUS, E)
    prop1.SetValue(CROSS_AREA, A)
    prop1.SetValue(DENSITY, rho)

    prop2 = model.model_part.Properties[2]
    prop2.SetValue(YOUNG_MODULUS, alpha*E)
    prop2.SetValue(CROSS_AREA, A)
    prop2.SetValue(DENSITY, beta*rho)

    model.model_part.CreateNewElement("CrisfieldTrussElement3D2N", 1, [1, 3], prop1)
    model.model_part.CreateNewElement("CrisfieldTrussElement3D2N", 2, [1, 2], prop2)
    model.model_part.CreateNewCondition("PointForce2D", 1, [3], prop1)

    if output:
        print(model.model_part)

    model.model_part.Nodes[1].Fix(DISPLACEMENT_Y)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_X)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_Y)
    model.model_part.Nodes[3].Fix(DISPLACEMENT_X)
    model.model_part.Nodes[3].Fix(DISPLACEMENT_Y)

    model.model_part.Nodes[1].Fix(DISPLACEMENT_Z)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_Z)
    model.model_part.Nodes[3].Fix(DISPLACEMENT_Z)

    model.model_part.Elements[1].Initialize(model.model_part.ProcessInfo)
    model.model_part.Elements[2].Initialize(model.model_part.ProcessInfo)

    model.solver.solver.set_displacement_callback = SetDisplacement

    time = 0.0
    if output:
        ifile = open("output.txt", "w")
        ifile.write("time\tux/b2\tuy/h\tux\tuy\tlambda\n")
        ifile.write(str(time) + '\t0.0\t0.0\n')

    delta_time = 1.0
    nsteps = 300
    for i in range(0, nsteps):
        time = time + delta_time

        model.SolveModel(time)
        # model.WriteOutput(time)

        if output:
            ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
            uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
            print("u step: " + str(ux) + " " + str(uy))
            lmbda = model.solver.solver.arc_length_control_process.GetLambda()
            print("lambda: " + str(lmbda))
            ifile.write(str(time) + '\t' + str(ux/B2) + '\t' + str(uy/H) + '\t' + str(ux) + '\t' + str(uy) + "\t" + str(lmbda) + '\n')

    if output:
        ifile.close()

    ######### pytesting results #########
    ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
    uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
    ref_ux = 0.768382678744
    ref_uy = -3.36209118097
    assert(abs(ux/B2 - ref_ux) < 1e-11)
    assert(abs(uy/H - ref_uy) < 1e-11)
    #####################################

def test():
    main(False, False)

if __name__ == '__main__':
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
