##################################################################
## This file is generated on Di 19. Sep 13:01:59 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
## Reproducing the test from
##   Griffiths et al, Observations on the extended Matsuoka–Nakai failure criterion
##################################################################
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

def WriteLog2(ifile, time, disp, element, process_info):
    stress = element.CalculateOnIntegrationPoints(STRESSES, process_info)
    ifile.write("%.6e\t%.10e\t%.10e\n" % (time, disp, stress[0][1]))
    ifile.flush()

def main(output=True, logging=True):
    model = quad4_include.Model('quad4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Y0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("monitoring.log", "w")
        # ifile.write("time(us)\t\tdisplacement\t\t\t\t\treaction\n")
        ifile.write("time(us)\t\te1\t\t\t\t\to1\n")

    time = 0.0
    model.SolveModel(time)
    if logging:
        # WriteLog(ifile, time*1e6, 0.0, prescribed_nodes)
        WriteLog2(ifile, time*1e6, 0.0, model.model_part.Elements[1], model.model_part.ProcessInfo)
    if output:
        model.WriteOutput(time)

    print("*******LOADING STARTED**********")

    disp = 0.0
    time = 0.0
    delta_disp = 1e-5
    delta_time = delta_disp
    nsteps = 80
    for i in range(0, nsteps):
        disp -= delta_disp
        time += delta_time
        print("*******LOAD STEP disp = " + str(disp))
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, -delta_disp)
        model.SolveModel(time)
        if output:
            model.WriteOutput(time*1e6)
        if logging:
            # WriteLog(ifile, time*1e6, disp, prescribed_nodes)
            WriteLog2(ifile, time*1e6, disp, model.model_part.Elements[1], model.model_part.ProcessInfo)

    if logging:
        ifile.close()

    return model

def tag():
    return "matsuoka-nakai"

def print_tag():
    print("Tag(s): " + tag())

def test():
    model = main(output=False, logging=False)

    ######### pytesting results #########
    stress = model.model_part.Elements[1].CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
    sigma_1 = stress[0][1]
    ref_sigma_1 = -4.2844440202e+04
    print("%.10e" % (sigma_1))
    assert(abs(sigma_1 - ref_sigma_1) / abs(ref_sigma_1) < 1e-10)
    #####################################
    print("Test passed")

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
