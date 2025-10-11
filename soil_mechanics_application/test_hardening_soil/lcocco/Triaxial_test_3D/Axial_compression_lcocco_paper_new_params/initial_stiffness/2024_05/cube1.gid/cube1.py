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
## This file is generated on Fr 28. Jan 19:26:47 CET 2022
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./cube1.gid')
import cube1_include
from cube1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True, paxial = -450.0):

    model = cube1_include.Model('cube1',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)

    time = 0.0
    model.SolveModel(time)

    print("## isotropic compression ##")
    piso = -100.0
    nsteps = 100
    dp = piso / nsteps

    time = 0.0
    delta_time = abs(dp)
    for step in range(0, nsteps):
        p = (step+1)*dp
        time = time + delta_time

        print("###Isotropic loading step " + str(step+1) + ", p = " + str(p) + "###")

        for cond_id in model.layer_cond_sets['loadx']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)

        for cond_id in model.layer_cond_sets['loady']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)

        for cond_id in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        elem = model.model_part.Elements[1]
        strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
        # print("strain", strain)
        stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
        # print("stress", stress)
        last_ezz = strain[0][2]
        last_sxx = stress[0][0]
        last_syy = stress[0][1]
        last_szz = stress[0][2]

    # print(last_ezz)
    # sys.exit(0)
    print("## axial compression ##")

    if logging:
        ifile = open("results.txt", "w")
        ifile.write("epsilon_zz\tsigma_xx\tsigma_yy\tsigma_zz\n")
        ifile.write("0.0\t" + str(last_sxx) + "\t" + str(last_syy) + "\t" + str(last_szz) + "\n")
        ifile.flush()

    nsteps = 200
    dp = (paxial - piso) / nsteps
    time = 2*time
    delta_time = abs(dp)
    for step in range(0, nsteps):
        p = piso + (step+1)*dp
        time = time + delta_time

        print("###Axial loading step " + str(step+1) + ", p = " + str(p) + "###")

        for cond_id in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            elem = model.model_part.Elements[1]
            strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
            # print("strain", strain)
            stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
            # print("stress", stress)
            ifile.write(str(strain[0][2] - last_ezz) + "\t" + str(stress[0][0]) + "\t" + str(stress[0][1]) + "\t" + str(stress[0][2]) + "\n")
            ifile.flush()

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False, paxial = -300.0)

    ###### pytesting results
    last_ezz = -0.00530945443058
    ref_dezz = -6.724445205495796e-03
    elem = model.model_part.Elements[1]
    strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
    dezz = strain[0][2] - last_ezz
    print("dezz: %.15e" % dezz)
    assert(abs(dezz - ref_dezz) / abs(ref_dezz) < 1e-10)
    ########################

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
