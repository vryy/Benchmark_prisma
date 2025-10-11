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
import cube1_include
from cube1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):

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
    # sys.exit(0)

    piso = -100.0 # final value of isotropic loading state

    # setting the prestress
    # values = [[-piso, -piso, -piso, 0.0, 0.0, 0.0]]*8
    values = [[-piso, -piso, -piso, 0.0, 0.0, 0.0]]*1
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(PRESTRESS, values, 6, model.model_part.ProcessInfo)
        elem.ResetConstitutiveLaw()
    stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
    print("stress", stress)
    # sys.exit(0)

    print("## isotropic compression ##")
    nsteps = 100
    dp = piso / nsteps

    time = 0.0
    delta_time = abs(dp)
    for step in range(0, 1):
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
        print("strain", strain)
        stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
        print("stress", stress)
        last_ezz = strain[0][2]
        last_sxx = stress[0][0]
        last_syy = stress[0][1]
        last_szz = stress[0][2]

        # ifile.write(str(strain[0][2]) + "\t" + str(stress[0][0]) + "\t" + str(stress[0][1]) + "\t" + str(stress[0][2]) + "\n")

    # sys.exit(0)

    print("last_ezz:", last_ezz)
    print("last_sxx:", last_sxx)
    print("last_syy:", last_syy)
    print("last_szz:", last_szz)
    # print("last_Eur:", last_Eur)
    for node in model.model_part.Nodes:
        print(node.GetSolutionStepValue(DISPLACEMENT))
    # sys.exit(0)

    if logging:
        ifile = open("results.txt", "w")
        ifile.write("epsilon_zz\tsigma_xx\tsigma_yy\tsigma_zz\tp\tdp\n")
        # ifile.write("0.0\t0.0\t0.0\t0.0\n")

    print("## axial compression ##")

    if logging:
        ifile.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (0.0, last_sxx, last_syy, last_szz, piso, 0.0))
        ifile.flush()

    nsteps = 200
    paxial = -500.0 # -450.0
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

        converged = model.SolveModel(time)
        print("converged", converged)
        if not converged:
            raise Exception("The constitutive law does not converge")
        if output:
            model.WriteOutput(time)

        if logging:
            elem = model.model_part.Elements[1]
            strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
            # print("strain", strain)
            stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
            # print("stress", stress)
            ifile.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (strain[0][2] - last_ezz, stress[0][0], stress[0][1], stress[0][2], p, dp))
            ifile.flush()

    if logging:
        ifile.close()

    return model

def test():
    main(logging=False, output=False)

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
