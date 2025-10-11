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

def write_header(ifile, data):
    for i in range(0, len(data)-1):
        ifile.write("%-*s" % (20, data[i]))
    ifile.write("%s\n" % (data[-1]))

def write_values(ifile, data):
    for i in range(0, len(data)-1):
        ifile.write("%-*.6e" % (20, data[i]))
    ifile.write("%.6e\n" % (data[-1]))

def main(output=True, logging=True, paxial = -450.0, nsteps_axial=40):

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

    ## preconsolidation pressure
    pp = [100.0]*8
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, pp, model.model_part.ProcessInfo)

    time = 0.0
    model.SolveModel(time)
    # sys.exit(0)

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

    if logging:
        elem = model.model_part.Elements[1]
        strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
        # print("strain", strain)
        stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
        pp = elem.CalculateOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, model.model_part.ProcessInfo)
        # print("stress", stress)
        last_ezz = strain[0][2]
        last_sxx = stress[0][0]
        last_syy = stress[0][1]
        last_szz = stress[0][2]
        young_modulus = elem.CalculateOnIntegrationPoints(YOUNG_MODULUS, model.model_part.ProcessInfo)
        last_Eur = young_modulus[0][0]
        last_pp = pp[0][0]

        print("last_ezz: %.10e" % last_ezz)
        print("last_sxx: %.10e" % last_sxx)
        print("last_syy: %.10e" % last_syy)
        print("last_szz: %.10e" % last_szz)
        print("last_Eur: %.10e" % last_Eur)
        print("last_pp: %.10e" % last_pp)
    for node in model.model_part.Nodes:
        print(node.GetSolutionStepValue(DISPLACEMENT))
    # sys.exit(0)
    print("## axial compression ##")

    if logging:
        rad_to_deg = 180.0/math.pi

        elem = model.model_part.Elements[1]

        theta = elem.CalculateOnIntegrationPoints(LODE_ANGLE, model.model_part.ProcessInfo)
        phim = elem.CalculateOnIntegrationPoints(MOBILIZED_FRICTION_ANGLE, model.model_part.ProcessInfo)
        psim = elem.CalculateOnIntegrationPoints(MOBILIZED_DILATANCY_ANGLE, model.model_part.ProcessInfo)
        gammapss = elem.CalculateOnIntegrationPoints(PLASTICITY_INDICATOR, model.model_part.ProcessInfo)
        pp = elem.CalculateOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, model.model_part.ProcessInfo)
        pi = elem.CalculateOnIntegrationPoints(PRESSURE_P, model.model_part.ProcessInfo)
        qi = elem.CalculateOnIntegrationPoints(PRESSURE_Q, model.model_part.ProcessInfo)
        Ei = elem.CalculateOnIntegrationPoints(PRIMARY_LOADING_STIFFNESS, model.model_part.ProcessInfo)
        beta = elem.CalculateOnIntegrationPoints(HS_BETA_FACTOR, model.model_part.ProcessInfo)

        list_of_values = ["epsilon_zz", "sigma_xx", "sigma_yy", "sigma_zz", "p", "dp", "pi", "qi", "Ei", "Eur", "theta", "phim", "psim", "gammapss", "pp", "beta", "plastic_mode"]

        ifile = open("results.txt", "w")
        write_header(ifile, list_of_values)
        values = [0.0, last_sxx, last_syy, last_szz, piso, 0.0, pi[0][0], qi[0][0], Ei[0][0], last_Eur, theta[0][0]*rad_to_deg, phim[0][0]*rad_to_deg, psim[0][0]*rad_to_deg, gammapss[0][0], pp[0][0], beta[0][0], 0]
        write_values(ifile, values)
        ifile.flush()

    nsteps = nsteps_axial
    # nsteps = 20
    # nsteps = 40
    # nsteps = 80
    # nsteps = 200
    dp = (paxial - piso) / nsteps
    time = 2*time
    timeold = time
    delta_time = abs(dp)
    pold = piso
    step = 0
    dp_list = []
    converged_dp_list = []
    adaptive = True
    while True:
        if (pold - paxial) < abs(dp):
            dp = paxial - pold
            delta_time = abs(dp)
        p = pold + dp
        time = timeold + delta_time

        print("###Axial loading step " + str(step+1) + ", p = " + str(p) + "###")

        for cond_id in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)

        converged = model.SolveModel(time)
        if converged:
            converged_dp_list.append(dp)
        if adaptive and not converged:
            dp_list.append(dp)
            dp = dp*0.5
            delta_time = abs(dp)
            print("The analysis does not converged at loading step p = %f, repeat with smaller time step dp = %f" % (p, dp))
            continue

        if output:
            model.WriteOutput(time)

        if logging:
            elem = model.model_part.Elements[1]
            strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
            # print("strain", strain)
            stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
            # print("stress", stress)
            young_modulus = elem.CalculateOnIntegrationPoints(YOUNG_MODULUS, model.model_part.ProcessInfo)
            theta = elem.CalculateOnIntegrationPoints(LODE_ANGLE, model.model_part.ProcessInfo)
            phim = elem.CalculateOnIntegrationPoints(MOBILIZED_FRICTION_ANGLE, model.model_part.ProcessInfo)
            psim = elem.CalculateOnIntegrationPoints(MOBILIZED_DILATANCY_ANGLE, model.model_part.ProcessInfo)
            gammapss = elem.CalculateOnIntegrationPoints(PLASTICITY_INDICATOR, model.model_part.ProcessInfo)
            pp = elem.CalculateOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, model.model_part.ProcessInfo)
            pi = elem.CalculateOnIntegrationPoints(PRESSURE_P, model.model_part.ProcessInfo)
            qi = elem.CalculateOnIntegrationPoints(PRESSURE_Q, model.model_part.ProcessInfo)
            Ei = elem.CalculateOnIntegrationPoints(PRIMARY_LOADING_STIFFNESS, model.model_part.ProcessInfo)
            beta = elem.CalculateOnIntegrationPoints(HS_BETA_FACTOR, model.model_part.ProcessInfo)
            pm = elem.CalculateOnIntegrationPoints(PLASTIC_MODE, model.model_part.ProcessInfo)
            values = [strain[0][2] - last_ezz, stress[0][0], stress[0][1], stress[0][2], p, dp, pi[0][0], qi[0][0], Ei[0][0], young_modulus[0][0], theta[0][0]*rad_to_deg, phim[0][0]*rad_to_deg, psim[0][0]*rad_to_deg, gammapss[0][0], pp[0][0], beta[0][0], pm[0][0]]
            write_values(ifile, values)
            ifile.flush()

        # if converged, then restore the last (converged) increment
        if adaptive:
            if (len(dp_list) > 0):
                dp = dp_list.pop()
            delta_time = abs(dp)

        pold = p
        timeold = time
        step += 1

        # run until the load reach the expected load step
        if p <= paxial:
            print("Stop, p = %e, dp = %e" % (p, dp))
            break

        # if step == 15:
        #     break

    if logging:
        ifile.close()

    for node in model.model_part.Nodes:
        print(node.GetSolutionStepValue(DISPLACEMENT))

    print("list of converged steps: ", converged_dp_list)

    return model

def test():
    model = main(logging=False, output=False, paxial=-450.0, nsteps_axial=40)

    ###### pytesting results
    last_ezz = -4.2384959002e-03
    ref_dezz = -2.4779123219e-02
    elem = model.model_part.Elements[1]
    strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
    dezz = strain[0][2] - last_ezz
    print("%.10e" % dezz)
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
