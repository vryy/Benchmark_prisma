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
    # sys.exit(0)

    ## DEBUGGING
    pp = [1.0e10]*8
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, pp, model.model_part.ProcessInfo)
    ## END DEBUGGING

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
            cond.SetValue(POSITIVE_FACE_PRESSURE, 0.0)

        for cond_id in model.layer_cond_sets['loady']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)
            cond.SetValue(POSITIVE_FACE_PRESSURE, 0.0)

        for cond_id in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)
            cond.SetValue(POSITIVE_FACE_PRESSURE, 0.0)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

    if logging:
        elem = model.model_part.Elements[1]
        strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
        # print("strain", strain)
        stress = elem.CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
        # print("stress", stress)
        last_ezz = strain[0][2]
        last_sxx = stress[0][0]
        last_syy = stress[0][1]
        last_szz = stress[0][2]
        young_modulus = elem.CalculateOnIntegrationPoints(YOUNG_MODULUS, model.model_part.ProcessInfo)
        last_Eur = young_modulus[0][0]

    print("last_ezz:", last_ezz)
    print("last_sxx:", last_sxx)
    print("last_syy:", last_syy)
    print("last_szz:", last_szz)
    print("last_Eur:", last_Eur)
    for node in model.model_part.Nodes:
        print(node.GetSolutionStepValue(DISPLACEMENT))
    # sys.exit(0)
    print("## axial compression ##")

    if logging:
        rad_to_deg = 1 #180.0/math.pi

        elem = model.model_part.Elements[1]

        strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
        eep = elem.CalculateOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_ELASTIC_STRAIN, model.model_part.ProcessInfo)
        eeq = elem.CalculateOnIntegrationPoints(EQUIVALENT_DEVIATORIC_ELASTIC_STRAIN, model.model_part.ProcessInfo)
        epp = elem.CalculateOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN, model.model_part.ProcessInfo)
        epq = elem.CalculateOnIntegrationPoints(EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN, model.model_part.ProcessInfo)
        ep = elem.CalculateOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_STRAIN, model.model_part.ProcessInfo)
        eq = elem.CalculateOnIntegrationPoints(EQUIVALENT_DEVIATORIC_STRAIN, model.model_part.ProcessInfo)

        theta = elem.CalculateOnIntegrationPoints(LODE_ANGLE, model.model_part.ProcessInfo)
        phim = elem.CalculateOnIntegrationPoints(MOBILIZED_FRICTION_ANGLE, model.model_part.ProcessInfo)
        psim = elem.CalculateOnIntegrationPoints(MOBILIZED_DILATANCY_ANGLE, model.model_part.ProcessInfo)
        gammapss = elem.CalculateOnIntegrationPoints(PLASTICITY_INDICATOR, model.model_part.ProcessInfo)
        pp = elem.CalculateOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, model.model_part.ProcessInfo)
        pi = elem.CalculateOnIntegrationPoints(PRESSURE_P, model.model_part.ProcessInfo)
        qi = elem.CalculateOnIntegrationPoints(PRESSURE_Q, model.model_part.ProcessInfo)
        Ei = elem.CalculateOnIntegrationPoints(PRIMARY_LOADING_STIFFNESS, model.model_part.ProcessInfo)

        ifile = open("results.txt", "w")
        ifile.write("epsilon_zz\tsigma_xx\tsigma_yy\tsigma_zz\tpres\t\t\tdp\t\t\tpi\t\t\tqi\t\t\tEi\t\t\t\tEur\t\t\t\ttheta\t\tphim\t\tpsim\t\tgammapss\tpp\n")
        ifile.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (0.0, last_sxx, last_syy, last_szz, piso, 0.0, pi[0][0], qi[0][0], Ei[0][0], last_Eur, theta[0][0], phim[0][0], psim[0][0], gammapss[0][0], pp[0][0]))
        ifile.flush()

        monfile = open("sigma_xx.txt", "w")
        monfile.write("p0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\n")
        monfile.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n" % (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

        monfile1 = open("plastic_mode.txt", "w")
        monfile1.write("p0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\n")
        monfile1.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (0, 0, 0, 0, 0, 0, 0, 0))

        monfile2 = open("epsilon.txt", "w")
        monfile2.write("ee_p\t\t\tee_q\t\t\tep_p\t\t\tep_q\t\t\te_p\t\t\te_q\t\t\te_xx\te_yy\te_zz\te_xy\te_yz\te_xz\n")
        monfile2.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n" % (eep[0][0], eeq[0][0], epp[0][0], epq[0][0], ep[0][0], eq[0][0], strain[0][0], strain[0][1], strain[0][2], strain[0][3], strain[0][4], strain[0][5]))
        monfile2.flush()

    nsteps = 50 #45 #20 #10 #4 #5 #20
    dp = (paxial - piso) / nsteps
    time = 2*time
    timeold = time
    delta_time = abs(dp)
    pold = piso
    step = 0
    dp_list = []
    while p > paxial:
        p = pold + dp
        time = timeold + delta_time

        print("###Axial loading step " + str(step+1) + ", p = " + str(p) + ", time = " + str(time) + "###")

        for cond_id in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[cond_id]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, p)

        converged = model.SolveModel(time)
        # if not converged:
        #     dp_list.append(dp)
        #     dp = dp*0.5
        #     delta_time = abs(dp)
        #     print("The analysis does not converge at loading step p = %f, repeat with smaller time step dp = %f" % (p, dp))
        #     continue

        if output:
            model.WriteOutput(time)

        if logging:
            elem = model.model_part.Elements[1]
            strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
            eep = elem.CalculateOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_ELASTIC_STRAIN, model.model_part.ProcessInfo)
            eeq = elem.CalculateOnIntegrationPoints(EQUIVALENT_DEVIATORIC_ELASTIC_STRAIN, model.model_part.ProcessInfo)
            epp = elem.CalculateOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN, model.model_part.ProcessInfo)
            epq = elem.CalculateOnIntegrationPoints(EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN, model.model_part.ProcessInfo)
            ep = elem.CalculateOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_STRAIN, model.model_part.ProcessInfo)
            eq = elem.CalculateOnIntegrationPoints(EQUIVALENT_DEVIATORIC_STRAIN, model.model_part.ProcessInfo)
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
            cs = elem.CalculateOnIntegrationPoints(CONVERGENCE_STATE, model.model_part.ProcessInfo)
            ifile.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (strain[0][2] - last_ezz, stress[0][0], stress[0][1], stress[0][2], p, dp, pi[0][0], qi[0][0], Ei[0][0], young_modulus[0][0], theta[0][0], phim[0][0], psim[0][0], gammapss[0][0], pp[0][0]))
            ifile.flush()

            monfile.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n" % (stress[0][0], stress[1][0], stress[2][0], stress[3][0], stress[4][0], stress[5][0], stress[6][0], stress[7][0]))
            monfile.flush()

            pm = elem.CalculateOnIntegrationPoints(PLASTIC_MODE, model.model_part.ProcessInfo)
            monfile1.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (pm[0][0], pm[1][0], pm[2][0], pm[3][0], pm[4][0], pm[5][0], pm[6][0], pm[7][0]))
            monfile1.flush()

            monfile2.write("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n" % (eep[0][0], eeq[0][0], epp[0][0], epq[0][0], ep[0][0], eq[0][0], strain[0][0], strain[0][1], strain[0][2], strain[0][3], strain[0][4], strain[0][5]))
            monfile2.flush()

            print("Convergence state of constitutive law: ", cs)

        # # if converged, then restore the last (converged) increment
        # if (len(dp_list) > 0):
        #     dp = dp_list.pop()
        # delta_time = abs(dp)

        pold = p
        timeold = time
        step += 1

        # if step == 130:
        # if step == 47:
        # if step == 20:
        # if step == 15:
        # if step == 17:
        # if step == 43:
        # if step == 2:
        #     break

    if logging:
        ifile.close()
        monfile.close()
        monfile1.close()
        monfile2.close()

    for node in model.model_part.Nodes:
        print(node.GetSolutionStepValue(DISPLACEMENT))

    return model

# def test():
#     model = main(logging=False, output=False, paxial = -300.0)

#     ###### pytesting results
#     last_ezz = -0.00530946088843
#     ref_dezz = -0.00672446661269 #-0.00677878711178
#     elem = model.model_part.Elements[1]
#     strain = elem.CalculateOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
#     dezz = strain[0][2] - last_ezz
#     print(dezz)
#     assert(abs(dezz - ref_dezz) / abs(ref_dezz) < 1e-10)
#     ########################

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
