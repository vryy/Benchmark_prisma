##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
import math
##################################################################
##################################################################
import mesh_05_include
from mesh_05_include import *

sys.path.append("../../../")
import analytical_solution

##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = mesh_05_include.Model('mesh_05',os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    P = 10.0
    ri = 1.0

    load_fn_x = LoadFunctionR3RnPlateWithTheHoleX(P, ri)
    load_fn_x.Assign(LOAD_FUNCTION, load_fn_x, model.model_part.Properties[2])
    model.model_part.Properties[2].SetValue(INTEGRATION_ORDER, 1)

    load_fn_y = LoadFunctionR3RnPlateWithTheHoleY(P, ri)
    load_fn_y.Assign(LOAD_FUNCTION, load_fn_y, model.model_part.Properties[3])
    model.model_part.Properties[3].SetValue(INTEGRATION_ORDER, 1)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    ###COMPUTE STRAIN ENERGY###
    strain_energy = 0.0
    for element in model.model_part.Elements:
        strain_energy_at_ip = element.GetValuesOnIntegrationPoints(STRAIN_ENERGY, model.model_part.ProcessInfo)
    #    print(strain_energy)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
    #    print(J0)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
    #    print(W)
        for i in range(0, len(strain_energy_at_ip)):
            strain_energy = strain_energy + W[i][0] * J0[i][0] * strain_energy_at_ip[i][0]
    print("strain_energy:", strain_energy)
    ###########################

    tol = 1.0e-6
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
    G = E/(2*(1+nu))
    kappa = 3-4*nu # plane strain

    ###COMPUTE THE STRAIN ENERGY USING ANALYTICAL SOLUTION ###
    a = 4.0
    f_se = analytical_solution.get_strain_energy_function(P, ri, G, kappa)

    #print(f_se.GetValue(0.0, 4.0))
    #print(f_se.GetValue(2.0, 2.0))
    #print(f_se.GetValue(0.3, 1.2))

    strain_energy_exact = 0.0
    for element in model.model_part.Elements:
        strain_energy_exact = strain_energy_exact + f_se.Integrate(element, 2)
    print("strain_energy_exact:", strain_energy_exact)

    ana_sol = analytical_solution.PlaneStrainSolution(P, ri, G, nu)

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
    nom = 0.0
    denom = 0.0
    for element in model.model_part.Elements:
        u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        for i in range(0, len(u)):
            ana_u = ana_sol.get_displacement(Q[i][0], Q[i][1], 0.0)
            nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2)) * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
    l2_error = math.sqrt(nom / denom)
    print("Global displacement (L2) error: %.16e" % l2_error)

    ###COMPUTE GLOBAL STRESS (H1) ERROR###
    nom = 0.0
    denom = 0.0
    for element in model.model_part.Elements:
        o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        for i in range(0, len(o)):
            ana_o = ana_sol.get_stress_3d(Q[i][0], Q[i][1], 0.0)
            nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
    h1_error = math.sqrt(nom / denom)
    print("Global stress (H1) error: %.16e" % h1_error)

    ####COMPUTE X-DISPLACEMENT AT (A,0) ###
    #a = 4.0
    #for node in model.model_part.Nodes:
    #    if abs(node.X0 - a) < tol and abs(node.Y0 - 0.0) < tol:
    #        disp_x_a_0 = node.GetSolutionStepValue(DISPLACEMENT_X)

    #        ana_disp = analytical_solution.get_displacement(node.X0, node.Y0, P, ri, G, kappa)
    #        print("x-displacement (a,0):", disp_x_a_0)
    #        print("analytical x-displacement (a,0):", ana_disp[0])
    #        print("error x-displacement (a,0):", abs((disp_x_a_0 - ana_disp[0]) / ana_disp[0]))

    ####COMPUTE Y-DISPLACEMENT AT (0,A) ###
    #for node in model.model_part.Nodes:
    #    if abs(node.Y0 - a) < tol and abs(node.X0 - 0.0) < tol:
    #        disp_y_0_a = node.GetSolutionStepValue(DISPLACEMENT_Y)

    #        ana_disp = analytical_solution.get_displacement(node.X0, node.Y0, P, ri, G, kappa)
    #        print("y-displacement (0,a):", disp_y_0_a)
    #        print("analytical y-displacement (0,a):", ana_disp[1])
    #        print("error y-displacement (0,a):", abs((disp_y_0_a - ana_disp[1]) / ana_disp[1]))

    return model, l2_error, h1_error

def test():
    model, l2_error, h1_error = main(output=False, logging=False)

    l2_error_ref = 1.7049624806358438e-02
    h1_error_ref = 6.7102735170833747e-02

    assert(abs(l2_error - l2_error_ref) < 1e-10)
    assert(abs(h1_error - h1_error_ref) < 1e-10)

    #####################
    print("Test passed")

def tag():
    return "std"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True)
