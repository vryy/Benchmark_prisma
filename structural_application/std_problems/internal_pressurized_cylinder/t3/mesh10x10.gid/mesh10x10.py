##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import mesh10x10_include
from mesh10x10_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

sys.path.append("../..")
import analytical_solution

def run(model, output=True):

    # boundary condition

    P = 1.0e2
    load = Vector(3)
    for node in model.model_part.Nodes:
        x = node.X0
        y = node.Y0
        r = math.sqrt(x*x + y*y)
        t = math.acos(x/r)
        load[0] = P * math.cos(t)
        load[1] = P * math.sin(t)
        load[2] = 0.0
        node.SetSolutionStepValue(FACE_LOAD, load)

    tol = 1e-06
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)

    # analysis

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print("node " + str(node.Id) + " DISPLACEMENT: " + str(node.GetSolutionStepValue(DISPLACEMENT)))

    tol = 1.0e-6
    a = 100.0
    b = 200.0
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)

    ##COMPARE THE SOLUTION AT 1 POINT
    for node in model.model_part.Nodes:
        if abs(node.X0 - b) < tol and abs(node.Y0 - 0.0) < tol:
            ux = node.GetSolutionStepValue(DISPLACEMENT_X)
            uy = node.GetSolutionStepValue(DISPLACEMENT_Y)
            u_ref = analytical_solution.get_displacement(node.X0, node.Y0, a, b, P, 0.0, E, nu)
            error_disp_b_0 = math.sqrt(math.pow(ux-u_ref[0], 2) + math.pow(uy-u_ref[1], 2))
            disp_b_0_norm_ana = math.sqrt(math.pow(u_ref[0], 2) + math.pow(u_ref[1], 2))
            print("rel_error_disp_b_0: %.16e" % (error_disp_b_0/disp_b_0_norm_ana))

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
    nom = 0.0
    denom = 0.0
    for element in model.model_part.Elements:
        u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        for i in range(0, len(u)):
            ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1], a, b, P, 0.0, E, nu)
            nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2)) * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
    # print("nom:", nom, "denom:", denom)
    model.l2_error = math.sqrt(nom / denom)
    print("Global displacement (L2) error: %.16e" % model.l2_error)

    ###COMPUTE GLOBAL STRESS (H1) ERROR###
    nom = 0.0
    denom = 0.0
    for element in model.model_part.Elements:
        o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        for i in range(0, len(o)):
            ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], a, b, P, 0.0, E, nu)
            nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
    model.h1_error = math.sqrt(nom / denom)
    print("Global stress (H1) error: %.16e" % model.h1_error)

    return model

def main(logging=True, output=True):

    model = mesh10x10_include.Model('mesh10x10',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    model = run(model, output=output)

    return model

def test():

    model = main(logging=False, output=False)
    l2_error_ref = 1.7467908086117562e-02
    h1_error_ref = 9.1725497539060563e-02
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    print("Test passed")

def tag():
    return "kinematic_linear"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=False, logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
