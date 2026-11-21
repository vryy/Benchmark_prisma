##################################################################
## Testing the B-bar method for near incompressible problem (nu ~= 0.5)
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
sys.path.append(os.getcwd())
import mesh8x4_include
try:
    from mesh8x4_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

sys.path.append(os.getcwd() + "/..")
import analytical_solution

####COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
def compute_L2_error(model_part, a, b, P, E, nu):
    nom = 0.0
    denom = 0.0
    for element in model_part.Elements:
       if element.GetValue(IS_INACTIVE) == False:
           u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model_part.ProcessInfo)
           J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
           Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
           W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
           for i in range(0, len(u)):
               ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1], a, b, P, 0.0, E, nu)
               nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2)) * W[i][0] * J0[i][0]
               denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
    return math.sqrt(nom/denom)

####COMPUTE GLOBAL STRESS (H1) ERROR###
def compute_H1_error(model_part, a, b, P, E, nu):
    nom = 0.0
    denom = 0.0
    for element in model_part.Elements:
       if element.GetValue(IS_INACTIVE) == False:
           o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model_part.ProcessInfo)
           J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
           Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
           W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
           for i in range(0, len(o)):
               ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], a, b, P, 0.0, E, nu)
               nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
               denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
    return math.sqrt(nom/denom)

def main(logging=True, output=True):

    model = mesh8x4_include.Model('mesh8x4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if (abs(node.Z0) < tol) or (abs(node.Z0-10.) < tol):
            node.Fix(DISPLACEMENT_Z)

    P = 1e3
    for node_id in model.layer_nodes_sets["inner"]:
        node = model.model_part.Nodes[node_id]
        L = math.sqrt(node.X0*node.X0 + node.Y0*node.Y0)
        cosa = node.X0 / L
        sina = node.Y0 / L
        node.SetSolutionStepValue(FACE_LOAD_X, P*cosa)
        node.SetSolutionStepValue(FACE_LOAD_Y, P*sina)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    a = 100.
    b = 200.
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
    l2_error = compute_L2_error(model.model_part, a, b, P, E, nu)
    print("Global displacement (L2) error: %.16e" % (l2_error))
    h1_error = compute_H1_error(model.model_part, a, b, P, E, nu)
    print("Global stress (H1) error: %.16e" % (h1_error))

    return model, l2_error, h1_error

def test():
    model, l2_error, h1_error = main(logging=False, output=False)

    l2_error_ref = 2.1163561514303314e-04
    h1_error_ref = 4.3372083260983128e-03

    assert(abs(l2_error - l2_error_ref) < 1e-10)
    assert(abs(h1_error - h1_error_ref) < 1e-10)
    print("Test passed")

def tag():
    tags = "bbar"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
