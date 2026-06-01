##################################################################
import sys
import os
import math
##################################################################
import curved_Lshape_include
from curved_Lshape_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def compute_L2_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
            for i in range(0, len(u)):
                ana_u = solution.GetTemperatureAt(Q[i][0], Q[i][1], Q[i][2])
                nom = nom + pow(u[i][0] - ana_u, 2) * W[i][0] * J0[i][0]
                denom = denom + pow(ana_u, 2) * W[i][0] * J0[i][0]
    print("nom:", nom)
    print("denom:", denom)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

def main(logging=True, output=True):
    model = curved_Lshape_include.Model('curved_Lshape',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    util = HeatStdProblem1Solution()

    ###### DEBUGGING
    #for node in model.model_part.Nodes:
    #    if abs(node.Y0) < 0.04:
    #        ana_temp = util.CalculateTemperatureOnNode(node)
    #sys.exit(0)
    ###### END DEBUGGING

    for cond in model.model_part.Conditions:
        cond.SetValue(ACTIVATION_LEVEL, -1)

    for node_id in model.layer_nodes_sets['boundary']:
        node = model.model_part.Nodes[node_id]
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, 0.0)

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    # if output:
    #     model.WriteOutput(time)

    print(model.layer_nodes_sets['boundary'])
    for node_id in model.layer_nodes_sets['boundary']:
        node = model.model_part.Nodes[node_id]
        ana_temp = util.CalculateTemperatureOnNode(node)
        node.SetSolutionStepValue(TEMPERATURE, ana_temp)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)

    for node in model.model_part.Nodes:
        ana_temp = util.CalculateTemperatureOnNode(node)
        temp = node.GetSolutionStepValue(TEMPERATURE)
    #    if abs(ana_temp) < 1.0e-10:
    #        error = temp-ana_temp
    #    else:
    #        error = (temp-ana_temp)/abs(ana_temp)
        error = temp - ana_temp
        node.SetSolutionStepValue(TEMPERATURE_ERROR, error)
        node.SetSolutionStepValue(REFERENCE_TEMPERATURE, ana_temp)

    error = compute_L2_error(model.model_part.Elements, util, model.model_part.ProcessInfo)
    model.l2_error = error
    print("Global L2 error: %.16e" % error)

    # eigen_solver = SpectraEigenvaluesSolver()
    # values = eigen_solver.SolveLargestSym(model.solver.solver.A, 3)
    # print("largest eigenvalues:", values)
    # values = eigen_solver.SolveSmallestSPD(model.solver.solver.A, SuperLUSolver(), 2)
    # print("smallest eigenvalues:", values)

    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)
    l2_error_ref = 9.7263405426734311e-04
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    print("Test passed")

def tag():
    return "thermal"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True)
