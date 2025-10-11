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
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./square10x10.gid')
import square10x10_include
from square10x10_include import *
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

def main(output=True, logging=True):

    model = square10x10_include.Model('square10x10',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 0.0) < tol) or (abs(node.X0 - 1.0) < tol):
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0)
        if (abs(node.Y0 - 0.0) < tol) or (abs(node.Y0 - 1.0) < tol):
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    solution = HeatStdProblem2Solution(1.6, 2.4)
    for node in model.model_part.Nodes:
        ana_temp = solution.CalculateTemperatureOnNode(node)
        temp = node.GetSolutionStepValue(TEMPERATURE)
    #    if abs(ana_temp) < 1.0e-10:
    #        error = temp-ana_temp
    #    else:
    #        error = (temp-ana_temp)/abs(ana_temp)
        error = temp - ana_temp
        node.SetSolutionStepValue(TEMPERATURE_ERROR, error)
        node.SetSolutionStepValue(REFERENCE_TEMPERATURE, ana_temp)

    ###### pytesting results
    error = compute_L2_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
    print("Global L2 error:", error)
    assert(abs(error - 0.021869281302720332) < 1e-10)
    ########################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
