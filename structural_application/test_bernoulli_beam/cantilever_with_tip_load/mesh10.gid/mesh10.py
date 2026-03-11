##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import mesh10_include
from mesh10_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = mesh10_include.Model('mesh10',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(ROTATION_Z)
        if abs(node.X0 - 1.0) < tol:
            tip_node = node
        if abs(node.X0 - 0.0) < tol:
            origin_node = node
        node.Fix(DISPLACEMENT_Z)
        node.Fix(ROTATION_X)
        node.Fix(ROTATION_Y)

    F = -1.0
    tip_node.SetSolutionStepValue(FORCE_Y, F)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    variable_transfer_util = VariableTransferUtility()
    variable_transfer_util.TransferVariablesToNodes(model.model_part, MOMENT)
    variable_transfer_util.TransferVariablesToNodes(model.model_part, FORCE)

    print("-----------------")
    print("Deflection at tip node: %.16e" % (tip_node.GetSolutionStepValue(DISPLACEMENT_Y)))
    print("Rotation at tip node: %.16e" % (tip_node.GetSolutionStepValue(ROTATION_Z)))
    print("Shear force at origin node: %.16e" % (origin_node.GetSolutionStepValue(FORCE_Y)))
    print("Shear force at tip node: %.16e" % (tip_node.GetSolutionStepValue(FORCE_Y)))
    print("Moment at origin node: %.16e" % (origin_node.GetSolutionStepValue(MOMENT_Z)))
    print("Moment at tip node: %.16e" % (tip_node.GetSolutionStepValue(MOMENT_Z)))
    print("-----------------")

    ## accurate solution, REF: https://www.engineeringtoolbox.com/cantilever-beams-d_1848.html
    L = tip_node.X0
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    I = model.model_part.Properties[1].GetValue(LOCAL_INERTIA)[0,0]
    d = F*L**3/(3*E*I)
    r = F*L**2/(2*E*I)
    print("Correct deflection: %.16e" % (d))
    print("Deflection error: %.16e" % (d - tip_node.GetSolutionStepValue(DISPLACEMENT_Y)))
    print("Correct rotation: %.16e" % (r))
    print("Rotation error: %.16e" % (r - tip_node.GetSolutionStepValue(ROTATION_Z)))
    print("-----------------")
    for node in model.model_part.Nodes:
        print("moment at node %d, X = %.6e: %.16e" % (node.Id, node.X0, node.GetSolutionStepValue(MOMENT_Z)))
        print(" shear force at node %d: %.16e" % (node.Id, node.GetSolutionStepValue(FORCE_Y)))

    # for elem in model.model_part.Elements:
    #     print(elem.GetValuesOnIntegrationPoints(MOMENT, model.model_part.ProcessInfo))
    print("-----------------")

    return model

def test():
    model = main(output=False, logging=False)

    ######### pytesting results #########

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            tip_node = node

    F = -1.0
    L = tip_node.X0
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    I = model.model_part.Properties[1].GetValue(LOCAL_INERTIA)[0,0]
    d_ref = F*L**3/(3*E*I)
    r_ref = F*L**2/(2*E*I)

    d = tip_node.GetSolutionStepValue(DISPLACEMENT_Y)
    r = tip_node.GetSolutionStepValue(ROTATION_Z)

    assert(abs(d - d_ref) < 1e-16)
    assert(abs(r - r_ref) < 1e-16)

    for node in model.model_part.Nodes:
        x = node.X0
        m = node.GetSolutionStepValue(MOMENT_Z)
        q = node.GetSolutionStepValue(FORCE_Y)
        assert(abs(q + 1.0) < 1e-11)        # shear force = constant = -1
        assert(abs(m - (1.0-x)) < 1e-12)    # moment is linear = 1-x
    #####################################
    print("Test passed")

def tag():
    tags = "beam,bernoulli"
    return tags

def print_tag():
    print("Tag(s): " + tag())

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
