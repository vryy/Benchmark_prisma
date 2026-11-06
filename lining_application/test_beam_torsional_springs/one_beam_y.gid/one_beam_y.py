##################################################################
import sys
import os
import math
##################################################################
##################################################################
sys.path.append('./one_beam_y.gid')
import one_beam_y_include
from one_beam_y_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = one_beam_y_include.Model('one_beam_y',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 10.0) < tol:
            load_node = node
        if abs(node.Y0 - 0.0) < tol:
            constraint_node = node
    load_node.SetSolutionStepValue(FORCE_X, 10.0)

    constraint_node.Fix(DISPLACEMENT_X)
    constraint_node.Fix(DISPLACEMENT_Y)
    constraint_node.Fix(DISPLACEMENT_Z)
    constraint_node.Fix(ROTATION_X)
    # constraint_node.Fix(ROTATION_Y)
    constraint_node.Fix(ROTATION_Z)

    last_cond_id = len(model.model_part.Conditions) + 1
    prop = model.model_part.Properties[2]
    prop[TORTIONAL_STIFFNESS] = 1.0e5
    for elem in model.model_part.Elements:
        nodes = elem.GetNodes()
        new_cond = BeamTorsionalSprings(last_cond_id, nodes[0], nodes[1], prop)
        last_cond_id = last_cond_id + 1
        model.model_part.AddCondition(new_cond)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 10.0) < tol:
            load_node = node

    ######### pytesting results #########
    u_tip = load_node.GetSolutionStepValue(DISPLACEMENT_X)
    theta_tip = load_node.GetSolutionStepValue(ROTATION_Z)
    u_ref = 1.5873015873015872e-02
    theta_ref = -2.3809523809523799e-03
    print("disp on tip: %.16e" % u_tip)
    print("rotation on tip: %.16e" % theta_tip)
    assert(abs(u_tip - u_ref) < 1e-10)
    assert(abs(theta_tip - theta_ref) < 1e-10)
    #####################################
    print("Test passed")

def tag():
    tags = ""
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
