# Test for reduced integration
# Ref:
# + Nakshatrala, On finite element formulations for nearly incompressible linear elasticity
# + Chennakesava, Mixed Galerkin and least-squares formulations for isogeometric analysis.
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import mesh_8x8_include
from mesh_8x8_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True, integration_order=1):

    model = mesh_8x8_include.Model('mesh_8x8',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel(integration_order=integration_order)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)

        if abs(node.X0-48.0) < tol:
            node.SetSolutionStepValue(FACE_LOAD_Y, 100.0/16.0)

        if abs(node.X0-48.0) < tol and abs(node.Y0-60.0) < tol:
            top_node = node

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    print("Vertical displacement of top-right node: %.16e" % (top_node.GetSolutionStepValue(DISPLACEMENT_Y)))

    return model

def test1():
    model = main(logging=False, output=False, integration_order=1)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0-48.0) < tol and abs(node.Y0-60.0) < tol:
            top_node = node

    u_top = top_node.GetSolutionStepValue(DISPLACEMENT_Y)
    # the reduced integration significantly improve the tip displacement
    u_top_ref = 8.0062373809037357e+00

    assert(abs(u_top - u_top_ref) < 1e-10)
    print("Test passed")

def test2():
    model = main(logging=False, output=False, integration_order=2)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0-48.0) < tol and abs(node.Y0-60.0) < tol:
            top_node = node

    u_top = top_node.GetSolutionStepValue(DISPLACEMENT_Y)
    # due to voluemtric locking because of nearly incompressible condition (nu ~= 0.5),
    # the displacement at the tip is significantly smaller
    u_top_ref = 2.2236640506946190e+00

    assert(abs(u_top - u_top_ref) < 1e-10)
    print("Test passed")

def test():
    test1()
    test2()

def tag():
    tags = "EAS"
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
