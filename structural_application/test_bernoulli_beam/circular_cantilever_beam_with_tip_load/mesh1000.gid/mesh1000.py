##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import mesh1000_include
from mesh1000_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = mesh1000_include.Model('mesh1000',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(ROTATION_Z)
        if abs(node.X0 + 5.0) < tol and abs(node.Y0 - 5.0) < tol:
            tip_node = node
        node.Fix(DISPLACEMENT_Z)
        node.Fix(ROTATION_X)
        node.Fix(ROTATION_Y)

    Fx = -1.0
    Fy = 1.0
    M = 0.0
    tip_node.SetSolutionStepValue(FORCE_X, Fx)
    tip_node.SetSolutionStepValue(FORCE_Y, Fy)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    print("Deflection x at tip node: %.16e" % (tip_node.GetSolutionStepValue(DISPLACEMENT_X)))
    print("Deflection y at tip node: %.16e" % (tip_node.GetSolutionStepValue(DISPLACEMENT_Y)))
    print("Rotation at tip node: %.16e" % (tip_node.GetSolutionStepValue(ROTATION_Z)))

    ## accurate solution, REF: Young et al, Roark's Formula for Stress and Strains, Section 9.2
    R = abs(tip_node.Y0)
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    I = model.model_part.Properties[1].GetValue(LOCAL_INERTIA)[0,0]
    dx = (0.5*Fy*R**3 + (0.75*math.pi-2.0)*Fx*R**3 + (math.pi/2-1)*M*R**2) / (E*I)
    dy = (0.25*math.pi*Fy*R**3 + 0.5*Fx*R**3 + M*R**2) / (E*I)
    theta = (Fy*R**2 + (0.5*math.pi-1.0)*Fx*R**2 + 0.5*math.pi*M*R) / (E*I)
    print("Correct deflection: %.16e, %.16e" % (dx, dy))
    print("Correct rotation: %.16e" % (theta))

    print("Error deflection x: %.16e" % (tip_node.GetSolutionStepValue(DISPLACEMENT_X) - dx))
    print("Error deflection y: %.16e" % (tip_node.GetSolutionStepValue(DISPLACEMENT_Y) - dy))
    print("Error rotation: %.16e" % (tip_node.GetSolutionStepValue(ROTATION_Z) + theta))

    return model

def test():
    model = main(output=False, logging=False)

    ######### pytesting results #########

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 + 5.0) < tol and abs(node.Y0 - 5.0) < tol:
            tip_node = node

    ## accurate solution, REF: Young et al, Roark's Formula for Stress and Strains, Section 9.2
    Fx = -1.0
    Fy = 1.0
    M = 0.0
    R = abs(tip_node.Y0)
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    I = model.model_part.Properties[1].GetValue(LOCAL_INERTIA)[0,0]
    dx_ref = (0.5*Fy*R**3 + (0.75*math.pi-2.0)*Fx*R**3 + (math.pi/2-1)*M*R**2) / (E*I)
    dy_ref = (0.25*math.pi*Fy*R**3 + 0.5*Fx*R**3 + M*R**2) / (E*I)
    theta_ref = (Fy*R**2 + (0.5*math.pi-1.0)*Fx*R**2 + 0.5*math.pi*M*R) / (E*I)

    dx = tip_node.GetSolutionStepValue(DISPLACEMENT_X)
    dy = tip_node.GetSolutionStepValue(DISPLACEMENT_Y)
    theta = tip_node.GetSolutionStepValue(ROTATION_Z)
    assert(abs(dx - dx_ref) < 4e-6)
    assert(abs(dy - dy_ref) < 4e-6)
    assert(abs(theta + theta_ref) < 8e-9)

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
