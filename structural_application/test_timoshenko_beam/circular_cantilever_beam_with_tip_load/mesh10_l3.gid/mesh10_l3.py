##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025 by Hoang-Giang Bui (UoB)                 #####
#####          2026 by Hoang-Giang Bui (DU)                  #####
##### all rights reserved                                    #####
##################################################################
## This file is generated on Wed Mar 11 06:01:37 PM GMT 2026
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import mesh10_l3_include as simulation_include
try:
    from mesh10_l3_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'mesh10_l3'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def get_analytical_solution(Fx, Fy, M, E, I, R):
    ## accurate solution, REF: Young et al, Roark's Formula for Stress and Strains, Section 9.2
    dx = (0.5*Fy*R**3 + (0.75*math.pi-2.0)*Fx*R**3 + (math.pi/2-1)*M*R**2) / (E*I)
    dy = (0.25*math.pi*Fy*R**3 + 0.5*Fx*R**3 + M*R**2) / (E*I)
    theta = (Fy*R**2 + (0.5*math.pi-1.0)*Fx*R**2 + 0.5*math.pi*M*R) / (E*I)
    return dx, dy, theta

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
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

    R = abs(tip_node.Y0)
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    I = model.model_part.Properties[1].GetValue(INERTIA_Y)
    dx, dy, theta = get_analytical_solution(Fx, Fy, M, E, I, R)
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

    Fx = -1.0
    Fy = 1.0
    M = 0.0
    R = abs(tip_node.Y0)
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    I = model.model_part.Properties[1].GetValue(INERTIA_Y)
    dx_ref, dy_ref, theta_ref = get_analytical_solution(Fx, Fy, M, E, I, R)

    dx = tip_node.GetSolutionStepValue(DISPLACEMENT_X)
    dy = tip_node.GetSolutionStepValue(DISPLACEMENT_Y)
    theta = tip_node.GetSolutionStepValue(ROTATION_Z)
    assert(abs(dx - dx_ref) < 2e-4)
    assert(abs(dy - dy_ref) < 2e-4)
    assert(abs(theta + theta_ref) < 7e-5)

    #####################################
    print("Test passed")

def tag():
    tags = "beam,timoshenko"
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
