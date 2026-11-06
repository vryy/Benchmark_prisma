##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import full_beam_include
from full_beam_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

timer = Timer()
start_time = time_module.time()

tunnel_radius = 4.7
cx = 0.0
cy = -20.0
thick = 0.4
width = 2.0 # to match with FEMSL

g = 9.81
k0_value = 0.5
sigma_v = 0.34
sigma_h = k0_value*sigma_v

ground_modulus = 150.0
ground_density = 17.0e-3 / g
poisson_ratio = 0.3
beta_factor = 2.0
tangent_spring_factor = 1.0/3
cohesion = 0.005
friction_angle = 37/57.297

def main(output=True, logging=True, number_of_elements = 100*4):

    print(tunnel_radius)

    lining_util = LiningUtility()
    area = thick*width
    inertia = Matrix(3, 3)
    inertia[0, 0] = width*(thick**3)/12
    inertia[0, 1] = 0.0
    inertia[0, 2] = 0.0
    inertia[1, 0] = 0.0
    inertia[1, 1] = width*(thick**3)/12
    inertia[1, 2] = 0.0
    inertia[2, 0] = 0.0
    inertia[2, 1] = 0.0
    inertia[2, 2] = 0.0
    body_force = Vector(3)
    body_force[0] = 0.0
    body_force[1] = 0.0
    body_force[2] = 0.0

    ## Generate the mesh

    model_part = ModelPart("Beam_HRM_" + str(number_of_elements))

    ## Beam properties

    model_part.Properties[1].SetValue(YOUNG_MODULUS, 35.0e3)
    model_part.Properties[1].SetValue(POISSON_RATIO, 0.15)
    model_part.Properties[1].SetValue(THICKNESS, width)
    model_part.Properties[1].SetValue(BODY_FORCE, body_force)

    import structural_solver_advanced
    structural_solver_advanced.AddVariables(model_part)
    model_part.AddNodalSolutionStepVariable(HYPERSTATIC_FORCE)
    beam_elements = lining_util.GenerateCircularMesh(model_part, "BeamElement3D2N", cx, cy, tunnel_radius, number_of_elements, model_part.Properties[1])

    for elem in beam_elements:
        elem.SetValue(AREA, area)
        elem.SetValue(INERTIA, inertia)

    model = full_beam_include.Model(model_part,os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## Boundary condition
    tol = 1.0e-6
    for node in model.model_part.Nodes:
    #    if abs(node.Y0 - 0.0) < tol:
    #        node.Fix(DISPLACEMENT_Y)
    #        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        node.Fix(DISPLACEMENT_Z)
        node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        node.Fix(ROTATION_X)
        node.SetSolutionStepValue(ROTATION_X, 0.0)
        node.Fix(ROTATION_Y)
        node.SetSolutionStepValue(ROTATION_Y, 0.0)

    ## setup geotechnical load

    load_prop = model.model_part.Properties[1]

    sample_cond_name = "GeotechnicalLoadI2D2N"
    last_cond_id = len(model.model_part.Conditions) + 1

    new_load_conds = []
    for elem in model.model_part.Elements:
        new_cond = lining_util.CreateCondition(elem, sample_cond_name, last_cond_id, load_prop)
        new_cond.Initialize(model.model_part.ProcessInfo)
        lining_util.SetLocalFrame(new_cond, cx, cy, model.model_part.ProcessInfo)
        new_cond.SetValue(VERTICAL_LOAD, sigma_v)
        new_cond.SetValue(HORIZONTAL_LOAD, sigma_h)
        last_cond_id = last_cond_id + 1
        new_load_conds.append(new_cond)

    for cond in new_load_conds:
        model.model_part.AddCondition(cond)

    load_prop.SetValue(INTEGRATION_ORDER, 1)
    load_prop.SetValue(G_CONSTANT, g)
    load_prop.SetValue(K0, k0_value)
    load_prop.SetValue(SOIL_DENSITY, 0.017 / g)

    ## setup hyperstatic bedding
    bedding_prop = model.model_part.Properties[2]
    bedding_prop.SetValue(THICKNESS, width)
    #bedding_prop.SetValue(THICKNESS, width*math.cos(0.25*math.pi))
    #bedding_prop.SetValue(THICKNESS, width/math.cos(0.25*math.pi))

    sample_cond_name = "HyperstaticBeddingSecondOrder2D2N"
    last_cond_id = len(model.model_part.Conditions) + 1

    cosf = math.cos(friction_angle)
    sinf = math.sin(friction_angle)
    pclim = 2*cohesion*cosf/(1.0-sinf)+(1.0+sinf)/(1.0-sinf)*sigma_v*(1.0+k0_value)/2.0*(poisson_ratio/(1.0-poisson_ratio))
    print("pclim", pclim)
    f = math.sqrt(sigma_v**2 + sigma_h**2) # good one #1
    #f = sigma_v*cosf + sigma_h*sinf # good one #2
    hrm_calib_factor = 2*f/(pclim*math.log(pclim**2 / (pclim - f)) + f)
    print("hrm_calib_factor", hrm_calib_factor)
    bedding_prop.SetValue(HRM_CALIBRATION_FACTOR, hrm_calib_factor)

    new_bedding_conds = []
    for cond in model.model_part.Conditions:
        new_cond = lining_util.CreateCondition(cond, sample_cond_name, last_cond_id, bedding_prop)
        new_cond.Initialize(model.model_part.ProcessInfo)
        lining_util.SetLocalFrame(new_cond, cx, cy, model.model_part.ProcessInfo)
        lining_util.SetHyperstaticParams(new_cond, ground_modulus, poisson_ratio, tunnel_radius, beta_factor, tangent_spring_factor, cohesion, friction_angle, sigma_v, sigma_h, model.model_part.ProcessInfo)
        last_cond_id = last_cond_id + 1
        new_bedding_conds.append(new_cond)

    for cond in new_bedding_conds:
        model.model_part.AddCondition(cond)
    print(str(len(new_bedding_conds)) + " conditions of type " + sample_cond_name + " are added successfully")

    ## static solution
    time = 0.0
    delta_time = 1.0
    time = time + delta_time
    model.Solve(time, 0, 0, 0, 0)

    #for node in model.model_part.Nodes:
    #    print("node " + str(node.Id) + " displ: " + str(node.GetSolutionStepValue(DISPLACEMENT_X)) + ", " + str(node.GetSolutionStepValue(DISPLACEMENT_Y)))

    if output:
        # transfer the moment and force to node
        variable_transfer_util = VariableTransferUtility()
        variable_transfer_util.TransferVariablesToNodes(model.model_part, MOMENT)
        variable_transfer_util.TransferVariablesToNodes(model.model_part, FORCE)
        variable_transfer_util.TransferVariablesToNodesForConditions(model.model_part, new_bedding_conds, HYPERSTATIC_FORCE)

        model.WriteOutput(time)

    if logging:
        ## export to data file for plotting
        import hrm_post_utility
        mag_factor = 500.0
        hrm_post_utility.PostDeformedLining("displacement_profile.txt", model.model_part, cx, cy, mag_factor)
        hrm_post_utility.PostRotation("rotation_profile.txt", model.model_part, cx, cy, -math.pi/2)
        hrm_post_utility.PostMoment("moment_profile.txt", model.model_part, cx, cy, -math.pi/2)
        hrm_post_utility.PostForces("forces_profile.txt", model.model_part, cx, cy, -math.pi/2)

    return model, hrm_calib_factor

def test():
    model, hrm_calib_factor = main(logging=False, output=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.Y0 - (cy-tunnel_radius)) < tol:
            bottom_node = node
        if abs(node.Y0 - (cy+tunnel_radius)) < tol:
            top_node = node

    u_top = top_node.GetSolutionStepValue(DISPLACEMENT_Y)
    u_bottom = bottom_node.GetSolutionStepValue(DISPLACEMENT_Y)

    print("displacement at top: %.16e" % u_top)
    print("displacement at bottom: %.16e" % u_bottom)
    print("hrm_calib_factor: %.16e" % hrm_calib_factor)

    ######### pytesting results #########
    u_top_ref = -6.8286555408725129e-03
    u_bottom_ref = 6.8286555408756744e-03
    hrm_calib_factor_ref = 9.1683555240891812e-01
    assert(abs(u_top - u_top_ref) < 1e-10)
    assert(abs(u_bottom - u_bottom_ref) < 1e-10)
    assert(abs(hrm_calib_factor - hrm_calib_factor_ref) < 1e-10)
    #####################################
    print("Test passed")

def tag():
    tags = "vhrm,beam"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)

stop_time = time_module.time()
print(timer)
print("Analysis time: " + str(stop_time - start_time))