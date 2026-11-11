##################################################################
##################################################################
import sys
import os
import math
import copy
##################################################################
##################################################################
sys.path.append('./mesh_4_320.gid')
import mesh_4_320_include
from mesh_4_320_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

g = 9.81
k0_value = 0.5
cy = 0.0
cz = -20.0
sigma_v = 0.34
sigma_h = k0_value*sigma_v

ground_modulus = 150.0
ground_density = 17.0e-3 / g
poisson_ratio = 0.3
tunnel_radius = 4.7
thickness = 0.4
beta_factor = 2.0
tangent_spring_factor = 1.0/3
cohesion = 0.005
friction_angle = 37/57.297

def main(output=True, logging=True):
    model = mesh_4_320_include.Model('mesh_4_320',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## Boundary condition
    tol = 1.0e-6
    # for node in model.model_part.Nodes:
    #     node.Fix(DISPLACEMENT_X)
    #     node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
    # for node in model.model_part.Nodes:
    #     if abs(node.X0) < tol:
    #         node.Fix(DISPLACEMENT_X)
    #         node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)

    ## setup geotechnical load

    lining_util = LiningUtility()

    load_prop = model.model_part.Properties[1]

    load_prop.SetValue(INTEGRATION_ORDER, 3)
    load_prop.SetValue(G_CONSTANT, g)
    load_prop.SetValue(K0, k0_value)
    load_prop.SetValue(SOIL_DENSITY, 0.017 / g)

    for cond_id in model.layer_cond_sets['outer']:
        cond = model.model_part.Conditions[cond_id]
        cond.Initialize(model.model_part.ProcessInfo)
        lining_util.SetLocalFrame(cond, cy, cz, model.model_part.ProcessInfo)
        cond.SetValue(VERTICAL_LOAD, sigma_v)
        cond.SetValue(HORIZONTAL_LOAD, sigma_h)

    ## setup hyperstatic bedding
    bedding_prop = model.model_part.Properties[2]
    bedding_prop.SetValue(INTEGRATION_ORDER, 3)

    sample_cond_name = "HyperstaticBeddingSecondOrder3D9N"
    last_cond_id = len(model.model_part.Conditions) + 1

    restricting_axis = Vector(3)
    restricting_axis[0] = 1.0
    restricting_axis[1] = 0.0
    restricting_axis[2] = 0.0

    cosf = math.cos(friction_angle)
    sinf = math.sin(friction_angle)
    pclim = 2*cohesion*cosf/(1.0-sinf)+(1.0+sinf)/(1.0-sinf)*sigma_v*(1.0+k0_value)/2.0*(poisson_ratio/(1.0-poisson_ratio))
    print("pclim", pclim)
    f = math.sqrt(sigma_v**2 + sigma_h**2) # good one #1
    #f = sigma_v*cosf + sigma_h*sinf # good one #2
    # f = math.sqrt(sigma_v**2 + sigma_h**2 + (poisson_ratio*sigma_h)**2) # good one #3
    hrm_calib_factor = 2*f/(pclim*math.log(pclim**2 / (pclim - f)) + f)
    #hrm_calib_factor = 0.85
    #hrm_calib_factor = 1.0
    print("hrm_calib_factor", hrm_calib_factor)
    bedding_prop.SetValue(HRM_CALIBRATION_FACTOR, hrm_calib_factor)

    new_bedding_conds = []
    reverse = False
    for cond_id in model.layer_cond_sets['outer']:
        cond = model.model_part.Conditions[cond_id]
        new_cond = lining_util.CreateCondition(cond, sample_cond_name, last_cond_id, bedding_prop)
        new_cond.Initialize(model.model_part.ProcessInfo)
        #lining_util.SetLocalFrame(new_cond, cy, cz, model.model_part.ProcessInfo)
        lining_util.SetLocalFrame(new_cond, model.model_part.ProcessInfo, reverse)
    #    lining_util.SetHyperstaticParams(new_cond, ground_modulus, ground_density, g, poisson_ratio, tunnel_radius, beta_factor, tangent_spring_factor, cohesion, friction_angle, k0_value, model.model_part.ProcessInfo)
        lining_util.SetHyperstaticParams(new_cond, ground_modulus, poisson_ratio, tunnel_radius, beta_factor, tangent_spring_factor, cohesion, friction_angle, sigma_v, sigma_h, model.model_part.ProcessInfo)
        # new_cond.SetValue(RESTRICTING_AXIS, restricting_axis)
        last_cond_id = last_cond_id + 1
        new_bedding_conds.append(new_cond)

    for cond in new_bedding_conds:
        model.model_part.AddCondition(cond)

    ## static solution
    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    ## post processing
    for node in model.model_part.Nodes:
        if abs(node.Z0 - (cz-tunnel_radius)) < tol:
            bottom_node = node
        if abs(node.Z0 - (cz+tunnel_radius)) < tol:
            top_node = node

    print("displacement at top: %.16e" % top_node.GetSolutionStepValue(DISPLACEMENT_Z))
    print("displacement at bottom: %.16e" % bottom_node.GetSolutionStepValue(DISPLACEMENT_Z))

    if logging:

        # transfer the stress to node
        variable_transfer_util = VariableTransferUtility(MKLPardisoSolver())
        variable_transfer_util.TransferVariablesToNodes(model.model_part, STRESSES)

        n_div = 320
        # compute the lining forces and moment
        import numpy
        lining_post_util = LiningPostUtility_ReferenceFrame()
        integration_method = 2
        values = {}
        shift_angle = -math.pi/2
        #angles = []
        #normal_forces = []
        x_axis = Vector(3)
        x_axis[0] = 1.0
        x_axis[1] = 0.0
        x_axis[2] = 0.0
        for i in range(0, n_div):
            cond_id = model.layer_cond_sets['side_' + str(i+1)][0]
            cond = model.model_part.Conditions[cond_id]
            center = lining_post_util.ComputeCenter(cond)
        #    print(cond.Id)
        #    print("forces", str(forces))
        #    print("moment", str(moment))
            y = center[1] - cy
            z = center[2] - cz
        #    print("center", y, z)
            angle = math.atan2(z, y)
        #    angles.append(angle)
            forces = [0.0, 0.0, 0.0]
            moment = [0.0, 0.0, 0.0]
            center[1] = (tunnel_radius-0.5*thickness)*math.cos(angle) + cy
            center[2] = (tunnel_radius-0.5*thickness)*math.sin(angle) + cz
            for cond_id in model.layer_cond_sets['side_' + str(i+1)]:
                cond = model.model_part.Conditions[cond_id]
                cond_forces = lining_post_util.IntegrateForces(cond, integration_method)
                forces[0] = forces[0] + cond_forces[0]
                forces[1] = forces[1] + cond_forces[1]
                forces[2] = forces[2] + cond_forces[2]
                cond_moment = lining_post_util.IntegrateMoment(cond, center, x_axis, integration_method)
                moment[0] = moment[0] + cond_moment[0]
                moment[1] = moment[1] + cond_moment[1]
                moment[2] = moment[2] + cond_moment[2]
            shear_force = (forces[1]*math.cos(angle) + forces[2]*math.sin(angle))
            normal_force = (-forces[1]*math.sin(angle) + forces[2]*math.cos(angle))
        #    normal_forces.append(normal_force)
        #    print("angle", (angle)/math.pi*180.0)
        #    print("shear_force", shear_force)
        #    print("normal_force", normal_force)
            values[cond.Id] = [(angle + shift_angle)%(2*math.pi), normal_force, shear_force, moment[0]]
        #    print("##################################################")
        #print("angles:", angles)
        #print("normal_forces:", normal_forces)

        #export the lining forces
        ifile = open('forces_profile.txt', 'w')
        ifile.write("cond_id\tangle\tnormal_force\tshear_force\tmoment\n")
        values = sorted(values.items(), key = lambda x: x[1][0])
        first_value = copy.deepcopy(values[0])
        first_value[1][0] = 2*math.pi
        values.append(first_value)
        #print(values)
        for value in values:
            ifile.write(str(value[0]))
            ifile.write('\t' + str(value[1][0]))
            ifile.write('\t' + str(value[1][1]))
            ifile.write('\t' + str(value[1][2]))
            ifile.write('\t' + str(value[1][3]))
            ifile.write('\n')
        ifile.close()

        # check the element stress-x
        min_sigma = [1.0e99, 1.0e99, 1.0e99, 1.0e99, 1.0e99, 1.0e99]
        max_sigma = [-1.0e99, -1.0e99, -1.0e99, -1.0e99, -1.0e99, -1.0e99]
        for elem in model.model_part.Elements:
            stresses = elem.GetValuesOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
            # print("stresses:", stresses)
            for stress in stresses:
                for i in range(0, 6):
                    if stress[i] < min_sigma[i]:
                        min_sigma[i] = stress[i]

                    if stress[i] > max_sigma[i]:
                        max_sigma[i] = stress[0]

        print("min_sigma_xx: " + str(min_sigma[0]))
        print("min_sigma_xy: " + str(min_sigma[3]))
        print("min_sigma_xz: " + str(min_sigma[5]))
        print("max_sigma_xx: " + str(max_sigma[0]))
        print("max_sigma_xy: " + str(max_sigma[3]))
        print("max_sigma_xz: " + str(max_sigma[5]))

        print("min_sigma_yy: " + str(min_sigma[1]))
        print("min_sigma_zz: " + str(min_sigma[2]))
        print("max_sigma_yy: " + str(max_sigma[1]))
        print("max_sigma_zz: " + str(max_sigma[2]))

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.Z0 - (cz-tunnel_radius)) < tol:
            bottom_node = node
        if abs(node.Z0 - (cz+tunnel_radius)) < tol:
            top_node = node

    u_top = top_node.GetSolutionStepValue(DISPLACEMENT_Z)
    u_bottom = bottom_node.GetSolutionStepValue(DISPLACEMENT_Z)

    ######### pytesting results #########
    u_top_ref = -6.5574838534407202e-03
    u_bottom_ref = 6.5574838534407081e-03
    assert(abs(u_top - u_top_ref) < 1e-10)
    assert(abs(u_bottom - u_bottom_ref) < 1e-10)
    #####################################
    print("Test passed")

def tag():
    tags = "vhrm,solid"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
