##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on So 15. Mar 13:24:47 CET 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./ring_and_foundation.gid')
import ring_and_foundation_include
from ring_and_foundation_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = ring_and_foundation_include.Model('ring_and_foundation',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    # fix base nodes
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            #node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            #node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    prescribed_nodes = []
    for node_id in model.layer_nodes_sets['ring_upper_edges_nodes']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)    # restrain x movement of the point
        node.Fix(DISPLACEMENT_Y)
        #node.Fix(DISPLACEMENT_Z)
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)  # 0 is the current step, or the value to be assigned to DISPLACEMENT_X, I'm not sure
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        #node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        prescribed_nodes.append(node)

    # define some parameters for contact solution strategy
    model.solver.solver.contact_tying_indices = {}
    model.solver.solver.contact_tying_indices[10] = "contact_link_kinematic_linear_augmented_lagrange"
    model.solver.solver.Parameters['lambda'] = {10: 0.0}
    model.solver.solver.Parameters['penalty'] = {10: 1.0e4}
    model.solver.solver.Parameters['friction_coefficient'] = {10: 0.0}
    model.solver.solver.Parameters['dimension'] = 2
    model.solver.solver.Parameters['gap_tolerance'] = 1.0e99 # this is required for Augmented Lagrangian to include all proximity segments
    model.solver.solver.Parameters['penetration_tolerance'] = 1.0e-5
    model.solver.solver.Parameters['tying_integration_order'] = 2
    model.solver.solver.Parameters['predict_local_point_method'] = 0
    model.solver.solver.Parameters['maximal_detection_distance'] = 1e-3
    model.solver.solver.Parameters['polygon_offset'] = 0.01
    model.solver.solver.Parameters['compute_min_max_gap'] = True
    model.solver.solver.Parameters['compute_min_max_penalty'] = True
    model.solver.solver.Parameters['max_newton_raphson_iter'] = 30
    model.solver.solver.Parameters['max_uzawa_iter'] = 20
    model.solver.solver.Parameters['update_penalty'] = False
    model.solver.solver.Parameters['stop_uzawa_if_not_converged'] = True
    model.solver.solver.Parameters['stop_Newton_Raphson_if_not_converged'] = True
    model.solver.solver.Parameters['decouple_build_and_solve'] = True
    model.solver.solver.Parameters['mortar_echo_level'] = 1
    #model.solver.solver.Parameters['test_linearization'] = False
    #model.solver.solver.Parameters['test_linearization_disp'] = 1.0e-7
    #model.solver.solver.Parameters['test_linearization_tol'] = 1.0e-6
    #model.solver.solver.Parameters['visualize_contact_pairs'] = False
    model.solver.solver.Parameters['compute_contact_force'] = True
    #model.solver.solver.Parameters['active_set_rel_tol'] = 1.0e-8
    #model.solver.solver.Parameters['active_set_abs_tol'] = 1.0e-12
    #model.solver.solver.Parameters['max_active_set_iter'] = 10
    model.solver.solver.InitializeContact()

    # solve zero displacement
    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    #model.WriteOutput(time)

    print("#####################################################")
    print("#####################################################")
    print("###########First Increment###########################")
    print("#####################################################")
    print("#####################################################")

    # applying y=-20 displacement in 1 step

    time = 0.0
    disp = -20.0
    first_load_num_increment = 20
    delta_time = abs(disp) / first_load_num_increment

    for i in range(0, first_load_num_increment):
        time = time + delta_time
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, (i+1)*disp/first_load_num_increment)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)

    print("#####################################################")
    print("#####################################################")
    print("###########WARNING: HITTING HAPPEN###################")
    print("#####################################################")
    print("#####################################################")

    # applying y=-40 displacement in 80 steps

    #delta_time = 0.5
    #delta_disp = -0.5

    #for i in range(10):
    #    time = time + delta_time
    #    disp = disp + delta_disp
    #    for node in prescribed_nodes:
    #        node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

    #    model.Solve(time, 0, 0, 0, 0)
    #    model.WriteOutput(time)
    #    print("#####################################################")
    #    print("########LOAD STEP " + str(time) + " IS SOLVED########")
    #    print("#####################################################")

    # load: hit the foundation
    delta_time = 0.5
    delta_disp = -0.5
    for i in range(0, 80):
        time = time + delta_time

        print("#####################################################")
        print("########LOAD STEP " + str(time) + " STARTS###########")
        print("#####################################################")

        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#####################################################")
        print("########LOAD STEP " + str(time) + " IS SOLVED########")
        print("#####################################################")
    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False)

    monitoring_node = model.model_part.Nodes[144]
    disp_y = monitoring_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("%.16e" % (disp_y))
    ref_disp_y = -5.8110010716140668e+00
    # ref_disp_y = -5.8110010716138181e+00
    # assert(abs(disp_y - ref_disp_y) < 1e-8)
    assert(abs(disp_y - ref_disp_y) < 1e-3)

    print("Test passed")

def tag():
    return "contact,neo-hookean,ring_and_foundation"

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
