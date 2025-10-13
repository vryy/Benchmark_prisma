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
##################################################################
##################################################################
import ironing2_include
from ironing2_include import *
##################################################################
###  SIMULATION  #################################################
# This simulation is inspired by the paper:
# Puso et al, A mortar segment-to-segment contact method for large deformation solid mechanics
##################################################################

def main(output=True, logging=True, npush=90, nslide=300):
    model = ironing2_include.Model('ironing2',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    prescribed_nodes = []
    for node_id in model.layer_nodes_sets['load']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        prescribed_nodes.append(node)

    # define some parameters for contact solution strategy
    model.solver.solver.contact_tying_indices = {}
    model.solver.solver.contact_tying_indices[10] = "contact_link_kinematic_linear_augmented_lagrange"
    model.solver.solver.Parameters['penalty'] = {10: 1.0e2}
    model.solver.solver.Parameters['friction_coefficient'] = {10: 0.0}
    model.solver.solver.Parameters['gap_tolerance'] = 1.0e99 # this is required for Augmented Lagrangian to include all proximity segments
    model.solver.solver.Parameters['maximal_detection_distance'] = 1e-3
    model.solver.solver.Parameters['tying_integration_order'] = 3
    model.solver.solver.Parameters['predict_local_point_method'] = 2
    model.solver.solver.Parameters['test_linearization'] = False
    model.solver.solver.Parameters['test_linearization_disp'] = 1.0e-7
    model.solver.solver.Parameters['test_linearization_tol'] = 1.0e-6
    model.solver.solver.Parameters['visualize_contact_pairs'] = False
    model.solver.solver.Parameters['compute_contact_force'] = True
    model.solver.solver.Parameters['compute_min_max_gap'] = True
    model.solver.solver.Parameters['stop_Newton_Raphson_if_not_converged'] = True
    model.solver.solver.Parameters['update_master_in_each_iteration'] = True
    model.solver.solver.InitializeContact()

    # solve zero displacement
    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    #model.WriteOutput(time)

    print("#####################################################")
    print("#####################################################")
    print("###########LOAD STEP START###########################")
    print("#####################################################")
    print("#####################################################")

    disp = 0.0

    nramp = 1
    delta_disp = 0.5 / nramp
    delta_time = delta_disp
    for i in range(0, nramp):
        time = time + delta_time
        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Z, -disp)

        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#####################################################")
        print("########LOAD STEP " + str(time) + " IS SOLVED########")
        print("#####################################################")

    print("#####################################################")
    print("#####################################################")
    print("###########PENETRATION START#########################")
    print("#####################################################")
    print("#####################################################")

    delta_time = 0.01
    delta_disp = 0.01
    for i in range(0, npush):
        time = time + delta_time
        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Z, -disp)

        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#####################################################")
        print("#####PENETRATION STEP " + str(i+1) + " IS SOLVED#####")
        print("#####time: " + str(time) + "#########################")
        print("#####################################################")

    print("#####################################################")
    print("#####################################################")
    print("###############SLIDING START#########################")
    print("#####################################################")
    print("#####################################################")

    time = 1000.0

    disp = 0.0
    delta_time = 0.01
    delta_disp = 0.01
    for i in range(0, nslide):
        time = time + delta_time
        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, disp)

        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#####################################################")
        print("######SLIDING STEP " + str(i+1) + " IS SOLVED########")
        print("######time: " + str(time) + "########################")
        print("#####################################################")

    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False, npush=10, nslide=20)

    monitoring_node = model.model_part.Nodes[1258]
    disp_x = monitoring_node.GetSolutionStepValue(DISPLACEMENT_X)
    disp_z = monitoring_node.GetSolutionStepValue(DISPLACEMENT_Z)
    print("%.16e" % (disp_x))
    print("%.16e" % (disp_z))
    ref_disp_x = 2.0000135799920685e-01
    ref_disp_z = -5.9999509819717622e-01
    assert(abs(disp_x - ref_disp_x) < 1e-10)
    assert(abs(disp_y - ref_disp_z) < 1e-10)

    print("Test passed")

def tag():
    return "contact"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
