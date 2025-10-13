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
sys.path.append('./ironing.gid')
import ironing_include
from ironing_include import *

def main(logging=True, output=True):
    model = ironing_include.Model('ironing',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## boundary conditions
    for node_id in model.layer_nodes_sets['base']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

    prescribed_nodes = []
    for node_id in model.layer_nodes_sets['load']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        prescribed_nodes.append(node)

    mortar_util = MortarUtility()
    master_conds = mortar_util.ExtractMasterMortarEntities(model.model_part.Conditions, 10)
    FindNeighbourConditionsInMortarSegmentsProcess(master_conds).Execute()

    ## contact parameters
    model.solver.solver.contact_tying_indices = {}
    model.solver.solver.contact_tying_indices[10] = "contact_gptms_c0n_link_kinematic_linear_augmented_lagrange"
    model.solver.solver.Parameters['penalty'] = {10: 1.0e4}
    model.solver.solver.Parameters['friction_coefficient'] = {10: 0.0}
    model.solver.solver.Parameters['stop_active_set_if_not_converged'] = True
    model.solver.solver.Parameters['update_master_in_each_iteration'] = True
    model.solver.solver.Parameters['dimension'] = 2
    model.solver.solver.Parameters['gap_tolerance'] = 1.0e99 # this is required for Augmented Lagrangian to include all proximity segments
    model.solver.solver.Parameters['tying_integration_order'] = 4
    model.solver.solver.Parameters['predict_local_point_method'] = 0
    model.solver.solver.Parameters['maximal_detection_distance'] = 1e-3
    model.solver.solver.Parameters['polygon_offset'] = 0.01
    model.solver.solver.Parameters['query_tools'] = {}
    model.solver.solver.Parameters['query_tools'][10] = NodalNormalProjectionQueryTool2D()
    #model.solver.solver.Parameters['test_linearization'] = False
    #model.solver.solver.Parameters['test_linearization_disp'] = 1.0e-7
    #model.solver.solver.Parameters['test_linearization_tol'] = 1.0e-6
    #model.solver.solver.Parameters['visualize_contact_pairs'] = False
    #model.solver.solver.Parameters['compute_contact_force'] = True
    #model.solver.solver.Parameters['active_set_rel_tol'] = 1.0e-8
    #model.solver.solver.Parameters['active_set_abs_tol'] = 1.0e-12
    #model.solver.solver.Parameters['max_active_set_iter'] = 10
    model.solver.solver.InitializeContact()
    # model.solver.solver.write_output_callback = model.WriteOutput

    #################################
    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    #model.WriteOutput(time)

    #################################

    disp = 0.0
    delta_disp = 0.1
    delta_time = 0.1

    disp = disp + delta_disp
    time = time + delta_time
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Y, -disp)

    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    #################################

    print("########PUSHING STARTS##########")

    delta_disp = 0.01
    delta_time = 0.01

    for step in range(0, 64):
        print("#######Pushing Step " + str(step+1) + " starts#######")
        disp = disp + delta_disp
        time = time + delta_time
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, -disp)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#######Pushing Step " + str(step+1) + " completed#######")

    print("########SLIDING STARTS##########")

    time = 100.0
    delta_disp = 0.01
    delta_time = 0.01

    for step in range(0, 400):
        print("#######Sliding Step " + str(step+1) + " starts#######")
        disp = disp + delta_disp
        time = time + delta_time
        for node in prescribed_nodes:
    #        node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#######Sliding Step " + str(step+1) + " completed#######")

    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False)

    monitoring_node = model.model_part.Nodes[545]
    disp_x = monitoring_node.GetSolutionStepValue(DISPLACEMENT_X)
    disp_y = monitoring_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("%.16e" % (disp_x))
    print("%.16e" % (disp_y))
    ref_disp_x = 4.0002216130642481e+00
    ref_disp_y = -6.3433324454162221e-01
    assert(abs(disp_x - ref_disp_x) < 1e-10)
    assert(abs(disp_y - ref_disp_y) < 1e-10)

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
