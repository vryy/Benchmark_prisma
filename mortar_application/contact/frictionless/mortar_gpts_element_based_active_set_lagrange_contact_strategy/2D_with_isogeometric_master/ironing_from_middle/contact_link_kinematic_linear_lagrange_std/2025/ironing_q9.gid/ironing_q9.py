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
## This file is generated on Fr 19. Nov 08:16:50 CET 2021
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./ironing_q9.gid')
import ironing_q9_include
from ironing_q9_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
mpatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
mortar_util = MortarUtility()

mesh_util = FiniteCellMeshUtility()

import geometry_factory

def CreateMultiPatch():
    # patch for master line
    master_patch_ptr = geometry_factory.CreateLine([10.0, 2.0, 0.0], [0.0, 2.0, 0.0], order=2)

    ######create multipatch
    cnt = 0
    mpatch = MultiPatch1D()
    mpatch.AddPatch(master_patch_ptr)
    master_patch = master_patch_ptr.GetReference()
    master_patch.Id = 1

    # elevate the degree
    multipatch_refine_util.DegreeElevate(mpatch[1], [1, 0, 0])

    mpatch.Enumerate()

    return mpatch

def RefineMaster(mpatch, nsampling_u=10):
    print("###############MASTER REFINEMENT###############")
    ins_knots_u = []
    for i in range(1, nsampling_u):
        ins_knots_u.append(float(i)/nsampling_u)
    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots_u])
    return mpatch

def CreateModelPart(mpatch):
    mpatch_util = MultiPatchUtility()
    mortar_condition_name = "LineMortarConditionBezier2D"

    mpatch_mp = MultiPatchModelPart1D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( model_part )

    mpatch_mp.CreateNodes()

    #problem data

    contact_indices = IntegerVector(1)
    contact_indices[0] = 10

    print("################################################")
    print("############ GENERATE MODEL_PART ###############")
    print("################################################")

    ref_point = Point3D()
    ref_point[0] = 5.0
    ref_point[1] = 3.0
    ref_point[2] = 0.0

    prop = model_part.Properties[1]

    patch_ids = [1]
    for sid in patch_ids:
        patch_ptr = mpatch[sid]
        patch = patch_ptr.GetReference()

        ## add master mortar conditions
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        patch = patch_ptr.GetReference()
        patch.Id = 100
        master_mortar_conds = mpatch_mp.AddConditions(patch, mortar_condition_name, last_cond_id+1, prop)
        for cond in master_mortar_conds:
            cond.SetValue(MASTER_INDEX_SET, contact_indices)
            cond.Initialize(model_part.ProcessInfo)
            print("Check normal for master condition " + str(cond.Id))
            i = mortar_util.CheckNormal(cond, ref_point, False)
            if i != -1:
                raise Exception("WARNING: The master normal does not point outward at integration point " + str(i))

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    return mpatch_mp

def CreateModel(logging=True):
    model = ironing_q9_include.Model('ironing_q9',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    # create multipatch
    mpatch = CreateMultiPatch()
    mpatch = RefineMaster(mpatch, nsampling_u=40)
    #    print(mpatch)
    if logging:
        mpatch_export.Export(mpatch, "ironing_master.m")

    # create model_part for the master segment
    mpatch_mp = CreateModelPart(mpatch)
    master_model_part = mpatch_mp.GetModelPart()

    imported_nodes = mesh_util.ImportNodes(model.model_part, master_model_part, 0)
    model.AddDofsForNodes(imported_nodes)
    imported_master_conds = mesh_util.ImportConditions(model.model_part, master_model_part.Conditions, 0, 0)
    model.layer_cond_sets['isogeometric_master'] = []
    master_index_set = IntegerVector(2)
    master_index_set[0] = 10 # for contact
    master_index_set[1] = 20 # for tying
    for cond in imported_master_conds:
        model.layer_cond_sets['isogeometric_master'].append(cond.Id)
        cond.SetValue(MASTER_INDEX_SET, master_index_set)

    print(model.layer_cond_sets['isogeometric_master'])

    return model

def main(logging=True, output=True, n_push=64, n_slide=400):

    model = CreateModel(logging=logging)

    # model.WriteOutput(0.0)

    # sys.exit(0)

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

    # setup tying
    tying_util = MortarTyingUtility2D()
    tying_util.SetDefaultSearchParameters()
    tying_util.SetValue(PREDICT_LOCAL_POINT_METHOD, 0)
    tying_util.SetQueryTool(20, IsogeometricProjectionQueryTool2D())
    tying_util.InitializeSearchTree(model.model_part, 20)
    tying_util.SetEchoLevel(1)
    mortar_links = tying_util.SetupTyingLinkElementBased(model.model_part, 20, "tying_link_kinematic_linear_penalty")
    penalty = 1e8
    for cond in mortar_links:
       cond.SetValue(INITIAL_PENALTY, penalty)
    # sys.exit(0)

    ## contact parameters
    model.solver.solver.contact_tying_indices = {}
    model.solver.solver.contact_tying_indices[10] = "contact_link_kinematic_linear_lagrange_std"
    model.solver.solver.Parameters['friction_coefficient'] = {10: 0.0}
    model.solver.solver.Parameters['update_master_in_each_iteration'] = False
    model.solver.solver.Parameters['stop_active_set_if_not_converged'] = True
    model.solver.solver.Parameters['stop_Newton_Raphson_if_not_converged'] = True
    model.solver.solver.Parameters['dimension'] = 2
    model.solver.solver.Parameters['gap_tolerance'] = 1.0e-10
    model.solver.solver.Parameters['tying_integration_order'] = 4
    model.solver.solver.Parameters['slave_integration_order'] = {10: 4}
    model.solver.solver.Parameters['predict_local_point_method'] = 0
    model.solver.solver.Parameters['maximal_detection_distance'] = 1e-3
    model.solver.solver.Parameters['polygon_offset'] = 0.01
    model.solver.solver.Parameters['mortar_echo_level'] = 1
    model.solver.solver.Parameters['max_active_set_iter'] = 10
    model.solver.solver.Parameters['max_newton_raphson_iter'] = 20
    model.solver.solver.Parameters['query_tool'] = IsogeometricProjectionQueryTool2D()
    model.solver.solver.InitializeContact()
    # model.solver.solver.write_output_callback = model.WriteOutput

    #################################
    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    #model.WriteOutput(time)
    # sys.exit(0)

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
    # sys.exit(0)

    #################################

    print("########PUSHING STARTS##########")

    delta_disp = 0.01
    delta_time = 0.01

    for step in range(0, n_push):
        disp = disp + delta_disp
        time = time + delta_time
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, -disp)
            # node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, -delta_disp)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#######Pushing Step " + str(step) + " completed#######")

    print("########PUSHING COMPLETED##########")
    # sys.exit(0)

    print("########SLIDING STARTS##########")

    time = 100.0
    delta_disp = 0.01
    delta_time = 0.01

    for step in range(0, n_slide):
        print("#######Sliding Step " + str(step) + " starts#######")
        disp = disp + delta_disp
        time = time + delta_time
        for node in prescribed_nodes:
    #        node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("#######Sliding Step " + str(step) + " completed#######")

    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False, n_push=64, n_slide=400)

    monitoring_node = model.model_part.Nodes[2013]
    disp_x = monitoring_node.GetSolutionStepValue(DISPLACEMENT_X)
    disp_y = monitoring_node.GetSolutionStepValue(DISPLACEMENT_Y)
    lmbda = monitoring_node.GetSolutionStepValue(LAMBDA)
    print("%.16e" % (disp_x))
    print("%.16e" % (disp_y))
    print("%.16e" % (lmbda))
    ref_disp_x = 4.0001754270634846e+00
    ref_disp_y = -6.3261108492590368e-01
    ref_lmbda = -4.4578903066321736e+03
    assert(abs(disp_x - ref_disp_x) < 1e-10)
    assert(abs(disp_y - ref_disp_y) < 1e-10)
    assert(abs(lmbda - ref_lmbda) / abs(ref_lmbda) < 1e-10)

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

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
