##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
import math
import time as time_module
start_time = time_module.time()
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MortarApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.IsogeometricMortarApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
kernel = Kernel()   #defining kernel

import model_iga_include
import mortar_gpts_element_based_active_set_penalty_contact_strategy

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
mpatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
mortar_util = MortarUtility()

import geometry_factory

x_start = 2.0

def CreateMultiPatch():
    ## create ring - penetrator

    # part 1 of the ring
    axis = 'z'
    r = 1.0
    cen = [x_start, 3.1, 0.0]

    slave_patches = geometry_factory.CreateHalfCircle4(cen, axis, r, 180.0)

    master_patch_ptr = geometry_factory.CreateRectangle([0.0, 0.0, 0.0], [10.0, 2.0, 0.0])

    ######create multipatch
    cnt = 0
    mpatch = MultiPatch2D()
    for patch_ptr in slave_patches:
        patch = patch_ptr.GetReference()
        mpatch.AddPatch(patch_ptr)
        cnt = cnt + 1
        patch.Id = cnt
        patch.LayerIndex = 1
    mpatch.AddPatch(master_patch_ptr)
    master_patch = master_patch_ptr.GetReference()
    master_patch.Id = 5
    master_patch.LayerIndex = 2

    # # elevate the degree
    multipatch_refine_util.DegreeElevate(mpatch[1], [0, 1, 0])
    multipatch_refine_util.DegreeElevate(mpatch[5], [1, 1, 0])

    mpatch.Enumerate()
    #    print(mpatch)
    mpatch_export.Export(mpatch, "ironing.m")

    print("##################")

    return mpatch

def RefineSlave(mpatch, ins_knots_u, ins_knots_v):
    print("###############SLAVE REFINEMENT###############")
    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots_u, ins_knots_v])
    multipatch_refine_util.InsertKnots(mpatch[2], [ins_knots_u, []])
    return mpatch

def RefineMaster(mpatch, ins_knots_u, ins_knots_v):
    print("###############MASTER REFINEMENT###############")
    multipatch_refine_util.InsertKnots(mpatch[5], [ins_knots_u, ins_knots_v])
    return mpatch

def CreateModelPart(mpatch):
    mpatch_util = MultiPatchUtility()
    element_name = "KinematicLinearBezier2D"
#    element_name = "TotalLagrangianBezier2D"
    mortar_condition_name = "LineMortarConditionBezier2D"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( model_part )
    model_part.AddNodalSolutionStepVariable(CONTACT_FORCE)
    model_part.AddNodalSolutionStepVariable(CONTACT_PRESSURE)
    model_part.AddNodalSolutionStepVariable(GAP)
    model_part.AddNodalSolutionStepVariable(THREED_STRESSES)
    model_part.AddNodalSolutionStepVariable(TANGENTIAL_STRESS)

    mpatch_mp.CreateNodes()

    #problem data
    body_force = ZeroVector(3)
    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0 #-9.81
    gravity[2] = 0.0
    ring_prop = model_part.Properties[1]
    ring_prop.SetValue(GRAVITY, gravity )
    ring_prop.SetValue(BODY_FORCE, body_force )
    ring_prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    ring_prop.SetValue(INTEGRATION_ORDER, 2)
    ring_prop.SetValue(DENSITY,            1.0 )
    ring_prop.SetValue(YOUNG_MODULUS, 3.0e4 )
    ring_prop.SetValue(POISSON_RATIO, 0.3 )
    ring_prop.SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
    ring_prop.SetValue(THICKNESS, 1)

    foundation_prop = model_part.Properties[2]
    foundation_prop.SetValue(GRAVITY, gravity )
    foundation_prop.SetValue(BODY_FORCE, body_force )
    foundation_prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    foundation_prop.SetValue(INTEGRATION_ORDER, 2)
    foundation_prop.SetValue(DENSITY,            1.0 )
    foundation_prop.SetValue(YOUNG_MODULUS, 1.0e4 )
    foundation_prop.SetValue(POISSON_RATIO, 0.3 )
    foundation_prop.SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
    foundation_prop.SetValue(THICKNESS, 1)

    contact_indices = IntegerVector(1)
    contact_indices[0] = 10

    print("################################################")
    print("############ GENERATE MODEL_PART ###############")
    print("################################################")

    patch_ids = [1, 2, 3]
    ref_point = Point3D()
    ref_point[0] = x_start
    ref_point[1] = 3.1
    ref_point[2] = 0.0
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, ring_prop)

        ## add slave mortar conditions
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        slave_mortar_conds = mpatch_mp.AddConditions(patch, BoundarySide2D.V0, mortar_condition_name, last_cond_id+1, ring_prop)
#        print("len(slave_mortar_conds):", len(slave_mortar_conds))
        for cond in slave_mortar_conds:
            cond.SetValue(SLAVE_INDEX_SET, contact_indices)
            cond.Initialize(model_part.ProcessInfo)
#            print("condition nodes:")
#            for node in cond.GetNodes():
#                print(" ", node.Id, ": ", node.X0, " ", node.Y0, " ", node.Z0)
#            print("condition integration_points:")
#            mortar_util.DumpIntegrationPoints(cond)
#            print("Check normal for slave condition " + str(cond.Id))
            i = mortar_util.CheckNormalReverse(cond, ref_point, False)
            if i != -1:
                raise Exception("WARNING: The slave normal does not point outward at integration point " + str(i) + " of condition " + str(cond.Id) + ", patch " + str(sid))

    patch_ids = [4]
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, ring_prop)

    patch_ids = [5]
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, foundation_prop)

        ## add master mortar conditions
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        bpatch_ptr = patch.ConstructBoundaryPatch(BoundarySide2D.V1)
        bpatch = bpatch_ptr.GetReference()
        bpatch.Id = 100
        bsplines_patch_util.Reverse(bpatch, 0)
        # print(bpatch)
        # sys.exit(0)
        master_mortar_conds = mpatch_mp.AddConditions(bpatch, mortar_condition_name, last_cond_id+1, foundation_prop)
        for cond in master_mortar_conds:
            cond.SetValue(MASTER_INDEX_SET, contact_indices)
            cond.Initialize(model_part.ProcessInfo)
            # print("Check normal for master condition " + str(cond.Id))
            i = mortar_util.CheckNormal(cond, ref_point, False)
            if i != -1:
                raise Exception("WARNING: The master normal does not point outward at integration point " + str(i))

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    return mpatch_mp

def ExtractLayer(mpatch):
    layer_nodes_sets = {}
    layer_nodes_sets['foundation_bottom'] = []
    layer_nodes_sets['ring_top'] = []

    ring1_ptr = mpatch[1]
    ring1 = ring1_ptr.GetReference()
    ring1_top = ring1.FESpace().BoundaryFunctionIndices(BoundarySide2D.U0)
    for eq_id in ring1_top:
        layer_nodes_sets['ring_top'].append(eq_id + 1)

    ring3_ptr = mpatch[3]
    ring3 = ring3_ptr.GetReference()
    ring3_top = ring3.FESpace().BoundaryFunctionIndices(BoundarySide2D.U1)
    for eq_id in ring3_top:
        layer_nodes_sets['ring_top'].append(eq_id + 1)

    ring4_ptr = mpatch[4]
    ring4 = ring4_ptr.GetReference()
    ring4_top = ring4.FESpace().BoundaryFunctionIndices(BoundarySide2D.V1)
    for eq_id in ring4_top:
        layer_nodes_sets['ring_top'].append(eq_id + 1)

    foundation_ptr = mpatch[5]
    foundation = foundation_ptr.GetReference()
    foundation_top = foundation.FESpace().BoundaryFunctionIndices(BoundarySide2D.V0)
    for eq_id in foundation_top:
        layer_nodes_sets['foundation_bottom'].append(eq_id + 1)

    return layer_nodes_sets

class WriteOutputHelper:
    def __init__(self, mpatch_mp, dim, params, model_part):
        self.mpatch_mp = mpatch_mp
        self.mpatch = self.mpatch_mp.GetMultiPatch()
        self.dim = dim
        self.params = params
        self.model_part = model_part

    def WriteOutput(self, time):
        self.mpatch_mp.SynchronizeBackward(DISPLACEMENT)
        self.mpatch_mp.SynchronizeBackward(THREED_STRESSES)
        import model_iga_include
        post_model_part = model_iga_include.CreatePostModelPart(self.mpatch, self.dim, self.params)
        # # transfer the mortar links from model_part to post_model_part
        mortar_util.TransferMortarIPLink(self.model_part, post_model_part, "LineMortarCondition3D2N", post_model_part.Properties[2])
        # print(post_model_part)
        model_iga_include.WriteGiD(post_model_part, time, self.params)

def CreateModel(logging=True):
    mpatch = CreateMultiPatch()

    ins_knots_u = []
    nsampling_u = 80
    for i in range(1, nsampling_u):
        ins_knots_u.append(float(i)/nsampling_u)
    ins_knots_v = []
    nsampling_v = 5
    for i in range(1, nsampling_v):
        ins_knots_v.append(float(i)/nsampling_v)
    mpatch = RefineSlave(mpatch, ins_knots_u, ins_knots_v)

    ins_knots_u = []
    nsampling_u = 40
    for i in range(1, nsampling_u):
        u = float(i)/nsampling_u
        ins_knots_u.append((math.exp(u)-1) / (math.exp(1)-1))
    ins_knots_v = []
    nsampling_v = 5
    for i in range(1, nsampling_v):
        ins_knots_v.append(float(i)/nsampling_v)
    mpatch = RefineMaster(mpatch, ins_knots_u, ins_knots_v)

    mpatch.Enumerate()
    #    print(mpatch)
    mpatch_export.Export(mpatch, "ironing_refined.m")
    # sys.exit(0)

    mpatch_mp = CreateModelPart(mpatch)

    layer_nodes_sets = ExtractLayer(mpatch)
#    sys.exit(0)

    #############ANALYSIS MODEL#######################################
    model_part = mpatch_mp.GetModelPart()
    params = model_iga_include.StaticParameters()
    params["builder_and_solver_type"] = "residual-based block"
    params["log_residuum"] = logging
    model = model_iga_include.Model('ironing', os.getcwd()+"/", model_part, params)
    # change to contact solver
    model.solver = mortar_gpts_element_based_active_set_penalty_contact_strategy.SampleSolver(model_part, model.abs_tol, model.rel_tol, model.analysis_parameters)
    model.solver.structure_linear_solver = MKLPardisoSolver()
    model.solver.Initialize()
    (model.solver.solver).SetEchoLevel(model.analysis_parameters['echo_level'])
    (model.solver.solver).max_iter = model.analysis_parameters['max_iter'] #control the maximum iterations of Newton Raphson loop
    (model.solver.solver).CalculateReactionsFlag = False
    (model.solver.solver).MoveMeshFlag = False
    model.InitializeModel()

    ## boundary condition
    for node_id in layer_nodes_sets['foundation_bottom']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)

    # tol = 1.0e-6
    # for node in model.model_part.Nodes:
    #     if abs(node.X0 - 0.0) < tol:
    #         node.Fix(DISPLACEMENT_X)

    prescribed_nodes = []
    for node_id in layer_nodes_sets['ring_top']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        prescribed_nodes.append(node)

    ## setup contact
    # define some parameters for contact solution strategy
    model.solver.solver.contact_tying_indices = {}
    # model.solver.solver.contact_tying_indices[10] = Isogeometric_Contact_Link_Kinematic_Linear_Penalty_Regularized_Coulomb_Friction_NonSym_2D()
    model.solver.solver.contact_tying_indices[10] = Isogeometric_Contact_Link_Kinematic_Linear_Penalty_Regularized_Coulomb_Friction_2D()
    model.solver.solver.Parameters['penalty'] = {10: 1.0e7}
    model.solver.solver.Parameters['penalty_t'] = {10: 4.0e5}
    model.solver.solver.Parameters['friction_coefficient'] = {10: 0.3}
    model.solver.solver.Parameters['integration_type'] = "element based"
    model.solver.solver.Parameters['dimension'] = 2
    model.solver.solver.Parameters['gap_tolerance'] = 1.0e-10
    model.solver.solver.Parameters['gap_threshold'] = 1.0e99
    model.solver.solver.Parameters['tying_integration_order'] = 2
    model.solver.solver.Parameters['predict_local_point_method'] = 0
    model.solver.solver.Parameters['maximal_detection_distance'] = 1e-2
    model.solver.solver.Parameters['polygon_offset'] = 0.01
    model.solver.solver.Parameters['compute_min_max_gap'] = True
    model.solver.solver.Parameters['decouple_build_and_solve'] = False
    model.solver.solver.Parameters['max_newton_raphson_iter'] = 20
    model.solver.solver.Parameters['stop_Newton_Raphson_if_not_converged'] = True
    model.solver.solver.Parameters['max_active_set_iter'] = 5
    model.solver.solver.Parameters['stop_active_set_if_not_converged'] = False
    model.solver.solver.Parameters['mortar_echo_level'] = 1 #3
    model.solver.solver.Parameters['transfer_tangential_stress'] = False
    model.solver.solver.Parameters['transfer_old_tangential_stress'] = True
    model.solver.solver.Parameters['transfer_old_master'] = True
    # model.solver.solver.Parameters['variable_transfer_linear_solver'] = SuperLUSolver()
    query_tool = IsogeometricBVHProjectionQueryTool2D()
    query_tool.SetValue(NUM_DIVISION_1, 5)
    model.solver.solver.Parameters['query_tool'] = query_tool
    #model.solver.solver.Parameters['test_linearization'] = False
    #model.solver.solver.Parameters['test_linearization_disp'] = 1.0e-7
    #model.solver.solver.Parameters['test_linearization_tol'] = 1.0e-6
    #model.solver.solver.Parameters['visualize_contact_pairs'] = False
    #model.solver.solver.Parameters['compute_contact_force'] = True
    #model.solver.solver.Parameters['active_set_rel_tol'] = 1.0e-8
    #model.solver.solver.Parameters['active_set_abs_tol'] = 1.0e-12
    #model.solver.solver.Parameters['max_active_set_iter'] = 10
    model.solver.solver.InitializeContact()

    return [mpatch_mp, model, layer_nodes_sets, prescribed_nodes]

def main(logging=True, output=True, npush=64, nslide=600):

    # adolc_util = AdolCMortarTapeUtility()
    # adolc_util.Register_Tape_Projection()
    # adolc_util.Register_Tape_ContactLinkKinematicLinearPenaltyContact()

    ## post processing parameters
    params_post = {}
    params_post['name'] = "ironing"
#    params_post['division mode'] = "uniform"
#    params_post['uniform division number'] = 20
    ##
    # params_post['division mode'] = "non-uniform"
    # params_post['division number u'] = 40
    # params_post['division number v'] = 5
    ##
    params_post['division mode'] = "non-uniform per patch"
    params_post['division number u'] = {1: 15, 2: 15, 3: 15, 4: 15, 5: 80}
    params_post['division number v'] = {1: 3, 2: 3, 3: 3, 4: 3, 5: 10}
    params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES]
    params_post['output format'] = "binary" # "ascii"
    dim = 2

    #############ANALYSIS MODEL#######################################
    [mpatch_mp, model, layer_nodes_sets, prescribed_nodes] = CreateModel()
    mpatch = mpatch_mp.GetMultiPatch()

    write_output_helper = WriteOutputHelper(mpatch_mp, dim, params_post, model.model_part)
    # model.solver.solver.write_output_callback = write_output_helper.WriteOutput

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)

    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())

    # write_output_helper.WriteOutput(time)
    ##################################################################

    print("#####################################################")
    print("#####################################################")
    print("###########SIMULATION################################")
    print("#####################################################")
    print("#####################################################")

    #################################

    disp = 0.0
    delta_disp = 0.1
    delta_time = delta_disp

    disp = disp + delta_disp
    time = time + delta_time
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Y, -disp)

    model.Solve(time, 0, 0, 0, 0)
    if output:
        write_output_helper.WriteOutput(time)
    #################################

    delta_disp = 0.01
    delta_time = delta_disp

    for step in range(0, npush): #64): #100
        print("########PUSHING STEP " + str(step+1) + " STARTS##########")
        disp = disp + delta_disp
        time = time + delta_time
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, -disp)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            write_output_helper.WriteOutput(time)

    time = time + 100.0
    disp = 0.0
    delta_disp = 0.01
    delta_time = 0.01

    for step in range(0, nslide): #600):
        print("########SLIDING STEP " + str(step+1) + " STARTS##########")
        disp = disp + delta_disp
        time = time + delta_time
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)
        model.Solve(time, 0, 0, 0, 0)
        if output:
            write_output_helper.WriteOutput(time)

    print("Analysis completed")
    timer = Timer()
    print(timer)

    return model

def test():
    model = main(logging=False, output=False, npush=10, nslide=20)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.Y0 - 2.1) < tol):
            monitoring_node = node # choose last node

    ux = monitoring_node.GetSolutionStepValue(DISPLACEMENT_X)
    uy = monitoring_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("ux: %.16e, uy: %.16e" % (ux, uy))
    ref_ux = 1.8595303002138053e-01
    ref_uy = -1.7969148618261876e-01
    assert(abs(ux - ref_ux) < 1e-10)
    assert(abs(uy - ref_uy) < 1e-10)

    print("Test passed")

def tag():
    return "contact,friction"

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
