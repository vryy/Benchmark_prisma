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
import operator
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
import pdb
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.MKLSolversApplication import *
#from KratosMultiphysics.MultigridSolversApplication import *
#from KratosMultiphysics.FiniteCellApplication import *
kernel = Kernel()   #defining kernel

import model_iga_include
from model_iga_include import *

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
hbsplines_patch_util = HBSplinesPatchUtility()
hbsplines_refinement_util = HBSplinesRefinementUtility()
hmpatch_export = MultiHBSplinesPatchMatlabExporter()

import geometry_factory

sys.path.append("../../../")
import analytical_solution
# this file can also be found in ${BENCHMARK_PRISMA}/structural_application/std_problems/internal_pressurized_cylinder
# where ${BENCHMARK_PRISMA} is the path of the benchmarking folder, which can be cloned from https://github.com/vryy/Benchmark_prisma.git

E = 2.1e5
nu = 0.3
P = 100.0
r1 = 100.0
r2 = 200.0
ana_sol = analytical_solution.Solution(r1, r2, P, 0.0, E, nu)

def CreateMultiPatch():
    ## create arc 1
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 0.0, 90.0)
    arc1 = arc1_ptr.GetReference()
    arc1.Id = 1

    ## create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r2, 0.0, 90.0)
    arc2 = arc2_ptr.GetReference()
    arc2.Id = 2

    ## create ring patch by connect the two arcs
    ring_patch_ptr = bsplines_patch_util.CreateLoftPatch(arc2, arc1)
    ring_patch = ring_patch_ptr.GetReference()
    ring_patch.Id = 1

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(ring_patch_ptr)

    #### elevate the degree
    multipatch_refine_util.DegreeElevate(mpatch[1], [0, 1])

    return mpatch

def Refine(mpatch, ins_knots):
    print("###############REFINEMENT###############")
    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots, ins_knots])

    return mpatch

def CreateHBMultiPatch(mpatch):
    # convert the NURBS patches to HB-Splines patches
    hmpatch = MultiPatch2D()
    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        hpatch_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch)
        hpatch = hpatch_ptr.GetReference()
        hpatch.Prefix = "HB_" + patch.Prefix
        hmpatch.AddPatch(hpatch_ptr)

    # generate the internal cells
    for hpatch_ptr in hmpatch.Patches():
        hpatch = hpatch_ptr.GetReference()
        hpatch.FESpace().UpdateCells()
#    print(hmpatch)

    return hmpatch

def ExtractLayer(hmpatch):
    layer_nodes_sets = {}
    layer_nodes_sets['left'] = []
    layer_nodes_sets['bottom'] = []
    layer_nodes_sets['inner'] = []
    layer_nodes_sets['outer'] = []

    hpatch1_ptr = hmpatch[1]
    hpatch1 = hpatch1_ptr.GetReference()
    boundary_basis_1_u0 = hpatch1.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.U0))
    for i in range(0, len(boundary_basis_1_u0)):
        bf = boundary_basis_1_u0[i]
        layer_nodes_sets['bottom'].append(bf.EquationId + 1)
    boundary_basis_1_u1 = hpatch1.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.U1))
    for i in range(0, len(boundary_basis_1_u1)):
        bf = boundary_basis_1_u1[i]
        layer_nodes_sets['left'].append(bf.EquationId + 1)
    boundary_basis_1_v0 = hpatch1.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.V0))
    for i in range(0, len(boundary_basis_1_v0)):
        bf = boundary_basis_1_v0[i]
        layer_nodes_sets['outer'].append(bf.EquationId + 1)
    boundary_basis_1_v1 = hpatch1.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.V1))
    for i in range(0, len(boundary_basis_1_v1)):
        bf = boundary_basis_1_v1[i]
        layer_nodes_sets['inner'].append(bf.EquationId + 1)

    return layer_nodes_sets

def CreateModel(mpatch, layer_nodes_sets):
    mpatch_util = MultiPatchUtility()
    element_name = "KinematicLinearBezier2D"
    load_condition_name = "LinePressureBezier2D"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( model_part )
    model_part.AddNodalSolutionStepVariable(THREED_STRESSES)

    mpatch_mp.CreateNodes()

    #problem data
    body_force = ZeroVector(2)
    gravity = Vector(2)
    gravity[0] = 0.0
    gravity[1] = 0.0
    prop = model_part.Properties[1]
    prop.SetValue(GRAVITY, gravity )
    prop.SetValue(BODY_FORCE, body_force )
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 1)
    prop.SetValue(DENSITY,            0 )
    prop.SetValue(YOUNG_MODULUS, E )
    prop.SetValue(POISSON_RATIO, nu )
    prop.SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
    prop.SetValue(THICKNESS, 1)

    patch_ids = [1]
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)

        ## add loading conditions
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        print("last_cond_id", last_cond_id)
        load_conds = mpatch_mp.AddConditions(patch, BoundarySide2D.V1, load_condition_name, last_cond_id+1, prop)
        print("load_conds:", load_conds)
        for cond in load_conds:
            cond.SetValue(PRESSURE, -P)

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    # fix displacement on the bottom
    for node_id in layer_nodes_sets['bottom']:
        node = model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_Y)
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

    # fix displacement on the left
    for node_id in layer_nodes_sets['left']:
        node = model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_X)
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)

    return mpatch_mp

def Solve(model_part, time, logging=True):
    #############ANALYSIS MODEL#######################################
    params = model_iga_include.StaticParameters()
    params["builder_and_solver_type"] = "residual-based block"
    params["log_residuum"] = logging
    model = model_iga_include.Model('internal_pressurized_cylinder', os.getcwd()+"/", model_part, params)
    model.InitializeModel()
    model.Solve(time, 0, 0, 0, 0)

    # transfer back the results to node
    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())

def main(logging=True, output=True):
    ## create the B-Splines multipatch
    mpatch = CreateMultiPatch()
#    mpatch = Refine(mpatch, [0.5])
    mpatch.Enumerate()
    #    print(mpatch)
    if logging:
        mpatch_export.Export(mpatch, "internal_pressurized_cylinder.m")

    ## create the hierarchical B-Splines multipatch
    hmpatch = CreateHBMultiPatch(mpatch)

    ## extract layer
    layer_nodes_sets = ExtractLayer(hmpatch)
    print("layer_nodes_sets:", layer_nodes_sets)
#    pdb.set_trace()

    ###############################
    params_post = {}
    params_post['name'] = "internal_pressurized_cylinder"
    params_post['division mode'] = "uniform"
    params_post['uniform division number'] = 40
    #params_post['division mode'] = "non-uniform"
    #params_post['division number u'] = 10
    #params_post['division number v'] = 10
    #params_post['division number w'] = 1
    params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES]
    dim = 2
    ###############################

    nsteps = 3
    time = 0.0
    delta_time = 1.0
    l2_errors = []
    h1_errors = []
    for i in range (0, nsteps):

        #############ANALYSIS MODEL#######################################

        hmpatch_mp = CreateModel(hmpatch, layer_nodes_sets)

        model_part = hmpatch_mp.GetModelPart()

        time = time + delta_time
        if logging:
            hmpatch_export.Export(hmpatch, "hb_internal_pressurized_cylinder_" + str(int(time)) + ".m")
        Solve(model_part, time, logging=logging)

        ######Synchronize back the results to multipatch
        hmpatch_mp.SynchronizeBackward(DISPLACEMENT)
        hmpatch_mp.SynchronizeBackward(THREED_STRESSES)
        ##################################################################

        if output:
            ## post processing
            model_iga_include.PostMultiPatch(hmpatch, dim, time, params_post)

        ################# COMPUTE ERROR ###############################

        hpatch1_ptr = hmpatch[1]
        hpatch1 = hpatch1_ptr.GetReference()
        echo_level = IsogeometricEchoFlags.ECHO_REFINEMENT + IsogeometricEchoFlags.ECHO_REFINEMENT_DETAIL

        ## compute error
        print("********Error Analysis at time " + str(time) + "*******")
        print("L2Error\tH1Error")
        H1ErrorList = {}
        for elem in model_part.Elements:
            L2Error = model_iga_include.ComputeL2errorOnElement(elem, ana_sol, model_part.ProcessInfo)
            H1Error = model_iga_include.ComputeH1errorOnElement(elem, ana_sol, model_part.ProcessInfo)
            H1ErrorList[elem.Id] = H1Error
            print("%.6e\t%.6e" % (L2Error, H1Error))
        l2_errors.append( model_iga_include.ComputeL2error(model_part, ana_sol) )
        h1_errors.append( model_iga_include.ComputeH1error(model_part, ana_sol) )
        print("*******************************************")
        print ("H1ErrorList", H1ErrorList)

        ################# HIERACHICAL REFINE ###############################

        if i == nsteps-1:
            break

        H1ErrorList_sorted = sorted(H1ErrorList.items(), key = operator.itemgetter(1), reverse=True)
        print("H1ErrorList_sorted:", H1ErrorList_sorted)

        refine_rate = 0.2
        refined_elems_count = int(refine_rate*len(model_part.Elements))
        if refined_elems_count == 0:
            refined_elems_count = 1

        for j in range(0, refined_elems_count):
            tmp = H1ErrorList_sorted[j]
            elem_id = tmp[0]
            elem = model_part.Elements[elem_id]
            left = elem.GetValue(KNOT_LEFT)
            right = elem.GetValue(KNOT_RIGHT)
            top = elem.GetValue(KNOT_TOP)
            bottom = elem.GetValue(KNOT_BOTTOM)
            print("Refine window according to element", elem_id)
            print("Left\tRight\tTop\tBottom")
            print(left, right, top, bottom)
            hbsplines_refinement_util.RefineWindow(hpatch1, [[left, right], [bottom, top]], echo_level)
            hbsplines_refinement_util.LinearDependencyRefine(hpatch1, 0, echo_level)

        # generate the internal cells
        for hpatch_ptr in hmpatch.Patches():
            hpatch = hpatch_ptr.GetReference()
            hpatch.FESpace().UpdateCells()

        ## extract layer
        layer_nodes_sets = ExtractLayer(hmpatch)
        print("layer_nodes_sets:", layer_nodes_sets)

    return l2_errors, h1_errors

def test():
    l2_errors, h1_errors = main(logging=False, output=False)

    print("l2_errors:", [f"{x:.16e}" for x in l2_errors])
    print("h1_errors:", [f"{x:.16e}" for x in h1_errors])

    l2_error_ref = 2.6917235390884519e-03
    h1_error_ref = 3.7016838456000387e-02

    assert(abs(l2_errors[-1] - l2_error_ref) < 1e-10)
    assert(abs(h1_errors[-1] - h1_error_ref) < 1e-10)

    print("Test passed")

def tag():
    return "IGA,hbsplines"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
