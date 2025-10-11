from __future__ import absolute_import
import sys
import os

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.IsogeometricThermalApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *

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

cx = 1.6
cy = 2.4

def CreateMultiPatch():

    #### create rectagle
    a = 1.0
    rec_ptr = geometry_factory.CreateRectangle([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
    rec = rec_ptr.GetReference()
    rec.Id = 1

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(rec_ptr)

    return mpatch

def Refine(mpatch, order, nsampling):
    print("###############REFINEMENT###############")
    multipatch_refine_util.DegreeElevate(mpatch[1], [order-1, order-1])

    ins_knots = []
    for i in range(1, nsampling):
        ins_knots.append(float(i)/nsampling)

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
    layer_nodes_sets['right'] = []
    layer_nodes_sets['top'] = []
    layer_nodes_sets['bottom'] = []

    for hpatch_ptr in hmpatch.Patches():
        hpatch = hpatch_ptr.GetReference()
    #hpatch1_ptr = hmpatch[1]
    #hpatch1 = hpatch1_ptr.GetReference()
        boundary_basis_1_u0 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.U0))
        for i in range(0, len(boundary_basis_1_u0)):
            bf = boundary_basis_1_u0[i]
            layer_nodes_sets['left'].append(bf.EquationId + 1)
        boundary_basis_1_u1 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.U1))
        for i in range(0, len(boundary_basis_1_u1)):
            bf = boundary_basis_1_u1[i]
            layer_nodes_sets['right'].append(bf.EquationId + 1)
        boundary_basis_1_v0 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.V0))
        for i in range(0, len(boundary_basis_1_v0)):
            bf = boundary_basis_1_v0[i]
            layer_nodes_sets['bottom'].append(bf.EquationId + 1)
        boundary_basis_1_v1 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.V1))
        for i in range(0, len(boundary_basis_1_v1)):
            bf = boundary_basis_1_v1[i]
            layer_nodes_sets['top'].append(bf.EquationId + 1)

    return layer_nodes_sets

def CreateModel(mpatch, layer_nodes_sets):
    mpatch_util = MultiPatchUtility()
    element_name = "LinearPoissonBezier2D"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)

    mpatch_mp.CreateNodes()

    #problem data
    prop = model_part.Properties[1]
    util = ThermalAuxiliaryUtility()
    prop.SetValue(THERMAL_CONDUCTIVITY, 1.0 )
    prop.SetValue(THICKNESS, 1)
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 2)
    util.SetValueForProperties(prop, HEAT_SOURCE, HeatSourceStdProblem2(cx, cy))

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

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    # fix temperature on the boundary
    tol = 1.0e-6
    for layer, nodes_id in layer_nodes_sets.items():
        for node_id in nodes_id:
            node = model_part.Nodes[node_id]
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0.0)

    return mpatch_mp

def compute_L2_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, process_info)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
            for i in range(0, len(u)):
                ana_u = solution.GetTemperatureAt(Q[i][0], Q[i][1], Q[i][2])
                nom = nom + pow(u[i][0] - ana_u, 2) * W[i][0] * J0[i][0]
                denom = denom + pow(ana_u, 2) * W[i][0] * J0[i][0]
    print("nom:", nom)
    print("denom:", denom)
#    if denom == 0.0:
#        if nom == 0.0:
#            return 0.0
#        else:
#            return float('nan');
#    else:
#        return math.sqrt(abs(nom / denom))
    return math.sqrt(nom)

def Compute_L2_error_OnElement(element, solution, process_info):
    error = 0.0
    u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
    J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
    Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, process_info)
    W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
    for i in range(0, len(u)):
        ana_u = solution.GetTemperatureAt(Q[i][0], Q[i][1], Q[i][2])
        error = error + (pow(u[i][0] - ana_u, 2)) * W[i][0] * J0[i][0]
    return error

def main(output=True, logging=True, nsampling=2, order=3, nsteps=2):
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch, order, nsampling)
    mpatch.Enumerate()
    print(mpatch)

    if output:
        name = "square_" + str(order) + "_" + str(nsampling) + ".m"
        mpatch_export.Export(mpatch, name)

    ## create the hierarchical B-Splines multipatch
    hmpatch = CreateHBMultiPatch(mpatch)
    hpatch = hmpatch[1].GetReference()
    hpatch.FESpace().SetMaxLevel(2)

    ## extract layer
    layer_nodes_sets = ExtractLayer(hmpatch)

    ########################################
    params_post = {}
    params_post['name'] = "linear_poisson_square"
    params_post['division mode'] = "uniform"
    params_post['uniform division number'] = 40
#    params_post['division mode'] = "non-uniform"
#    params_post['division number u'] = 10
#    params_post['division number v'] = 10
#    params_post['division number w'] = 1
    params_post['variables list'] = [TEMPERATURE]
    dim = 2
    ############################

    time = 0.0
    delta_time = 1.0
    l2_error_list = []
    solution = HeatStdProblem2Solution(cx, cy)
    for i in range (0, nsteps):

        hmpatch_mp = CreateModel(hmpatch, layer_nodes_sets)
        model_part = hmpatch_mp.GetModelPart()

        #############ANALYSIS MODEL#######################################
        params_sim = model_iga_include.StaticParameters()
        params_sim['builder_and_solver_type'] = "residual-based block"
        params_sim['log_residuum'] = logging
        model = model_iga_include.Model('linear_poisson_square', os.getcwd()+"/", model_part, params_sim)
        model.solver.solver.convergence_criteria = DisplacementCriteria(model.rel_tol, model.abs_tol)
        for node in model.model_part.Nodes:
            node.AddDof(TEMPERATURE)
            node.AddDof(LAGRANGE_TEMPERATURE)
        model.InitializeModel()
        time = time + delta_time
        if output:
            hmpatch_export.Export(hmpatch, "hb_linear_poisson_square_" + str(int(time)) + ".m")

        model.Solve(time, 0, 0, 0, 0)

#        for node in model.model_part.Nodes:
#            if node.Id == 33:
#                node.SetSolutionStepValue(TEMPERATURE, 0.000438039324001)

        ################# COMPUTE ERROR ###############################
        l2_error = compute_L2_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
        print("Global L2 error:", l2_error)
        l2_error_list.append(l2_error)

        ######Synchronize back the results to multipatch
        hmpatch_mp.SynchronizeBackward(TEMPERATURE)
        ##################################################################

#         ## CHECKING
#         for node in model.model_part.Nodes:
#             print(node.GetSolutionStepValue(TEMPERATURE), hpatch.FESpace()[node.Id-1].Weight(), hpatch.FESpace()[node.Id-1].Id, node.Id)



#         msh_grid = hpatch.GridFunction(CONTROL_POINT)
#         tmp_grid = hpatch.GridFunction(TEMPERATURE)
# #        xi = [0.034715922101487, 0.534715922101487] # ok
# #        xi = [0.982642038949257, 0.482642038949257] # not ok
#         xi = [0.832502369551893, 0.417497630448107]
#         print("control point at (", xi, "):", str(msh_grid.GetValue(xi)))
#         print("temperature at (", xi, "):", str(tmp_grid.GetValue(xi)))


        ##################################################################

        ## post processing
        if output:
            model_iga_include.PostMultiPatch(hmpatch, dim, time, params_post)

        ################# HIERARCHICAL REFINE ###############################
        echo_level = IsogeometricEchoFlags.ECHO_REFINEMENT + IsogeometricEchoFlags.ECHO_REFINEMENT_DETAIL

        for hpatch_ptr in hmpatch.Patches():
            hpatch = hpatch_ptr.GetReference()

            # if i == nsteps-1:
            #     break

            ### TEST 1: REFINE EVERYTHING
#            hbsplines_refinement_util.RefineWindow(hpatch, [[0.0, 1.0], [0.0, 1.0]], echo_level)
            ###########################

            ### TEST 2: REFINE ONLY LEFT PART
            hbsplines_refinement_util.RefineWindow(hpatch, [[0.0, 0.5], [0.0, 1.0]], echo_level)
            ###########################

#             ### TEST 3: REFINE MANUALLY
# #            refine_list = [3]
#             refine_list = [3, 4, 6, 11, 16, 17, 21]
#             for fid in refine_list:
#                 print(hpatch.FESpace()[fid])
#                 hbsplines_refinement_util.Refine(hpatch, fid, echo_level)
#             hbsplines_refinement_util.LinearDependencyRefine(hpatch, 0, echo_level)
#             ###########################

#            L2ErrorList = {}
#            for elem in model_part.Elements:
#                L2Error = Compute_L2_error_OnElement(elem, solution, model_part.ProcessInfo)
#                #H1Error = model_iga_include.ComputeH1errorOnElement(elem, ana_sol, model_part.ProcessInfo)
#                L2ErrorList[elem.Id] = L2Error
#                print L2Error

#            L2ErrorList_sorted = sorted(L2ErrorList.items(), key = operator.itemgetter(1), reverse=True)
#            print L2ErrorList_sorted

#            refined_elems_count = int(1*len(model_part.Elements))
#            if refined_elems_count == 0:
#                refined_elems_count = 1

#            for j in range(0, refined_elems_count):
#                tmp = L2ErrorList_sorted[j]
#                elem_id = tmp[0]
#                elem = model_part.Elements[elem_id]
#                left = elem.GetValue(KNOT_LEFT)
#                right = elem.GetValue(KNOT_RIGHT)
#                top = elem.GetValue(KNOT_TOP)
#                bottom = elem.GetValue(KNOT_BOTTOM)
#                print("Refine window according to element", elem_id)
#                print("Left\tRight\tTop\tBottom")
#                print left, right, top, bottom
#                hbsplines_refinement_util.RefineWindow(hpatch, [[left, right], [bottom, top]], echo_level)
##                hbsplines_refinement_util.LinearDependencyRefine(hpatch, 0, echo_level)

            # re-enumerate
            hmpatch.Enumerate()

            # generate the internal cells
            for hpatch_ptr in hmpatch.Patches():
                hpatch = hpatch_ptr.GetReference()
                hpatch.FESpace().UpdateCells()

        ## extract layer
        layer_nodes_sets = ExtractLayer(hmpatch)
        print("layer_nodes_sets:", layer_nodes_sets)

        print("l2_error_list:", l2_error_list)

    return model, l2_error_list

def test():
    model, error = main(output=False, logging=False, nsteps=2)

    ### pytesting results
    assert(abs(error[0] - 0.00011232830613820067) < 1e-10)
    assert(abs(error[1] - 7.304144369299667e-05) < 1e-10)
    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=True, logging=True)
