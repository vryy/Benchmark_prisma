import sys
import os
import operator
##################################################################
import pdb
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.IsogeometricThermalApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MultithreadedSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.VisualApplication import *
##################################################################

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

mpatch_import = MultiNURBSPatchGeoImporter2D()
mpatch_util = MultiPatchUtility()

def CreateMultiPatch():

    ######import the multipatch
    mpatch = mpatch_import.Import("../../geo_curvedL_3patches.txt")
    # original source: GeoPDEs-3.1.0/geopdes/inst/examples/geometry_files/multipatch/geo_curvedL_3patches.txt

    return mpatch

def Refine(mpatch, nsampling=2):
    print("###############REFINEMENT###############")
    multipatch_refine_util.DegreeElevate(mpatch[1], [1, 2])
    multipatch_refine_util.DegreeElevate(mpatch[3], [1, 2])

    ins_knots = []
    for i in range(1, nsampling):
        ins_knots.append(float(i)/nsampling)

    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots, ins_knots])
    multipatch_refine_util.InsertKnots(mpatch[2], [ins_knots, []])
    multipatch_refine_util.InsertKnots(mpatch[3], [[], ins_knots])

    return mpatch

def CreateHBMultiPatch(mpatch):
    # convert the NURBS patches to HB-Splines patches
    hmpatch = hbsplines_patch_util.CreateMultiPatchFromBSplines(mpatch)

    # print(hmpatch)
    return hmpatch

def ExtractLayer(hmpatch):
    layer_nodes_sets = {}
    layer_nodes_sets['left'] = []
    layer_nodes_sets['bottom'] = []
    layer_nodes_sets['inner'] = []
    layer_nodes_sets['outer'] = []

    for hpatch_ptr in hmpatch.Patches():
        hpatch = hpatch_ptr.GetReference()

        boundary_basis_1_u0 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.U0))
        for i in range(0, len(boundary_basis_1_u0)):
            bf = boundary_basis_1_u0[i]
            layer_nodes_sets['bottom'].append(bf.EquationId + 1)
        boundary_basis_1_u1 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.U1))
        for i in range(0, len(boundary_basis_1_u1)):
            bf = boundary_basis_1_u1[i]
            layer_nodes_sets['left'].append(bf.EquationId + 1)
        boundary_basis_1_v0 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.V0))
        for i in range(0, len(boundary_basis_1_v0)):
            bf = boundary_basis_1_v0[i]
            layer_nodes_sets['outer'].append(bf.EquationId + 1)
        boundary_basis_1_v1 = hpatch.FESpace().GetBoundaryBfs(multipatch_util.BoundaryFlag(BoundarySide2D.V1))
        for i in range(0, len(boundary_basis_1_v1)):
            bf = boundary_basis_1_v1[i]
            layer_nodes_sets['inner'].append(bf.EquationId + 1)

    return layer_nodes_sets

def CreateModel(mpatch, layer_nodes_sets):
    element_name = "LinearPoissonBezier2D"
    condition_name = "DummyConditionBezier1D2"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(ERROR_RATIO)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)

    mpatch_mp.CreateNodes()

    #problem data
    prop = model_part.Properties[1]
    prop.SetValue(THERMAL_CONDUCTIVITY, 1.0 )
    prop.SetValue(THICKNESS, 1)
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 2)
    util = ThermalAuxiliaryUtility()
    util.SetValueForProperties(prop, HEAT_SOURCE, HeatSourceStdProblem1())
    util.SetValueForProperties(prop, TEMPERATURE_FUNCTION, HeatStdProblem1Solution())

    # print(120)
    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        # print(patch)
        # sys.exit(0)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)
    # print(128)

    ## add constraint conditions
    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[1], BoundarySide2D.V1, condition_name, last_cond_id+1, prop)

    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[1], BoundarySide2D.V0, condition_name, last_cond_id+1, prop)
    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[1], BoundarySide2D.U0, condition_name, last_cond_id+1, prop)

    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[2], BoundarySide2D.V1, condition_name, last_cond_id+1, prop)
    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[2], BoundarySide2D.U1, condition_name, last_cond_id+1, prop)

    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[3], BoundarySide2D.V0, condition_name, last_cond_id+1, prop)
    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[3], BoundarySide2D.U0, condition_name, last_cond_id+1, prop)
    last_cond_id = mpatch_util.GetLastConditionId(model_part)
    load_conds = mpatch_mp.AddConditions(mpatch[3], BoundarySide2D.U1, condition_name, last_cond_id+1, prop)

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    return mpatch_mp

def compute_L2_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
            for i in range(0, len(u)):
                ana_u = solution.GetTemperatureAt(Q[i][0], Q[i][1], Q[i][2])
                nom = nom + pow(u[i][0] - ana_u, 2) * W[i][0] * J0[i][0]
                denom = denom + pow(ana_u, 2) * W[i][0] * J0[i][0]
    # print("nom:", nom)
    # print("denom:", denom)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(nom / denom)

def compute_H1_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            du = element.GetValuesOnIntegrationPoints(TEMPERATURE_GRADIENT, process_info)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
            for i in range(0, len(du)):
                ana_du = solution.GetTemperatureGradientAt(Q[i][0], Q[i][1], Q[i][2])
                dist = 0.0
                norm_ana_du = 0.0
                for j in range(0, 3):
                    dist += pow(du[i][j] - ana_du[j], 2)
                    norm_ana_du += pow(ana_du[j], 2)
                nom = nom + dist * W[i][0] * J0[i][0]
                denom = denom + norm_ana_du * W[i][0] * J0[i][0]
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

def Compute_L2_error_OnElement(element, solution, process_info):
    error = 0.0
    u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
    J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
    Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
    W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
    for i in range(0, len(u)):
        ana_u = solution.GetTemperatureAt(Q[i][0], Q[i][1], Q[i][2])
        error = error + (pow(u[i][0] - ana_u, 2)) * W[i][0] * J0[i][0]
    return error

def main(logging=True, output=True, nsteps=3):
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch)
    mpatch.Enumerate()
    mpatch_util.CheckInterfaces(mpatch)

    ## create the hierarchical B-Splines multipatch
    hmpatch = CreateHBMultiPatch(mpatch)
    mpatch_util.CheckInterfaces(hmpatch)

    ## extract layer
    layer_nodes_sets = ExtractLayer(hmpatch)
    print("layer_nodes_sets:", layer_nodes_sets)

    ########################################
    params_post = {}
    params_post['name'] = "linear_poisson_curved_Lshape"
    params_post['division mode'] = "uniform"
    params_post['uniform division number'] = 5 # 40
    params_post['variables list'] = [TEMPERATURE, DISPLACEMENT]
    dim = 2
    ############################

    time = 0.0
    delta_time = 1.0
    l2_error_list = []
    h1_error_list = []
    ndofs_list = []
    for i in range (0, nsteps):

        hmpatch_mp = CreateModel(hmpatch, layer_nodes_sets)
        model_part = hmpatch_mp.GetModelPart()

        ## set the required information for refinement
        system_size = hmpatch.EquationSystemSize()
        for idof in range(0, system_size):
            bf = hbsplines_patch_util.GetBfByEquationId(hmpatch, idof)
            node = model_part.Nodes[idof+1]
            node.SetValue(HIERARCHICAL_LEVEL, bf.Level())
            node.SetValue(BASIS_FUNCTION_INDEX, bf.Id)
            node.SetValue(EQUATION_INDEX, idof)

        #############ANALYSIS MODEL#######################################
        sim_params = model_iga_include.StaticParameters()
        sim_params['linear_solver'] = SuperLUMTSolver()
        #sim_params['linear_solver'] = SuperLUSolver()
        model = model_iga_include.Model('linear_poisson_curved_Lshape', os.getcwd()+"/", model_part, sim_params)
        model.solver.solver.convergence_criteria = DisplacementCriteria(model.rel_tol, model.abs_tol)
        for node in model.model_part.Nodes:
            node.AddDof(TEMPERATURE)
        model.InitializeModel()

        time = time + delta_time
        # hmpatch_export.Export(hmpatch, "hb_curved_Lshaped_" + str(int(time)) + ".m")

        # project the Dirichlet boundary values
        dirichlet_proc = DirichletL2ProjectionProcess(model_part, model_part.Conditions, SuperLUSolver())
        dirichlet_proc.SetTemperatureFunction(HeatStdProblem1Solution())
        dirichlet_proc.Execute()

        # solve the system
        model.Solve(time, 0, 0, 0, 0)

        solution = HeatStdProblem1Solution()

        ################# COMPUTE ERROR ###############################
        l2_error = compute_L2_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
        h1_error = compute_H1_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
        print("Global L2 error:", l2_error)
        print("Global H1 error:", h1_error)
        l2_error_list.append(l2_error)
        h1_error_list.append(h1_error)
        ndofs_list.append(model.solver.solver.builder_and_solver.GetEquationSystemSize())

        # if i == 0:
        for node in model.model_part.Nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Z, node.GetSolutionStepValue(TEMPERATURE))

        ######Synchronize back the results to multipatch
        hmpatch_mp.SynchronizeBackward(TEMPERATURE)
        hmpatch_mp.SynchronizeBackward(DISPLACEMENT)
        ##################################################################

        ##################################################################
        ## ERROR ESTIMATOR

        C0_coeff = 1.0
        error_estimator = ComputeErrorEstimateLaplaceProcessIsogeometric(model.model_part.Elements, model.model_part.ProcessInfo, C0_coeff)
        error_estimator.Execute()

        ##################################################################

        if output:
            ## post processing
            model_iga_include.PostMultiPatch(hmpatch, dim, time, params_post)

        ################# HIERARCHICAL REFINE ###############################
        # echo_level = IsogeometricEchoFlags.ECHO_REFINEMENT + IsogeometricEchoFlags.ECHO_REFINEMENT_DETAIL
        echo_level = 0

        if i == nsteps-1:
            break

        # ### TEST: REFINE EVERYTHING
        # for hpatch_ptr in hmpatch.Patches():
        #     hpatch = hpatch_ptr.GetReference()
        #     hbsplines_refinement_util.RefineWindow(hpatch, [[0.0, 1.0], [0.0, 1.0]], echo_level)
        # ###########################

        # ### TEST: REFINE A CORNER
        # for hpatch_ptr in hmpatch.Patches():
        #     hpatch = hpatch_ptr.GetReference()
        #     if hpatch.Id == 1:
        #         hbsplines_refinement_util.RefineWindow(hpatch, [[0.5, 1.0], [0.0, 0.5]], echo_level)

        # # re-enumerate
        # hmpatch.Enumerate()

        # # generate the internal cells
        # for hpatch_ptr in hmpatch.Patches():
        #     hpatch = hpatch_ptr.GetReference()
        #     hpatch.FESpace().UpdateCells()
        # ###########################

        ### TEST 4: Using maximum marking strategy
        mark_param = 0.5
        refine_list = []
        max_est = -1.0e99
        for node in model.model_part.Nodes:
            err = node.GetSolutionStepValue(ERROR_RATIO)
            if (err > max_est):
                max_est = err

        for node in model.model_part.Nodes:
            err = node.GetSolutionStepValue(ERROR_RATIO)
            if (err > mark_param*max_est):
                refine_list.append(node.GetValue(EQUATION_INDEX))

        # print("refine_list:", refine_list)
        refine_bf_list = []
        aux_list = []
        for eid in refine_list:
            for hpatch_ptr in hmpatch.Patches():
                hpatch = hpatch_ptr.GetReference()
                if hpatch.FESpace().HasBfByEquationId(eid):
#                    print(hpatch.FESpace()[fid])
                    bf = hpatch.FESpace().GetBfByEquationId(eid)
                    refine_bf_list.append([bf, hpatch])
                    break

        # print("refine_bf_list:", refine_bf_list)
        # print("len(refine_bf_list):", len(refine_bf_list))
        for tmp in refine_bf_list:
            bf = tmp[0]
            hpatch = tmp[1]
            hbsplines_refinement_util.Refine(hpatch, bf, echo_level)
            hbsplines_refinement_util.LinearDependencyRefine(hpatch, 0, echo_level)

            # re-enumerate
            hmpatch.Enumerate()

            # generate the internal cells
            for hpatch_ptr in hmpatch.Patches():
                hpatch = hpatch_ptr.GetReference()
                print("UpdateCells for patch %d begin" % (hpatch.Id), flush=True)
                hpatch.FESpace().UpdateCells()

            mpatch_util.CheckInterfaces(hmpatch)
            hbsplines_patch_util.ReportDuplicatedEquationId(hmpatch, True)

        ###########################

        ## extract layer
        layer_nodes_sets = ExtractLayer(hmpatch)
        # print("layer_nodes_sets:", layer_nodes_sets)

    print(f"l2_error_list: {l2_error_list}")
    print(f"h1_error_list: {h1_error_list}")
    print(f"ndofs_list: {ndofs_list}")

    ### COMPATIBILITY TEST 2
    lc = 0.1 # local reference point

    bpatch1 = hmpatch[1].GetReference().ConstructBoundaryPatch(BoundarySide2D.U1)
#        print("boundary fespace 1:")
#        print(bpatch1.GetReference().FESpace())
    cf1 = bpatch1.GetReference().GridFunction(CONTROL_POINT)
    tf1 = bpatch1.GetReference().GridFunction(TEMPERATURE)
    tf1big = hmpatch[1].GetReference().GridFunction(TEMPERATURE)
    print("boundary point 1: ", str(cf1.GetValue([lc])))
    print("boundary point 1 temperature: ", str(tf1.GetValue([lc])))
    print("boundary point 1 temperature by bulk: ", str(tf1big.GetValue([1.0, lc])))

    bpatch2 = hmpatch[2].GetReference().ConstructBoundaryPatch(BoundarySide2D.U0)
#        print("boundary fespace 2:")
#        print(bpatch2.GetReference().FESpace())
    cf2 = bpatch2.GetReference().GridFunction(CONTROL_POINT)
    tf2 = bpatch2.GetReference().GridFunction(TEMPERATURE)
    tf2big = hmpatch[2].GetReference().GridFunction(TEMPERATURE)
    print("boundary point 2: ", str(cf2.GetValue([lc])))
    print("boundary point 2 temperature: ", str(tf2.GetValue([lc])))
    print("boundary point 2 temperature by bulk: ", str(tf2big.GetValue([0.0, lc])))
    ########################

    model.ndofs_list = ndofs_list
    model.l2_error_list = l2_error_list
    model.h1_error_list = h1_error_list

    if logging:
        ifile = open("convergence.log", "w")
        ifile.write("%-*s%-*s%s\n" % (10, "ndofs", 20, "l2_error", "h1_error"))
        for i in range(0, len(ndofs_list)):
            ifile.write("%-*d%-*.10e%.10e\n" % (10, ndofs_list[i], 20, l2_error_list[i], h1_error_list[i]))
        ifile.close()

    # return the last model
    return model

def test():
    model = main(logging=False, output=False, nsteps=6)

    print("%.16e" % model.l2_error_list[-1])
    print("%.16e" % model.h1_error_list[-1])
    ref_l2_error = 4.9460358851197615e-05
    ref_h1_error = 4.9757316339100102e-03
    assert(abs(model.l2_error_list[-1] - ref_l2_error) < 1e-10)
    assert(abs(model.h1_error_list[-1] - ref_h1_error) < 1e-10)

    print("Test passed")

def tag():
    return "hbsplines"

def print_tag():
    print("Tags: " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(logging=True, output=True, nsteps=16)
