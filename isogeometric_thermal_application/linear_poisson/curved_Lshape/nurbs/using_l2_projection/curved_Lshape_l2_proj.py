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

kernel = Kernel()   #defining kernel

import model_iga_include
from model_iga_include import *

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

import geometry_factory

mpatch_import = MultiNURBSPatchGeoImporter2D()

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

def CreateModel(mpatch):
    mpatch_util = MultiPatchUtility()
    element_name = "LinearPoissonBezier2D"
    condition_name = "DummyConditionBezier1D2"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)

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

    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)

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

def main(logging=True, output=True, nsampling=2):
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch, nsampling=nsampling)
#    mpatch.Enumerate()
#    print(mpatch)

#    for patch_ptr in mpatch.Patches():
#        patch = patch_ptr.GetReference()
#        print(patch)
#    sys.exit(0)

    if output:
        mpatch_export.Export(mpatch, "curved_Lshape.m")

    mpatch_mp = CreateModel(mpatch)
    model_part = mpatch_mp.GetModelPart()

    #############ANALYSIS MODEL#######################################
    sim_params = model_iga_include.StaticParameters()
    sim_params['log_residuum'] = logging
    model = model_iga_include.Model('linear_poisson_curved_Lshape', os.getcwd()+"/", model_part, sim_params)
    model.solver.solver.convergence_criteria = DisplacementCriteria(model.rel_tol, model.abs_tol)
    for node in model.model_part.Nodes:
        node.AddDof(TEMPERATURE)
    model.InitializeModel()

    # project the Dirichlet boundary values
    dirichlet_proc = DirichletL2ProjectionProcess(model_part, model_part.Conditions, MKLPardisoSolver())
    dirichlet_proc.SetTemperatureFunction(HeatStdProblem1Solution())
    dirichlet_proc.Execute()

    # solve the system
    time = 1.0
    model.Solve(time, 0, 0, 0, 0)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(TEMPERATURE))

    solution = HeatStdProblem1Solution()
    l2_error = compute_L2_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
    h1_error = compute_H1_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
    print("Global L2 error: %.16e" % l2_error)
    print("Global H1 error: %.16e" % h1_error)
    model.l2_error = l2_error
    model.h1_error = h1_error

    ######Synchronize back the results to multipatch
    mpatch_mp.SynchronizeBackward(TEMPERATURE)
    ##################################################################

    if output:
        ## post processing
        params_post = {}
        params_post['name'] = "linear_poisson_curved_Lshape"
        params_post['division mode'] = "uniform"
        params_post['uniform division number'] = 40
        params_post['variables list'] = [TEMPERATURE]
        dim = 2
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    return model

def test():

    model = main(output=False, logging=False, nsampling=2)

    print("%.16e" % model.h1_error)
    l2_error_ref = 2.8552495988074197e-03
    h1_error_ref = 4.3472717999962079e-02
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    print("Test passed")

def study():

    ifile = open("convergence.log", "w")
    ifile.write("%-*s%-*s%s\n" % (10, "ndofs", 30, "l2_error", "h1_error"))

    nsampling = 2
    for i in range(0, 5):
        model = main(output=False, logging=False, nsampling=nsampling)

        ndofs = model.solver.solver.builder_and_solver.GetEquationSystemSize()

        ifile.write("%-*d%-*.16e%.16e\n" % (10, ndofs, 30, model.l2_error, model.h1_error))

        nsampling *= 2

    ifile.close()

def tag():
    return "nurbs"

def print_tag():
    print("Tags: " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(logging=True, output=True)
