import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']

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

def CreateModel(mpatch):
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
    for node in model_part.Nodes:
        if (abs(node.X0 - 0.0) < tol) or (abs(node.X0 - 1.0) < tol):
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0.0)

    for node in model_part.Nodes:
        if (abs(node.Y0 - 0.0) < tol) or (abs(node.Y0 - 1.0) < tol):
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0.0)

    return mpatch_mp

def compute_L2_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
#            print("u at element ", element.Id, u)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
#            print("J0 at element ", element.Id, J0)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, process_info)
#            print("Q at element ", element.Id, Q)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
#            print("W at element ", element.Id, W)
            for i in range(0, len(u)):
                ana_u = solution.GetTemperatureAt(Q[i][0], Q[i][1], Q[i][2])
                nom = nom + pow(u[i][0] - ana_u, 2) * W[i][0] * J0[i][0]
                denom = denom + pow(ana_u, 2) * W[i][0] * J0[i][0]
    print("nom:", nom)
    print("denom:", denom)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

def main(output=True, logging=True, nsampling=2, order=3):
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch, order, nsampling)
    mpatch.Enumerate()
    print(mpatch)

    if output:
        name = "square_" + str(order) + "_" + str(nsampling) + ".m"
        mpatch_export.Export(mpatch, name)

    mpatch_mp = CreateModel(mpatch)
    model_part = mpatch_mp.GetModelPart()

    #############ANALYSIS MODEL#######################################
    params_sim = model_iga_include.StaticParameters()
    params_sim['builder_and_solver_type'] = "residual-based block"
    params_sim['log_residuum'] = logging
    model = model_iga_include.Model('linear_poisson_square', os.getcwd()+"/", model_part, params_sim)
    model.solver.solver.convergence_criteria = DisplacementCriteria(model.rel_tol, model.abs_tol)
    for node in model.model_part.Nodes:
        node.AddDof(TEMPERATURE)
    model.InitializeModel()

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)

    for node in model.model_part.Nodes:
        print(node.GetSolutionStepValue(TEMPERATURE))

    ######Synchronize back the results to multipatch
    mpatch_mp.SynchronizeBackward(TEMPERATURE)
    ##################################################################

    ## compute error
    solution = HeatStdProblem2Solution(cx, cy)
    error = compute_L2_error(model.model_part.Elements, solution, model.model_part.ProcessInfo)
    print("Global L2 error:", error)

    ##################################################################

    if output:
        ## post processing
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
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    return model, error

def test():
    model, error = main(output=False, logging=False)

    ### pytesting results
    assert(abs(error - 0.01146298887382005) < 1e-10)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=True, logging=True)
