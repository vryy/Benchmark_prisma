import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
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

sys.path.append("../../")
import analytical_solution
# this file can also be found in ${BENCHMARK_PRISMA}/structural_application/std_problems/infinite_plate_with_hole
# where ${BENCHMARK_PRISMA} is the path of the benchmarking folder, which can be cloned from https://github.com/vryy/Benchmark_prisma.git

def CreateMultiPatch():

    #### create arc 1
    r1 = 1.0
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 0.0, 45.0)
    arc1 = arc1_ptr.GetReference()

    # create line 1
    b1 = 4.0
    line1_ptr = geometry_factory.CreateLine([b1, 0.0, 0.0], [b1, b1, 0.0], arc1.Order(0))
    line1 = line1_ptr.GetReference()

    # create patch 1
    patch1_ptr = bsplines_patch_util.CreateLoftPatch(line1, arc1)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1

    #### create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 45.0, 90.0)
    arc2 = arc2_ptr.GetReference()

    # create line 2
    line2_ptr = geometry_factory.CreateLine([b1, b1, 0.0], [0.0, b1, 0.0], arc2.Order(0))
    line2 = line2_ptr.GetReference()

    # create patch 2
    patch2_ptr = bsplines_patch_util.CreateLoftPatch(line2, arc2)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Right, patch2, BoundarySide.Left, BoundaryDirection.Forward)

    return mpatch

def Refine(mpatch):
    print("###############REFINEMENT###############")
##    patch1_ptr = mpatch[1]
    multipatch_refine_util.DegreeElevate(mpatch[1], [0, 1])

    ins_knots = []
    nsampling = 10
    for i in range(1, nsampling):
        ins_knots.append(float(i)/nsampling)

    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots, ins_knots])

##    patch2_ptr = mpatch[2]
    multipatch_refine_util.InsertKnots(mpatch[2], [ins_knots, []])

    return mpatch

def CreateModel(mpatch):
    mpatch_util = MultiPatchUtility()
    element_name = "KinematicLinearBezier2D"

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
    prop.SetValue(YOUNG_MODULUS, 1.0e5 )
    prop.SetValue(POISSON_RATIO, 0.3 )
    prop.SetValue(CONSTITUTIVE_LAW, PlaneStress() )
    prop.SetValue(THICKNESS, 1)

    patch_ids = [1, 2]
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        patch.LayerIndex = 1
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    # fix displacement on the bottom
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

    # fix displacement on the left
    for node in model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)

    b1 = 4.0
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.Y0 - b1) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)

    b1 = 4.0
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.X0 - b1) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_X)

    return mpatch_mp

def main(logging=True, output=True):
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch)
    mpatch.Enumerate()
    print(mpatch)

    if logging:
        mpatch_export.Export(mpatch, "example_plate_with_hole.m")

    mpatch_mp = CreateModel(mpatch)
    model_part = mpatch_mp.GetModelPart()

    #############ANALYSIS MODEL#######################################
    params = model_iga_include.StaticParameters()
    params['log_residuum'] = logging
    model = model_iga_include.Model('plate_with_hole', os.getcwd()+"/", model_part, params)
    model.InitializeModel()
    time = 0.0
    model.Solve(time, 0, 0, 0, 0)

    #analytical solutiion
    b1 = 4.0
    r1 = 1.0
    P = 10.0
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
    G = E/(2.0*(1.0+nu))
    kappa = (3-nu)/(1+nu) # plane stress
    ana_sol = analytical_solution.PlaneStressSolution(P, r1, G, nu)

    # prescribe displacement on the upper side
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.Y0 - b1) < tol:
            dv = ana_sol.get_displacement(node.X0, node.Y0, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_X, dv[0])
            node.SetSolutionStepValue(DISPLACEMENT_Y, dv[1])

    # prescribe displacement on the right side
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - b1) < tol:
            du = ana_sol.get_displacement(node.X0, node.Y0, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_X, du[0])
            node.SetSolutionStepValue(DISPLACEMENT_Y, du[1])

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)

    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())

    ######Synchronize back the results to multipatch
    mpatch_mp.SynchronizeBackward(DISPLACEMENT)
    mpatch_mp.SynchronizeBackward(THREED_STRESSES)
    ##################################################################

    if output:
        ## post processing
        params_post = {}
        params_post['name'] = "plate_with_hole"
        params_post['division mode'] = "uniform"
        params_post['uniform division number'] = 40
    #    params_post['division mode'] = "non-uniform"
    #    params_post['division number u'] = 10
    #    params_post['division number v'] = 10
    #    params_post['division number w'] = 1
        params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES]
        dim = 2
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    #    ## post processing
    #    r1 = 1.0
    #    b1 = 4.0
    #    disp_grid_function_1 = mpatch[1].GetReference().GridFunction(DISPLACEMENT)
    #    top_disp = disp_grid_function_1.GetValue([r1, 0.0, 0.0])
    #    print("displacement at (1.0, 0.0, 0.0): " + str(top_disp[2]))
    #    disp_grid_function_2 = mpatch[2].GetReference().GridFunction(DISPLACEMENT)
    #    left_side_disp = disp_grid_function_2.GetValue([0.0, r1, 0.0])
    #    print("displacement at (0.0, 1.0, 0.0): " + str(left_side_disp[2]))

    ## compute error
    model.l2_error = model_iga_include.ComputeL2error(model.model_part, ana_sol)
    model.h1_error = model_iga_include.ComputeH1error(model.model_part, ana_sol)

    return model

def test():
    model = main(logging=False, output=False)

    l2_error_ref = 5.9755404293107594e-04
    h1_error_ref = 1.1676942847371211e-02

    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)

    print("Test passed")

def tag():
    return "IGA"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
