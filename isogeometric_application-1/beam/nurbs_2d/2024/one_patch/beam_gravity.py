import sys
import os

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

rho = 2000
E = 1.91e8
nu = 0.3
g = 9.81
L = 10.0
h = 0.3

def CreateMultiPatch():

    rec_ptr = geometry_factory.CreateRectangle([0.0, 0.0, 0.0], [L, h, 0.0])
    rec = rec_ptr.GetReference()
    rec.Id = 1

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(rec_ptr)

    return mpatch

def Refine(mpatch):
    print("###############REFINEMENT###############")
##    patch1_ptr = mpatch[1]
    multipatch_refine_util.DegreeElevate(mpatch[1], [1, 1])

    ins_knots_u = []
    nsampling_u = 40
    for i in range(1, nsampling_u):
        ins_knots_u.append(float(i)/nsampling_u)

    ins_knots_v = []
    nsampling_v = 5
    for i in range(1, nsampling_v):
        ins_knots_v.append(float(i)/nsampling_v)

    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots_u, ins_knots_v])


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
    gravity[1] = -g

    prop = model_part.Properties[1]
    prop.SetValue(GRAVITY, gravity )
    prop.SetValue(BODY_FORCE, body_force )
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 1)
    prop.SetValue(DENSITY,       rho )
    prop.SetValue(YOUNG_MODULUS, E )
    prop.SetValue(POISSON_RATIO, nu )
    prop.SetValue(CONSTITUTIVE_LAW, PlaneStress() )
    prop.SetValue(THICKNESS, 1.0)

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

    # fix displacement on the left
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

    return mpatch_mp

def main(output=True, logging=True):
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch)
    mpatch.Enumerate()
    print(mpatch)

    if output:
        mpatch_export.Export(mpatch, "beam.m")

    mpatch_mp = CreateModel(mpatch)
    model_part = mpatch_mp.GetModelPart()

    #############ANALYSIS MODEL#######################################
    sim_params = model_iga_include.StaticParameters()
    sim_params['log_residuum'] = logging
    model = model_iga_include.Model('beam', os.getcwd()+"/", model_part, sim_params)
    model.InitializeModel()

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
        params_post['name'] = "beam_gravity"
        params_post['division mode'] = "uniform"
        # params_post['uniform division number'] = 40
        params_post['division mode'] = "non-uniform"
        params_post['division number u'] = 100
        params_post['division number v'] = 5
        params_post['division number w'] = 1
        params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES]
        dim = 2
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)
        ##################################################################

    return model

def test():
    model = main(logging = False, output = False)

    y_disp = 0.0
    ny = 0
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 10.0) < tol:
            y_disp = y_disp + node.GetSolutionStepValue(DISPLACEMENT_Y)
            # print("deflection at control point " + str(node.Id) + ": " + str(node.GetSolutionStepValue(DISPLACEMENT_Y)))
            ny = ny + 1
    y_disp = y_disp/ny

    b = 1.0
    p = rho*g*b*h
    I = (b*h**3)/12
    ana_disp = p*L**4/(8*E*I)

    print("Computed deflection: " + str(abs(y_disp)))
    print("Analytical deflection: " + str(ana_disp))
    print("Relative error: " + str(abs(abs(y_disp) - ana_disp) / ana_disp * 100.0) + " %")

    ref_disp = -17.1099674148
    assert(abs(y_disp - ref_disp) / abs(ref_disp) < 1e-10)
    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
