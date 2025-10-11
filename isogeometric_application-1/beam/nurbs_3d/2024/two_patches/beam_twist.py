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
b = 1.0

def Refine(mpatch, nsampling_u=40, nsampling_v=5, nsampling_w=5):
    print("###############REFINEMENT###############")
##    patch1_ptr = mpatch[1]
    multipatch_refine_util.DegreeElevate(mpatch[1], [1, 1, 1])
    multipatch_refine_util.DegreeElevate(mpatch[2], [1, 0, 0])

    ins_knots_u = []
    for i in range(1, nsampling_u):
        ins_knots_u.append(float(i)/nsampling_u)

    ins_knots_v = []
    for i in range(1, nsampling_v):
        t = float(i)/nsampling_v
        ins_knots_v.append(t*t) # this will make the negative deformation gradient. We don't know yet why.
        # ins_knots_v.append(t) # this will prevent negative deformation gradient

    ins_knots_w = []
    for i in range(1, nsampling_w):
        t = float(i)/nsampling_w
        ins_knots_w.append(t*t)
        # ins_knots_w.append(t)

    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots_u, ins_knots_v, ins_knots_w])
    multipatch_refine_util.InsertKnots(mpatch[2], [ins_knots_u, [], []])

    return mpatch

def CreateModel(mpatch, nlgeom=True):
    mpatch_util = MultiPatchUtility()
    if nlgeom:
        element_name = "FiniteStrainBezier3D"
    else:
        element_name = "KinematicLinearBezier3D"

    mpatch_mp = MultiPatchModelPart3D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( model_part )
    model_part.AddNodalSolutionStepVariable(THREED_STRESSES)

    mpatch_mp.CreateNodes()

    #problem data
    body_force = ZeroVector(2)
    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0
    gravity[2] = 0.0

    prop = model_part.Properties[1]
    prop.SetValue(GRAVITY, gravity )
    prop.SetValue(BODY_FORCE, body_force )
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 1)
    prop.SetValue(DENSITY,       rho )
    prop.SetValue(YOUNG_MODULUS, E )
    prop.SetValue(POISSON_RATIO, nu )
    if nlgeom:
        prop.SetValue(CONSTITUTIVE_LAW, HyperelasticFiniteStrainBridgingConstitutiveLaw(StVenantKirchhoff_Isotropic3D()) )
    else:
        prop.SetValue(CONSTITUTIVE_LAW, Isotropic3D() )

    patch_ids = [1, 2]
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
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    return mpatch_mp

def main(output=True, logging=True, setting=1, configuration=1, nlgeom=True, nsampling_l=40, nsampling_w=5, nsampling_h=5):
    map_dir = {0: 'u', 1: 'v', 2: 'w'}

    mpatch, param_dirs = geometry_factory.CreateTwoSlabs(setting=setting, configuration=configuration, L=0.5*L, w=b, h=h, dist=0.0)

    nsampling = {}
    nsampling[map_dir[param_dirs[0]]] = nsampling_l
    nsampling[map_dir[param_dirs[1]]] = nsampling_w
    nsampling[map_dir[param_dirs[2]]] = nsampling_h
    mpatch = Refine(mpatch, nsampling_u=nsampling['u'], nsampling_v=nsampling['v'], nsampling_w=nsampling['w'])

    mpatch.Enumerate()
    if output:
        mpatch_export.Export(mpatch, "beam.m")

    print("MultiPatch validation: %r" % (mpatch.Validate()))
    multipatch_util.CheckInterfaces(mpatch)
    print(mpatch)

    mpatch_mp = CreateModel(mpatch, nlgeom=nlgeom)
    model_part = mpatch_mp.GetModelPart()

    #############ANALYSIS MODEL#######################################
    sim_params = model_iga_include.StaticParameters()
    sim_params['log_residuum'] = logging
    sim_params['max_iter'] = 30
    sim_params['calculate_reaction'] = True
    model = model_iga_include.Model('beam', os.getcwd()+"/", model_part, sim_params)
    model.InitializeModel()

    time = 0.0
    model.SolveModel(time)

    time = 1.0

    tol = 1e-6
    du = 0.3
    for node in model.model_part.Nodes:
        if abs(node.X0 - L) < tol and abs(node.Y0) < tol and abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, du)
        if abs(node.X0 - L) < tol and abs(node.Y0 - b) < tol and abs(node.Z0 - h) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, -du)

    model.SolveModel(time)

    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())

    ######Synchronize back the results to multipatch
    mpatch_mp.SynchronizeBackward(DISPLACEMENT)
    mpatch_mp.SynchronizeBackward(REACTION)
    mpatch_mp.SynchronizeBackward(THREED_STRESSES)
    ##################################################################

    if output:
        ## post processing
        params_post = {}
        params_post['name'] = "beam_twist"
        params_post['division mode'] = "non-uniform per patch"
        params_post['division number u'] = {}
        params_post['division number v'] = {}
        params_post['division number w'] = {}
        params_post['division number u'][2] = 50
        params_post['division number v'][2] = 10
        params_post['division number w'][2] = 5
        params_post['division number ' + map_dir[param_dirs[0]]][1] = 50
        params_post['division number ' + map_dir[param_dirs[1]]][1] = 10
        params_post['division number ' + map_dir[param_dirs[2]]][1] = 5
        params_post['variables list'] = [DISPLACEMENT, REACTION, THREED_STRESSES]
        dim = 3
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)
        #
        params_post['name'] = "beam_twist_control_mesh"
        params_post['mesh type'] = "control mesh"
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)
        ##################################################################

    return model

def test():
    for setting in [1, 2, 3, 4, 5, 6]:
        for configuration in [1, 2, 3, 4]:
            model = main(logging = False, output = False, setting=setting, configuration=configuration, nlgeom=True, nsampling_l=0, nsampling_w=0, nsampling_h=0)

            tol = 1e-6
            ref_reac = [2.0275984823e+05, -4.4096323031e+04, -1.0293985948e+05]
            for node in model.model_part.Nodes:
                if abs(node.X0) < tol and abs(node.Y0) < tol and abs(node.Z0) < tol:
                    reaction = node.GetSolutionStepValue(REACTION)
                    for i in range(0, 3):
                        # print("reaction: %.10e" % (reaction[i]))
                        assert(abs(reaction[i] - ref_reac[i]) / abs(ref_reac[i]) < 1e-10)

    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True, setting=1, configuration=1, nlgeom=False, nsampling_l=5, nsampling_w=2, nsampling_h=2)
        # main(logging=True, output=True, setting=1, configuration=1, nlgeom=True, nsampling_l=5, nsampling_w=2, nsampling_h=2) # negative deformation gradient; deserve for more investigation
