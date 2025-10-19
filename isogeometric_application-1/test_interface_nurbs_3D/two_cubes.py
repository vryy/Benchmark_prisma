import sys
import os

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# from KratosMultiphysics.MKLSolversApplication import *

kernel = Kernel()   #defining kernel

import geometry_factory
import model_iga_include

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

def Refine(mpatch):
    print("###############REFINEMENT###############")
##    patch1_ptr = mpatch[1]

    ins_knots = []
    nsampling = 10
    for i in range(1, nsampling):
        ins_knots.append((float(i)/nsampling)**2)

    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots, ins_knots])

    return mpatch

def CreateModel(mpatch):
    mpatch_util = MultiPatchUtility()
    element_name = "KinematicLinearBezier3D"

    mpatch_mp = MultiPatchModelPart3D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( model_part )
    model_part.AddNodalSolutionStepVariable(THREED_STRESSES)

    mpatch_mp.CreateNodes()

    #problem data
    body_force = ZeroVector(3)
    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0
    gravity[2] = 0.0
    prop = model_part.Properties[1]
    prop.SetValue(GRAVITY, gravity )
    prop.SetValue(BODY_FORCE, body_force )
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 1)
    prop.SetValue(DENSITY,            0 )
    prop.SetValue(YOUNG_MODULUS, 1.0e5 )
    prop.SetValue(POISSON_RATIO, 0.3 )
    prop.SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
    prop.SetValue(THICKNESS, 1.0)

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

    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)

        if abs(node.X0 - 2.1) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)

        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    return mpatch_mp

def main(logging=True, output=True, setting=1, configuration=1, dist=0.1):
    mpatch, dummy = geometry_factory.CreateTwoSlabs(setting=setting, configuration = configuration, dist = dist, L=1.0, w=1.0, h=1.0)
    # mpatch = Refine(mpatch)
    mpatch.Enumerate()
    print(mpatch)

    if logging:
        mpatch_export.Export(mpatch, "two_cubes.m")

    multipatch_util.CheckInterfaces(mpatch, True)

    mpatch_mp = CreateModel(mpatch)
    model_part = mpatch_mp.GetModelPart()

    #############ANALYSIS MODEL#######################################
    params = model_iga_include.StaticParameters()
    params['log_residuum'] = logging
    model = model_iga_include.Model('two_cubes', os.getcwd()+"/", model_part, params)
    model.InitializeModel()
    time = 0.0
    model.SolveModel(time)

    # prescribe displacement on the right side
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - (2.0+dist)) < tol:
            node.SetSolutionStepValue(DISPLACEMENT_X, 1.0)

    time = 1.0
    model.SolveModel(time)

    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())

    ######Synchronize back the results to multipatch
    mpatch_mp.SynchronizeBackward(DISPLACEMENT)
    mpatch_mp.SynchronizeBackward(THREED_STRESSES)
    ##################################################################

    if output:
        ## post processing
        params_post = {}
        params_post['name'] = "two_cubes"
        params_post['division mode'] = "uniform"
        params_post['uniform division number'] = 10
    #    params_post['division mode'] = "non-uniform"
    #    params_post['division number u'] = 10
    #    params_post['division number v'] = 10
    #    params_post['division number w'] = 1
        params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES]
        dim = 3
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

        params_post['name'] = "two_cubes_control_mesh"
        params_post['mesh type'] = "control mesh"
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    return model, mpatch

def is_same(l1, l2):
    for i in range(0, len(l1)):
        if l1[i] != l2[i]:
            return False
    return True

def test():
    ref_interface_indices = {}
    ref_interface_indices[1] = {}
    ref_interface_indices[1][1] = [1, 3, 5, 7]
    ref_interface_indices[1][2] = [4, 6, 0, 2]
    ref_interface_indices[1][3] = [2, 0, 6, 4]
    ref_interface_indices[1][4] = [7, 5, 3, 1]
    ref_interface_indices[2] = {}
    ref_interface_indices[2][1] = [0, 4, 2, 6]
    ref_interface_indices[2][2] = [5, 1, 7, 3]
    ref_interface_indices[2][3] = [3, 7, 1, 5]
    ref_interface_indices[2][4] = [6, 2, 4, 0]
    ref_interface_indices[3] = {}
    ref_interface_indices[3][1] = [0, 1, 4, 5]
    ref_interface_indices[3][2] = [6, 7, 2, 3]
    ref_interface_indices[3][3] = [3, 2, 7, 6]
    ref_interface_indices[3][4] = [5, 4, 1, 0]
    ref_interface_indices[4] = {}
    ref_interface_indices[4][1] = [2, 6, 3, 7]
    ref_interface_indices[4][2] = [4, 0, 5, 1]
    ref_interface_indices[4][3] = [1, 5, 0, 4]
    ref_interface_indices[4][4] = [7, 3, 6, 2]
    ref_interface_indices[5] = {}
    ref_interface_indices[5][1] = [4, 5, 6, 7]
    ref_interface_indices[5][2] = [2, 3, 0, 1]
    ref_interface_indices[5][3] = [1, 0, 3, 2]
    ref_interface_indices[5][4] = [7, 6, 5, 4]
    ref_interface_indices[6] = {}
    ref_interface_indices[6][1] = [0, 2, 1, 3]
    ref_interface_indices[6][2] = [6, 4, 7, 5]
    ref_interface_indices[6][3] = [5, 7, 4, 6]
    ref_interface_indices[6][4] = [3, 1, 2, 0]

    settings = [1, 2, 3, 4, 5, 6]
    configurations = [1, 2, 3, 4]
    # settings = [1]
    # configurations = [2]
    for setting in settings:
        for configuration in configurations:
            model, mpatch = main(logging=False, output=False, setting=setting, configuration=configuration, dist=0.0)
            mpatch.Enumerate()
            multipatch_util.CheckInterfaces(mpatch, True) # check the indices on patch interface

            bpatch_ptr = mpatch[2].GetReference().ConstructBoundaryPatch(BoundarySide3D.U0)
            func_indices = bpatch_ptr.GetReference().FESpace().FunctionIndices()
            print(func_indices)
            print(setting)
            print(configuration)
            assert(is_same(func_indices, ref_interface_indices[setting][configuration]))

    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True, setting=6, configuration=4, dist=0.0)
