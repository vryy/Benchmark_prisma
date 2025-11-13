##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.BRepApplication import *
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
mpatch_export2 = MultiNURBSPatchGLVisExporter()

import geometry_factory

sys.path.append("../../")
import analytical_solution
# this file can also be found in ${BENCHMARK_PRISMA}/structural_application/std_problems/internal_pressurized_cylinder
# where ${BENCHMARK_PRISMA} is the path of the benchmarking folder, which can be cloned from https://github.com/vryy/Benchmark_prisma.git

class Model:

    def __init__(self, E, nu, P, r1, r2, order=2, nsampling=[10, 10], plinear_solver = SuperLUSolver()):
        self.E = E
        self.nu = nu
        self.P = P
        self.r1 = r1
        self.r2 = r2
        self.order = order
        self.nsampling = nsampling
        self.plinear_solver = plinear_solver

    def CreateMultiPatch(self, order):
        ## create arc 1
        arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', self.r1, 0.0, 90.0)
        arc1 = arc1_ptr.GetReference()
        arc1.Id = 1

        ## create arc 2
        arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', self.r2, 0.0, 90.0)
        arc2 = arc2_ptr.GetReference()
        arc2.Id = 2

        ## create ring patch by connect the two arcs
        ring_patch_ptr = bsplines_patch_util.CreateLoftPatch(arc2, arc1)
        ring_patch = ring_patch_ptr.GetReference()
        ring_patch.Id = 1
        ring_patch.LayerIndex = 1

        ######create multipatch
        mpatch = MultiPatch2D()
        mpatch.AddPatch(ring_patch_ptr)

        #### elevate the degree
        if order > 1:
            multipatch_refine_util.DegreeElevate(mpatch[1], [order-2, order-1])

        return mpatch

    def Refine(self, mpatch, nsampling):
        print("###############REFINEMENT###############")
        ins_knots_u = []
        for i in range(1, nsampling[0]):
            ins_knots_u.append(float(i)/nsampling[0])

        ins_knots_v = []
        for i in range(1, nsampling[1]):
            ins_knots_v.append(float(i)/nsampling[1])

        multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots_u, ins_knots_v])

        return mpatch

    def CreateModelPart(self, mpatch):
        mpatch_util = MultiPatchUtility()
        element_name = "KinematicLinearBezier2D"
        load_condition_name = "LinePressureBezier2D"

        mpatch_mp = MultiPatchModelPart2D(mpatch)

        mpatch_mp.BeginModelPart()
        model_part = mpatch_mp.GetModelPart()
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( model_part )
        model_part.AddNodalSolutionStepVariable(THREED_STRESSES)
        model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)

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
        prop.SetValue(YOUNG_MODULUS, self.E )
        prop.SetValue(POISSON_RATIO, self.nu )
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
            load_conds = mpatch_mp.AddConditions(patch, BoundarySide2D.V1, load_condition_name, last_cond_id+1, prop)
            for cond in load_conds:
                cond.SetValue(PRESSURE, -self.P)

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

        return mpatch_mp

    def CreateModel(self, logging=True):
        mpatch = self.CreateMultiPatch(self.order)
        mpatch = self.Refine(mpatch, self.nsampling)

        mpatch.Enumerate()
        #    print(mpatch)

        if logging:
            mpatch_export.Export(mpatch, "internal_pressurized_cylinder.m")
            mpatch_export2.Export(mpatch, "internal_pressurized_cylinder.mesh")

        mpatch_mp = self.CreateModelPart(mpatch)

        #############ANALYSIS MODEL#######################################
        model_part = mpatch_mp.GetModelPart()
        sim_params = model_iga_include.StaticParameters()
        sim_params["builder_and_solver_type"] = "residual-based block"
        sim_params["linear_solver"] = self.plinear_solver
        sim_params["log_residuum"] = logging
        model = model_iga_include.Model('internal_pressurized_cylinder', os.getcwd()+"/", model_part, sim_params)
        model.InitializeModel()

        return [mpatch_mp, model]

    def Run(self, logging=True, output=True,):
        #############ANALYSIS MODEL#######################################
        [mpatch_mp, model] = self.CreateModel(logging=logging)
        mpatch = mpatch_mp.GetMultiPatch()

        time = 0.0
        model.Solve(time, 0, 0, 0, 0)

        # transfer the values from the integration points to the control point. Note that we
        # don't do that with the DISPLACEMENt because it is already stored at the control point
        transfer_util = BezierPostUtility()
        transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())
        transfer_util.TransferVariablesToNodes(YOUNG_MODULUS, model.model_part, SuperLUSolver())

        ######Synchronize back the (model_part) results to multipatch
        mpatch_mp.SynchronizeBackward(DISPLACEMENT)
        mpatch_mp.SynchronizeBackward(THREED_STRESSES)
        mpatch_mp.SynchronizeBackward(YOUNG_MODULUS)
        ##################################################################

        ## compute error
        ana_sol = analytical_solution.Solution(self.r1, self.r2, self.P, 0.0, self.E, self.nu)
        print("********Error Analysis at time " + str(time) + "*******")
        l2_error = model_iga_include.ComputeL2error(model.model_part, ana_sol)
        h1_error = model_iga_include.ComputeH1error(model.model_part, ana_sol)
        print("*******************************************")

        if output:

            ## post processing
            params_post = {}
            params_post['name'] = "internal_pressurized_cylinder"
            params_post['division mode'] = "uniform"
            params_post['uniform division number'] = 20
            # params_post['output format'] = "ascii"
        #    params_post['division mode'] = "non-uniform"
        #    params_post['division number u'] = 10
        #    params_post['division number v'] = 10
        #    params_post['division number w'] = 1
            params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES, YOUNG_MODULUS]
            params_post['backend'] = ["GiD", "Glvis", "ParaView"]
            params_post['paraview_element_type'] = "Bezier"
            dim = 2
            model_iga_include.PostMultiPatch(mpatch, dim, time, params_post, model_part=model.model_part)

        #    ## post processing
        #    r1 = 1.0
        #    b1 = 4.0
        #    disp_grid_function_1 = mpatch[1].GetReference().GridFunction(DISPLACEMENT)
        #    top_disp = disp_grid_function_1.GetValue([r1, 0.0, 0.0])
        #    print("displacement at (1.0, 0.0, 0.0): " + str(top_disp[2]))
        #    disp_grid_function_2 = mpatch[2].GetReference().GridFunction(DISPLACEMENT)
        #    left_side_disp = disp_grid_function_2.GetValue([0.0, r1, 0.0])
        #    print("displacement at (0.0, 1.0, 0.0): " + str(left_side_disp[2]))

        return l2_error, h1_error
