import sys
import os
import math
kratos_root_path=os.environ['KRATOS_ROOT_PATH']

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.PlateAndShellApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.IsogeometricPlateAndShellApplication import *
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

import geometry_factory

class AnalyticalSolution:

    def __init__(self, q, G, nu, h, R):
        self.q = q
        self.G = G
        self.nu = nu
        self.h = h
        self.R = R

    def get_displacement(self, x, y):
        R2 = (self.R)**2
        R4 = R2**2
        h3 = (self.h)**3
        r = math.sqrt(x**2 + y**2)
        r2 = r**2
        wc = 0.3*self.q/(self.G*self.h)*(R2 - r2) + 3*self.q*(1-self.nu)/(32*self.G*h3)*(R2-r2)**2
        return wc

    def get_psi_r(self, x, y):
        R2 = (self.R)**2
        r = math.sqrt(x**2 + y**2)
        r2 = r**2
        psi = 3.0/8*(1.0-self.nu)*self.q/(self.G*(self.h)**3)*r*(R2 - r2)
        return psi

    def get_varphi_r(self, x, y):
        r = math.sqrt(x**2 + y**2)
        varphi = -0.6*self.q*r/(self.G*self.h)
        return varphi

    def get_kh_displacement(self, x, y):
        R2 = (self.R)**2
        R4 = R2**2
        h3 = (self.h)**3
        r = math.sqrt(x**2 + y**2)
        r2 = r**2
        wc = 3*self.q*(1-self.nu)/(32*self.G*h3)*(R2-r2)**2
        return wc

    def ExportUr(self, filename, npoints=100):
        ifile = open(filename, "w")
        ifile.write("r\tw\tw_kh\tpsi_r\tvarphir\n")
        for i in range(0, npoints):
            xi = float(i)/(npoints-1)
            r = xi*self.R
            ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (r, self.get_displacement(r, 0.0), self.get_kh_displacement(r, 0.0), self.get_psi_r(r, 0.0), self.get_varphi_r(r, 0.0)))
        ifile.close()

class Model:

    def __init__(self, E, nu, r, h, f, order=2, nsampling=10, plinear_solver = SuperLUSolver()):
        self.young_modulus = E
        self.poisson_ratio = nu
        self.radius = r
        self.thickness = h
        self.face_load = f
        self.order = order
        self.nsampling = nsampling
        self.plinear_solver = plinear_solver

    def CreateMultiPatch(self, radius, order=2):
        ####### create arc 1
        arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', radius, -45.0, 45.0)
        arc1 = arc1_ptr.GetReference()

        # create line 1
        b = radius/5
        line1_ptr = geometry_factory.CreateLine([b, -b, 0.0], [b, b, 0.0], arc1.Order(0))
        line1 = line1_ptr.GetReference()

        # create patch 1
        patch1_ptr = bsplines_patch_util.CreateLoftPatch(arc1, line1)
        patch1 = patch1_ptr.GetReference()
        patch1.Id = 1

        ####### create arc 2
        arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', radius, 45.0, 135.0)
        arc2 = arc2_ptr.GetReference()

        # create line 2
        line2_ptr = geometry_factory.CreateLine([b, b, 0.0], [-b, b, 0.0], arc2.Order(0))
        line2 = line2_ptr.GetReference()

        # create patch 2
        patch2_ptr = bsplines_patch_util.CreateLoftPatch(arc2, line2)
        patch2 = patch2_ptr.GetReference()
        patch2.Id = 2

        ####### create arc 3
        arc3_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', radius, 135.0, 225.0)
        arc3 = arc3_ptr.GetReference()

        # create line 3
        line3_ptr = geometry_factory.CreateLine([-b, b, 0.0], [-b, -b, 0.0], arc3.Order(0))
        line3 = line3_ptr.GetReference()

        # create patch 3
        patch3_ptr = bsplines_patch_util.CreateLoftPatch(arc3, line3)
        patch3 = patch3_ptr.GetReference()
        patch3.Id = 3

        ####### create arc 4
        arc4_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', radius, 225.0, 315.0)
        arc4 = arc4_ptr.GetReference()

        # create line 4
        line4_ptr = geometry_factory.CreateLine([-b, -b, 0.0], [b, -b, 0.0], arc4.Order(0))
        line4 = line4_ptr.GetReference()

        # create patch 4
        patch4_ptr = bsplines_patch_util.CreateLoftPatch(arc4, line4)
        patch4 = patch4_ptr.GetReference()
        patch4.Id = 4

        ####### create line 5
        line5_ptr = geometry_factory.CreateLine([-b, -b, 0.0], [-b, b, 0.0], arc3.Order(0))
        line5 = line5_ptr.GetReference()

        # create patch 5
        patch5_ptr = bsplines_patch_util.CreateLoftPatch(line1, line5)
        multipatch_refine_util.DegreeElevate(patch5_ptr, [0, 1])
        patch5 = patch5_ptr.GetReference()
        patch5.Id = 5

        # # print(patch2)
        # print(patch3)
        # print("line2:")
        # print(line2)

        ######create multipatch
        mpatch = MultiPatch2D()
        mpatch.AddPatch(patch1_ptr)
        mpatch.AddPatch(patch2_ptr)
        mpatch.AddPatch(patch3_ptr)
        mpatch.AddPatch(patch4_ptr)
        mpatch.AddPatch(patch5_ptr)
        bsplines_patch_util.MakeInterface(patch1, BoundarySide.Right, patch2, BoundarySide.Left, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch2, BoundarySide.Right, patch3, BoundarySide.Left, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch3, BoundarySide.Right, patch4, BoundarySide.Left, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch4, BoundarySide.Right, patch1, BoundarySide.Left, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch1, BoundarySide.Top, patch5, BoundarySide.Bottom, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch2, BoundarySide.Top, patch5, BoundarySide.Right, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch3, BoundarySide.Top, patch5, BoundarySide.Top, BoundaryDirection.Reversed)
        bsplines_patch_util.MakeInterface(patch4, BoundarySide.Top, patch5, BoundarySide.Left, BoundaryDirection.Reversed)

        if order >= 2:
            multipatch_refine_util.DegreeElevate(mpatch[1], [order-2, order-1])

        return mpatch

    def Refine(self, mpatch, nsampling):
        print("###############REFINEMENT###############")

        ins_knots = []
        for i in range(1, nsampling):
            ins_knots.append(float(i)/nsampling)

        multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots, ins_knots])
        multipatch_refine_util.InsertKnots(mpatch[2], [ins_knots, []])
        # multipatch_refine_util.InsertKnots(mpatch[3], [[], ins_knots])

        return mpatch

    def CreateModel(self, mpatch):
        mpatch_util = MultiPatchUtility()
        element_name = "ReissnerMindlinElasticPlateElementBezier2D"
        condition_name = "LineLoadBezier2D"

        mpatch_mp = MultiPatchModelPart2D(mpatch)

        mpatch_mp.BeginModelPart()
        model_part = mpatch_mp.GetModelPart()
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( model_part )
        model_part.AddNodalSolutionStepVariable(THREED_STRESSES)
        model_part.AddNodalSolutionStepVariable(TRUE_SHEAR_ROTATION)

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
        prop.SetValue(YOUNG_MODULUS,        self.young_modulus )
        prop.SetValue(POISSON_RATIO,          self.poisson_ratio )
        prop.SetValue(THICKNESS,            self.thickness )

        patch_elems = {}
        for patch_ptr in mpatch.Patches():
            patch = patch_ptr.GetReference()
            patch.LayerIndex = 1
            #        print(patch)

            ## add volume elements
            last_elem_id = mpatch_util.GetLastElementId(model_part)
            elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)
            # print(len(elems))
            patch_elems[patch.Id] = elems

        ## add conditions
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        mpatch_mp.AddConditions(mpatch[1], BoundarySide2D.V0, condition_name, last_cond_id+1, prop)
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        mpatch_mp.AddConditions(mpatch[2], BoundarySide2D.V0, condition_name, last_cond_id+1, prop)
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        mpatch_mp.AddConditions(mpatch[3], BoundarySide2D.V0, condition_name, last_cond_id+1, prop)
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        mpatch_mp.AddConditions(mpatch[4], BoundarySide2D.V0, condition_name, last_cond_id+1, prop)

        mpatch_mp.EndModelPart()
        #    print(mpatch_mp)

        return [patch_elems, mpatch_mp]

    def ExportUr(self, filename, mpatch, radius, npoints=100):
        ifile = open(filename, "w")
        ifile.write("r\tw\tpsi\tvarphi_pp\tvarphi\n")
        for i in range(0, npoints):
            xi = float(i)/(npoints-1)
            r = xi*radius
            # [patch_id, xi] = mpatch.LocalCoordinates([r, 0.0, 0.0], [xi, 0.0, 0.0])
            [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [r, 0.0, 0.0], [10, 10])
            patch = mpatch[patch_id].GetReference()
            cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
            point = cgf.GetValue(xi)
            # print("patch_id: " + str(patch_id) + ", xi: " + str(xi) + ", point: " + str(point))
            if patch_id < 0:
                print("patch_id: %d" % patch_id)
                print("r: %e" % r)
                raise Exception("Fail computing the local coordinates")
            dgf = patch.GridFunction(DISPLACEMENT)
            disp = dgf.GetValue(xi)
            pgf = patch.GridFunction(ROTATION)
            psi = pgf.GetValue(xi)
            vgf = patch.GridFunction(TRUE_SHEAR_ROTATION)
            varphi = vgf.GetValue(xi)
            ## TESTING
            # if patch_id == 1:
            #     # print(vgf.ControlGrid)
            #     # dpoint = cgf.GetDerivative(xi)
            #     # print("dpoint: " + str(dpoint[0]) + " " + str(dpoint[1]))
            #     du = multipatch_util.ComputeSpatialDerivatives(cgf, dgf, xi)
            #     # print("du: " + str(du))
            #     varphi_0 = psi[0] + du[2,0]
            #     varphi_1 = psi[1] + du[2,1]
            #     # print("psi: " + str(psi))
            #     # print("varphi (computed): " + str(varphi))
            #     # print("varphi (post-process): %.10e, %.10e" % (varphi_0, varphi_1))
            #     ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (r, disp[2], math.sqrt(psi[0]**2 + psi[1]**2), math.sqrt(varphi_0**2 + varphi_1**2), math.sqrt(varphi[0]**2 + varphi[1]**2)))
            # else:
            #     ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (r, disp[2], math.sqrt(psi[0]**2 + psi[1]**2), math.sqrt(varphi[0]**2 + varphi[1]**2), math.sqrt(varphi[0]**2 + varphi[1]**2)))
            ### END TESTING
            # print(vgf.ControlGrid)
            # dpoint = cgf.GetDerivative(xi)
            # print("dpoint: " + str(dpoint[0]) + " " + str(dpoint[1]))
            du = multipatch_util.ComputeSpatialDerivatives(cgf, dgf, xi)
            # print("du: " + str(du))
            varphi_0 = psi[0] + du[2,0]
            varphi_1 = psi[1] + du[2,1]
            # print("psi: " + str(psi))
            # print("varphi (computed): " + str(varphi))
            # print("varphi (post-process): %.10e, %.10e" % (varphi_0, varphi_1))
            ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (r, disp[2], psi[0], varphi_0, varphi[0]))
        ifile.close()

    def CheckDur(self, mpatch, r, dr):

        [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [r, 0.0, 0.0], [10, 10])
        patch = mpatch[patch_id].GetReference()
        cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
        point = cgf.GetValue(xi)
        dgf = patch.GridFunction(DISPLACEMENT)
        uz1 = dgf.GetValue(xi)[2]

        [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [r-dr, 0.0, 0.0], [10, 10])
        patch = mpatch[patch_id].GetReference()
        cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
        point = cgf.GetValue(xi)
        dgf = patch.GridFunction(DISPLACEMENT)
        uz2 = dgf.GetValue(xi)[2]

        print("dur = %.10e" % ((uz1 - uz2)/dr))

    def ComputeUError(self, model_part, analytical_solution):
        nom = 0.0
        denom = 0.0
        for element in model_part.Elements:
           if element.GetValue(IS_INACTIVE) == False:
               u = element.CalculateOnIntegrationPoints(DISPLACEMENT, model_part.ProcessInfo)
               J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
               Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
               W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
               for i in range(0, len(u)):
                   ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1])
                   nom = nom + (pow(u[i][2] - ana_u, 2)) * W[i][0] * J0[i][0]
                   denom = denom + (pow(ana_u, 2)) * W[i][0] * J0[i][0]
        return math.sqrt(nom/denom)

    def ComputePsiError(self, model_part, analytical_solution):
        nom = 0.0
        denom = 0.0
        for element in model_part.Elements:
           if element.GetValue(IS_INACTIVE) == False:
               psi = element.CalculateOnIntegrationPoints(ROTATION, model_part.ProcessInfo)
               J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
               Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
               W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
               for i in range(0, len(psi)):
                   ana_psi_r = analytical_solution.get_psi_r(Q[i][0], Q[i][1])
                   psi_r = math.sqrt(pow(psi[i][0], 2) + pow(psi[i][1], 2))
                   nom = nom + (pow(psi_r - ana_psi_r, 2)) * W[i][0] * J0[i][0]
                   denom = denom + (pow(ana_psi_r, 2)) * W[i][0] * J0[i][0]
        return math.sqrt(nom/denom)

    def Run(self, output=True):
        radius_bar = self.radius / self.thickness
        # create the scaled multipatch. On the scaled multipatch, the VARIABLE and TRUE_VARIABLE are not the same.
        mpatch = self.CreateMultiPatch(radius_bar, self.order)
        mpatch = self.Refine(mpatch, self.nsampling)
        # create the unscaled multipatch. On the unscaled multipatch, the VARIABLE and TRUE_VARIABLE are the same. They can be used interchangeably.
        mpatch_orig = self.CreateMultiPatch(self.radius, self.order)
        mpatch_orig = self.Refine(mpatch_orig, self.nsampling)
        print("########################################")
        mpatch.Enumerate()
        print(mpatch)
        # sys.exit(0)

        if output:
            mpatch_export.Export(mpatch, "multi_patch.m")

        patch_elems, mpatch_mp = self.CreateModel(mpatch)
        model_part = mpatch_mp.GetModelPart()

        #############ANALYSIS MODEL#######################################
        sim_params = model_iga_include.StaticParameters()
        sim_params["linear_solver"] = self.plinear_solver
        if output == False:
            sim_params["log_residuum"] = False
        model = model_iga_include.Model('multi_patch', os.getcwd()+"/", model_part, sim_params)
        model.InitializeModel()
        # time = 0.0
        # model.Solve(time, 0, 0, 0, 0)

        # fix displacement on the boundary
        tol = 1.0e-6
        for node in model.model_part.Nodes:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
        for cond in model.model_part.Conditions:
            for node in cond.GetNodes():
                node.Fix(DISPLACEMENT_Z)
                node.Fix(ROTATION_X)
                node.Fix(ROTATION_Y)

        # body force
        body_force = Vector(3)
        body_force[0] = 0
        body_force[1] = 0
        body_force[2] = self.face_load
        model.model_part.Properties[1].SetValue(BODY_FORCE, body_force)

        # analysis
        time = 1.0
        model.Solve(time, 0, 0, 0, 0)

        # preparation for post-processing
        transfer_util = BezierPostUtility()
        # transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())
        transfer_util.TransferVariablesToNodes(TRUE_SHEAR_ROTATION, model.model_part, SuperLUSolver())
        # true_rotation_control_values = transfer_util.TransferVariablesToNodalArray(TRUE_SHEAR_ROTATION, model.model_part, SuperLUSolver())

        true_rotation_patch_values = {}
        for patch_id, elems in patch_elems.items():
            true_rotation_control_values = transfer_util.TransferVariablesToNodalArray(TRUE_SHEAR_ROTATION, model.model_part, elems, SuperLUSolver())
            true_rotation_patch_values[patch_id] = true_rotation_control_values
        # print(true_rotation_patch_values)

        # ######Synchronize back the results ttransfer_util.TransferVariablesToNodeso multipatch
        # mpatch_mp.SynchronizeBackward(DISPLACEMENT)
        # mpatch_mp.SynchronizeBackward(ROTATION)
        # mpatch_mp.SynchronizeBackward(TRUE_SHEAR_ROTATION)
        # # mpatch_mp.SynchronizeBackward(THREED_STRESSES)
        ##################################################################

        if output:
            ## post processing

            patch_orig_elems, mpatch_mp_orig = self.CreateModel(mpatch_orig)
            model_part_orig = mpatch_mp_orig.GetModelPart()
            model_orig = model_iga_include.Model('multi_patch', os.getcwd()+"/", model_part_orig, model_iga_include.StaticParameters())
            model_orig.InitializeModel()

            for node_orig in model_part_orig.Nodes:
                node = model_part.Nodes[node_orig.Id]
                node_orig.SetSolutionStepValue(DISPLACEMENT, node.GetSolutionStepValue(DISPLACEMENT))
                node_orig.SetSolutionStepValue(ROTATION, node.GetSolutionStepValue(ROTATION)/self.thickness) # obtain the true rotation
                # node_orig.SetSolutionStepValue(TRUE_SHEAR_ROTATION, node.GetSolutionStepValue(TRUE_SHEAR_ROTATION))

            mpatch_mp_orig.SynchronizeBackward(DISPLACEMENT)
            mpatch_mp_orig.SynchronizeBackward(ROTATION)
            # mpatch_mp_orig.SynchronizeBackward(TRUE_SHEAR_ROTATION)
            mpatch_mp_orig.SynchronizeBackward(TRUE_SHEAR_ROTATION, true_rotation_patch_values)

            params_post = {}
            params_post['name'] = "multi_patch"
            params_post['division mode'] = "uniform"
            params_post['uniform division number'] = 40
        #    params_post['division mode'] = "non-uniform"
        #    params_post['division number u'] = 10
        #    params_post['division number v'] = 10
        #    params_post['division number w'] = 1
            params_post['variables list'] = [DISPLACEMENT, ROTATION, TRUE_SHEAR_ROTATION]
            dim = 2
            model_iga_include.PostMultiPatch(mpatch_orig, dim, time, params_post)

            ##################################################################

            E = self.young_modulus
            nu = self.poisson_ratio
            h = self.thickness
            G = E/(2*(1+nu))
            R = self.radius
            q = self.face_load

            analytical_solution = AnalyticalSolution(q, G, nu, h, R)

            wck = analytical_solution.get_kh_displacement(0.0, 0.0)
            wc = analytical_solution.get_displacement(0.0, 0.0)
            print("Analytical deflection at center: %.10e" % (wc))
            print("Kirchhoff solution: %.10e" % (wck))

            patch5 = mpatch_orig[5].GetReference()
            # cgf = patch5.GridFunction(CONTROL_POINT_COORDINATES)
            [stat, xi] = patch5.LocalCoordinates([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
            # print("stat: %d" % stat)
            dgf = patch5.GridFunction(DISPLACEMENT)
            disp = dgf.GetValue(xi)
            # print(disp)
            disp_error = abs(disp[2] - wc) / abs(wc) * 100
            print("Computed deflection at center: %.10e, error = %.10e %%" % (disp[2], disp_error))

            pgf = patch5.GridFunction(ROTATION)
            psi = pgf.GetValue(xi)
            print("Computed psi at center: %.10e, %.10e" % (psi[0], psi[1]))

            u_error = self.ComputeUError(model_part_orig, analytical_solution)
            print("U error: %.10e" % u_error)
            psi_error = self.ComputePsiError(model_part_orig, analytical_solution)
            print("Psi error: %.10e" % psi_error)

            self.ExportUr("ur_computed.txt", mpatch_orig, self.radius, npoints=20)
            analytical_solution.ExportUr("ur_analytical.txt", npoints=100)
            # self.CheckDur(mpatch_orig, self.radius, 1e-8)

        return model

    def TestSpatialDerivatives(self):
        radius_bar = 1.0
        mpatch = self.CreateMultiPatch(radius_bar, self.order)
        mpatch = self.Refine(mpatch, self.nsampling)
        print("########################################")
        mpatch.Enumerate()
        print(mpatch)
        # sys.exit(0)

        patch_elems, mpatch_mp = self.CreateModel(mpatch)
        model_part = mpatch_mp.GetModelPart()

        #############ANALYSIS MODEL#######################################
        sim_params = model_iga_include.StaticParameters()
        sim_params["linear_solver"] = self.plinear_solver
        model = model_iga_include.Model('multi_patch', os.getcwd()+"/", model_part, sim_params)
        model.InitializeModel()
        # time = 0.0
        # model.Solve(time, 0, 0, 0, 0)

        # fix displacement on the boundary
        tol = 1.0e-6
        for node in model.model_part.Nodes:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
        for cond in model.model_part.Conditions:
            for node in cond.GetNodes():
                node.Fix(DISPLACEMENT_Z)
                node.Fix(ROTATION_X)
                node.Fix(ROTATION_Y)

        # body force
        body_force = Vector(3)
        body_force[0] = 0
        body_force[1] = 0
        body_force[2] = self.face_load
        model.model_part.Properties[1].SetValue(BODY_FORCE, body_force)

        # analysis
        time = 1.0
        model.Solve(time, 0, 0, 0, 0)

        mpatch_mp.SynchronizeBackward(DISPLACEMENT)

        # [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [-0.1954919334,0.1954919334, 0.0], [10, 10])
        [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [-0.18, 0.1954919334, 0.0], [10, 10])
        patch = mpatch[patch_id].GetReference()
        cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
        point = cgf.GetValue(xi)
        print("patch_id: " + str(patch_id) + ", xi: " + str(xi) + ", point: " + str(point))
        if patch_id < 0:
            print("patch_id: %d" % patch_id)
            print("r: %e" % r)
            raise Exception("Fail computing the local coordinates")
        dgf = patch.GridFunction(DISPLACEMENT)
        disp = dgf.GetValue(xi)
        # pgf = patch.GridFunction(ROTATION)
        # psi = pgf.GetValue(xi)
        # vgf = patch.GridFunction(TRUE_SHEAR_ROTATION)
        # varphi = vgf.GetValue(xi)
        du = multipatch_util.ComputeSpatialDerivatives(cgf, dgf, xi)
        print("computed du: " + str(du))
        print(du[2,0])
        print(du[2,1])
