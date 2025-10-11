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

    def __init__(self, q, G, nu, h, L):
        self.q = q
        self.G = G
        self.nu = nu
        self.h = h
        self.L = L

        self.sigma = self.nu/(1-self.nu)
        self.Eb = (self.sigma + 1.0)/6
        self.Es = 5.0/6

    def get_displacement(self, x, y):
        xb = x / self.h
        Lb = self.L / self.h
        fb = self.q*self.h/self.G
        wc = fb/(self.Eb)*((xb**4)/24 - Lb*(xb**3)/6 + 0.5*((Lb**2)/2 + self.sigma/10)*(xb**2)) \
           + fb/(self.Es)*(-0.5*(xb**2) + Lb*xb) \
           - self.sigma/60*fb/self.Eb * (-0.5*(xb**2) + Lb*xb - (Lb**2)/2 - self.sigma/10)
        return wc

    def get_psi(self, x, y):
        xb = x / self.h
        Lb = self.L / self.h
        fb = self.q*self.h/self.G
        psi = fb/(self.Eb)*(-(xb**3)/6 + 0.5*Lb*(xb**2) - (0.5*(Lb**2)+self.sigma/10)*xb)
        return psi/self.h

    def get_varphi(self, x, y):
        xb = x / self.h
        Lb = self.L / self.h
        fb = self.q*self.h/self.G
        varphi = fb/self.Es*(Lb-xb)
        return varphi/self.h

    def get_be_displacement(self, x, y):
        xb = x / self.h
        Lb = self.L / self.h
        fb = self.q*self.h/self.G
        wc = fb/(self.Eb)*((xb**4)/24 - Lb*(xb**3)/6 + 0.5*(Lb**2)/2*(xb**2))
        return wc

    def get_be_psi(self, x, y):
        xb = x / self.h
        Lb = self.L / self.h
        fb = self.q*self.h/self.G
        psi = fb/(self.Eb)*(-(xb**3)/6 + 0.5*Lb*(xb**2) - (0.5*(Lb**2))*xb)
        return psi/self.h

    def ExportUx(self, filename, npoints=100):
        ifile = open(filename, "w")
        ifile.write("x\t\tw\tw_be\t\tpsi\t\tvarphi\t\tpsi_be\n")
        for i in range(0, npoints):
            xi = float(i)/(npoints-1)
            x = xi*self.L
            ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (x, self.get_displacement(x, 0.0), self.get_be_displacement(x, 0.0), self.get_psi(x, 0.0), self.get_varphi(x, 0.0), self.get_be_psi(x, 0.0)))
        ifile.close()

class Model:

    def __init__(self, E, nu, L, D, h, f, order=2, nsampling=[10, 10], plinear_solver = MKLPardisoSolver()):
        self.young_modulus = E
        self.poisson_ratio = nu
        self.length = L
        self.width = D
        self.thickness = h
        self.face_load = f
        self.order = order
        self.nsampling = nsampling
        self.plinear_solver = plinear_solver

    def CreateMultiPatch(self, l, w, order=2):
        # create patch 1
        p1 = [0.0, 0.0, 0.0]
        p2 = [l, 0.0, 0.0]
        p3 = [l, w, 0.0]
        p4 = [0.0, w, 0.0]
        patch1_ptr = geometry_factory.CreateParallelogram(p1, p2, p3, p4)
        patch1 = patch1_ptr.GetReference()
        patch1.Id = 1

        ######create multipatch
        mpatch = MultiPatch2D()
        mpatch.AddPatch(patch1_ptr)

        if order >= 2:
            multipatch_refine_util.DegreeElevate(mpatch[1], [order-1, order-1])

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

    def CreateModel(self, mpatch):
        mpatch_util = MultiPatchUtility()
        element_name = "ReissnerMindlinElasticPlateElementBezier2D"

        mpatch_mp = MultiPatchModelPart2D(mpatch)

        mpatch_mp.BeginModelPart()
        model_part = mpatch_mp.GetModelPart()
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( model_part )
        model_part.AddNodalSolutionStepVariable(THREED_STRESSES)
        model_part.AddNodalSolutionStepVariable(TRUE_SHEAR_ROTATION)
        model_part.AddNodalSolutionStepVariable(TRUE_AVERAGE_DISPLACEMENT)

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

        for patch_ptr in mpatch.Patches():
            patch = patch_ptr.GetReference()
            patch.LayerIndex = 1
            #        print(patch)

            ## add volume elements
            last_elem_id = mpatch_util.GetLastElementId(model_part)
            elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)
            # print(len(elems))

        mpatch_mp.EndModelPart()
        #    print(mpatch_mp)

        return mpatch_mp

    def ExportUx(self, filename, mpatch, L, npoints=100):
        ifile = open(filename, "w")
        ifile.write("x\tw\tpsi\tvarphi\txa\n")
        for i in range(0, npoints):
            xi = float(i)/(npoints-1)
            x = xi*L
            # [patch_id, xi] = mpatch.LocalCoordinates([x, 0.0, 0.0], [xi, 0.0, 0.0])
            [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [x, 0.5*self.width, 0.0], [10, 10])
            if patch_id < 0:
                print("patch_id: %d" % patch_id)
                print("x: %e" % x)
                raise Exception("Fail computing the local coordinates")
            patch = mpatch[patch_id].GetReference()
            dgf = patch.GridFunction(DISPLACEMENT)
            disp = dgf.GetValue(xi)
            pgf = patch.GridFunction(ROTATION)
            psi = pgf.GetValue(xi)
            vgf = patch.GridFunction(TRUE_SHEAR_ROTATION)
            varphi = vgf.GetValue(xi)
            xagf = patch.GridFunction(TRUE_AVERAGE_DISPLACEMENT)
            xa = xagf.GetValue(xi)
            ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (x, disp[2], psi[0], varphi[0], xa[2]))
        ifile.close()

    # def ComputeL2Error(self, model_part, analytical_solution):
    #     nom = 0.0
    #     denom = 0.0
    #     for element in model_part.Elements:
    #        if element.GetValue(IS_INACTIVE) == False:
    #            u = element.CalculateOnIntegrationPoints(DISPLACEMENT, model_part.ProcessInfo)
    #            J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
    #            Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
    #            W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
    #            for i in range(0, len(u)):
    #                ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1])
    #                nom = nom + (pow(u[i][2] - ana_u, 2)) * W[i][0] * J0[i][0]
    #                denom = denom + (pow(ana_u, 2)) * W[i][0] * J0[i][0]
    #     return math.sqrt(nom/denom)

    # def ComputeH1Error(self, model_part, analytical_solution):
    #     nom = 0.0
    #     denom = 0.0
    #     for element in model_part.Elements:
    #        if element.GetValue(IS_INACTIVE) == False:
    #            psi = element.CalculateOnIntegrationPoints(ROTATION, model_part.ProcessInfo)
    #            J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
    #            Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
    #            W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
    #            for i in range(0, len(psi)):
    #                ana_psi_r = analytical_solution.get_psi_r(Q[i][0], Q[i][1])
    #                psi_r = math.sqrt(pow(psi[i][0], 2) + pow(psi[i][1], 2))
    #                nom = nom + (pow(psi_r - ana_psi_r, 2)) * W[i][0] * J0[i][0]
    #                denom = denom + (pow(ana_psi_r, 2)) * W[i][0] * J0[i][0]
    #     return math.sqrt(nom/denom)

    def Run(self, output=True):
        L_bar = self.length / self.thickness
        D_bar = self.width / self.thickness
        # create the scaled multipatch. On the scaled multipatch, the VARIABLE and TRUE_VARIABLE are not the same.
        mpatch = self.CreateMultiPatch(L_bar, D_bar, self.order)
        mpatch = self.Refine(mpatch, self.nsampling)
        # create the unscaled multipatch. On the unscaled multipatch, the VARIABLE and TRUE_VARIABLE are the same. They can be used interchangeably.
        mpatch_orig = self.CreateMultiPatch(self.length, self.width, self.order)
        mpatch_orig = self.Refine(mpatch_orig, self.nsampling)
        print("########################################")
        mpatch.Enumerate()
        print(mpatch)
        # sys.exit(0)

        if output:
            mpatch_export.Export(mpatch, "one_patch.m")

        mpatch_mp = self.CreateModel(mpatch)
        model_part = mpatch_mp.GetModelPart()

        #############ANALYSIS MODEL#######################################
        sim_params = model_iga_include.StaticParameters()
        sim_params["linear_solver"] = self.plinear_solver
        if output == False:
            sim_params["log_residuum"] = False
        model = model_iga_include.Model('one_patch', os.getcwd()+"/", model_part, sim_params)
        model.InitializeModel()
        # time = 0.0
        # model.Solve(time, 0, 0, 0, 0)

        # fix displacement on the left side
        tol = 1.0e-6
        for node in model.model_part.Nodes:
            if abs(node.X0) < tol:
                node.Fix(DISPLACEMENT_Z)
                node.Fix(ROTATION_X)
                node.Fix(ROTATION_Y)

        # body force
        body_force = Vector(3)
        body_force[0] = 0
        body_force[1] = 0
        body_force[2] = self.face_load
        model.model_part.Properties[1].SetValue(BODY_FORCE, body_force)

        time = 1.0
        model.Solve(time, 0, 0, 0, 0)

        transfer_util = BezierPostUtility()
        # transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, MKLPardisoSolver())
        transfer_util.TransferVariablesToNodes(TRUE_SHEAR_ROTATION, model.model_part, MKLPardisoSolver())
        transfer_util.TransferVariablesToNodes(TRUE_AVERAGE_DISPLACEMENT, model.model_part, MKLPardisoSolver())

        ######Synchronize back the results to multipatch
        mpatch_mp.SynchronizeBackward(DISPLACEMENT)
        mpatch_mp.SynchronizeBackward(ROTATION)
        mpatch_mp.SynchronizeBackward(TRUE_SHEAR_ROTATION)
        mpatch_mp.SynchronizeBackward(TRUE_AVERAGE_DISPLACEMENT)
        # mpatch_mp.SynchronizeBackward(THREED_STRESSES)
        ##################################################################

        if output:
            ## post processing

            mpatch_mp_orig = self.CreateModel(mpatch_orig)
            model_part_orig = mpatch_mp_orig.GetModelPart()
            model_orig = model_iga_include.Model('one_patch', os.getcwd()+"/", model_part_orig, model_iga_include.StaticParameters())
            model_orig.InitializeModel()

            for node_orig in model_part_orig.Nodes:
                node = model_part.Nodes[node_orig.Id]
                node_orig.SetSolutionStepValue(DISPLACEMENT, node.GetSolutionStepValue(DISPLACEMENT))
                node_orig.SetSolutionStepValue(ROTATION, node.GetSolutionStepValue(ROTATION)/self.thickness) # obtain the true rotation
                node_orig.SetSolutionStepValue(TRUE_SHEAR_ROTATION, node.GetSolutionStepValue(TRUE_SHEAR_ROTATION))
                node_orig.SetSolutionStepValue(TRUE_AVERAGE_DISPLACEMENT, node.GetSolutionStepValue(TRUE_AVERAGE_DISPLACEMENT))

            mpatch_mp_orig.SynchronizeBackward(DISPLACEMENT)
            mpatch_mp_orig.SynchronizeBackward(ROTATION)
            mpatch_mp_orig.SynchronizeBackward(TRUE_SHEAR_ROTATION)
            mpatch_mp_orig.SynchronizeBackward(TRUE_AVERAGE_DISPLACEMENT)

            params_post = {}
            params_post['name'] = "one_patch"
            params_post['division mode'] = "uniform"
            params_post['uniform division number'] = 40
        #    params_post['division mode'] = "non-uniform"
        #    params_post['division number u'] = 10
        #    params_post['division number v'] = 10
        #    params_post['division number w'] = 1
            params_post['variables list'] = [DISPLACEMENT, ROTATION, TRUE_SHEAR_ROTATION, TRUE_AVERAGE_DISPLACEMENT]
            dim = 2
            model_iga_include.PostMultiPatch(mpatch_orig, dim, time, params_post)

            ##################################################################

            E = self.young_modulus
            nu = self.poisson_ratio
            h = self.thickness
            G = E/(2*(1+nu))
            q = self.face_load

            analytical_solution = AnalyticalSolution(q, G, nu, h, self.length)

            # wck = analytical_solution.get_kh_displacement(0.0, 0.0)
            wc = analytical_solution.get_displacement(self.length, 0.5*self.width)
            print("Analytical deflection at middle of right edge: %.10e" % (wc))
            # print("Kirchhoff solution: %.10e" % (wck))

            patch1 = mpatch_orig[1].GetReference()
            # cgf = patch5.GridFunction(CONTROL_POINT_COORDINATES)
            [stat, xi] = patch1.LocalCoordinates([self.length, 0.5*self.width, 0.0], [0.0, 0.0, 0.0])
            # print("stat: %d" % stat)
            dgf = patch1.GridFunction(DISPLACEMENT)
            disp = dgf.GetValue(xi)
            # print(disp)
            print("Computed deflection at middle of right edge: %.10e" % (disp[2]))

            # pgf = patch5.GridFunction(ROTATION)
            # psi = pgf.GetValue(xi)
            # print("Computed psi at center: %.10e, %.10e" % (psi[0], psi[1]))

            # l2_error = self.ComputeL2Error(model_part_orig, analytical_solution)
            # print("L2 error: %.10e" % l2_error)
            # h1_error = self.ComputeH1Error(model_part_orig, analytical_solution)
            # print("H1 error: %.10e" % h1_error)

            self.ExportUx("w_computed.txt", mpatch_orig, self.length, npoints=20)
            analytical_solution.ExportUx("w_analytical.txt", npoints=100)

        return model
