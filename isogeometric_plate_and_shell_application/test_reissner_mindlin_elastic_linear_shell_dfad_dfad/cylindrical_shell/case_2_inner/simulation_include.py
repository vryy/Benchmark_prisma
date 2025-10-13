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

# def ComputeL2Error(model_part, analytical_solution):
#     nom = 0.0
#     denom = 0.0
#     for element in model_part.Elements:
#         u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model_part.ProcessInfo)
#         J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
#         Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
#         W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
#         for i in range(0, len(u)):
#             a = math.atan2(Q[i][2], Q[i][1])
#             l = math.sqrt(Q[i][1]**2 + Q[i][2]**2)
#             c = Q[i][1] / l
#             s = Q[i][2] / l
#             ur = u[i][1]*c + u[i][2]*s
#             w = math.acos(c)
#             ana_u = analytical_solution.Interpolate(w)
#             nom += pow(ur - ana_u, 2) * W[i][0] * J0[i][0]
#             denom += pow(ana_u, 2) * W[i][0] * J0[i][0]
#     return nom / denom

def ComputeL2Error(model_part, analytical_solution):
    nom = 0.0
    denom = 0.0
    for element in model_part.Elements:
        u = element.GetValuesOnIntegrationPoints(TRUE_AVERAGE_DISPLACEMENT, model_part.ProcessInfo)
        # u = element.GetValuesOnIntegrationPoints(LOCAL_PROJECTED_DISPLACEMENT, model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
        for i in range(0, len(u)):
            a = math.atan2(Q[i][2], Q[i][1])
            l = math.sqrt(Q[i][1]**2 + Q[i][2]**2)
            c = Q[i][1] / l
            s = Q[i][2] / l
            ur = u[i][2]
            w = math.acos(c)
            ana_u = analytical_solution.Interpolate(w)
            nom += pow(ur - ana_u, 2) * W[i][0] * J0[i][0]
            denom += pow(ana_u, 2) * W[i][0] * J0[i][0]
    return math.sqrt(nom / denom)

class Model:

    def __init__(self, E, nu, R, L, h, p, order=2, nsampling=10, plinear_solver = SuperLUSolver(), drill_stiff=1e2):
        self.young_modulus = E
        self.poisson_ratio = nu
        self.radius = R
        self.length = L
        self.thickness = h
        self.pressure = p
        self.order = order
        self.nsampling = nsampling
        self.plinear_solver = plinear_solver
        self.drill_stiff = drill_stiff
        self.mode = 2

    def CreateMultiPatch(self, length, radius, order=3, mode=2):
        if mode == 1:
            ####### create arc 1
            arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'x', radius, 0.0, 180.0)
            arc1 = arc1_ptr.GetReference()

            ####### create arc 2
            arc2_ptr = geometry_factory.CreateSmallArc([length, 0.0, 0.0], 'x', radius, 0.0, 180.0)
            arc2 = arc2_ptr.GetReference()

        elif mode == 2:

            ####### create arc 1
            arc1_ptr = geometry_factory.CreateHalfCircle([0.0, 0.0, 0.0], 'x', radius)
            arc1 = arc1_ptr.GetReference()

            ####### create arc 2
            arc2_ptr = geometry_factory.CreateHalfCircle([length, 0.0, 0.0], 'x', radius)
            arc2 = arc2_ptr.GetReference()

        # create patch 2
        patch1_ptr = bsplines_patch_util.CreateLoftPatch(arc1, arc2)
        patch1 = patch1_ptr.GetReference()
        patch1.Id = 1

        # # print(patch2)
        # print(patch3)
        # print("line2:")
        # print(line2)

        ######create multipatch
        mpatch = MultiPatch2D()
        mpatch.AddPatch(patch1_ptr)

        if mode == 1:
            if order >= 2:
                multipatch_refine_util.DegreeElevate(mpatch[1], [order-2, order-1])
        elif mode == 2:
            if order >= 3:
                multipatch_refine_util.DegreeElevate(mpatch[1], [order-3, order-1])

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
        element_name = "ReissnerMindlinElasticLinearShellElementDFadDFadBezier2D3"

        mpatch_mp = MultiPatchModelPart2D(mpatch)

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
        prop.SetValue(YOUNG_MODULUS,        self.young_modulus )
        prop.SetValue(POISSON_RATIO,          self.poisson_ratio )
        prop.SetValue(THICKNESS,            self.thickness )
        prop.SetValue(DRILLING_STIFFNESS,            self.drill_stiff )

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

        mpatch_mp.EndModelPart()
        #    print(mpatch_mp)

        return [patch_elems, mpatch_mp]

    @staticmethod
    def ExportUr(filename, mpatch, radius, npoints=100):
        ifile = open(filename, "w")
        ifile.write("%-*s%s\n" % (20, "r", "w"))
        for i in range(0, npoints):
            alpha = float(i)/(npoints-1) * math.pi
            c = math.cos(alpha)
            s = math.sin(alpha)
            [patch_id, xi] = multipatch_util.LocalCoordinates(mpatch, [0.0, radius*c, radius*s], [10, 10])
            patch = mpatch[patch_id].GetReference()
            cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
            point = cgf.GetValue(xi)
            if patch_id < 0:
                print("patch_id: %d" % patch_id)
                print("r: %e" % r)
                raise Exception("Fail computing the local coordinates")
            dgf = patch.GridFunction(DISPLACEMENT)
            disp = dgf.GetValue(xi)
            ur = disp[1]*c + disp[2]*s
            ifile.write("%-*.10e%.10e\n" % (20, alpha, ur))
        ifile.close()

    @staticmethod
    def ExtractUrCheck(model_part, x_ref=None, export=True, filename="ur_check.txt"):
        all_values = []
        if x_ref == None:
            for elem in model_part.Elements:
                coords = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
                disp = elem.CalculateOnIntegrationPoints(TRUE_AVERAGE_DISPLACEMENT, model_part.ProcessInfo)
                rot = elem.CalculateOnIntegrationPoints(TRUE_ROTATION, model_part.ProcessInfo)
                # print(disp)
                for i in range(0, len(coords)):
                    all_values.append([coords[i], disp[i], rot[i]])
        else:
            for elem in model_part.Elements:
                coords = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
                disp = elem.CalculateOnIntegrationPoints(TRUE_AVERAGE_DISPLACEMENT, model_part.ProcessInfo)
                rot = elem.CalculateOnIntegrationPoints(TRUE_ROTATION, model_part.ProcessInfo)
                # print(disp)
                for i in range(0, len(coords)):
                    if abs(coords[i][0] - x_ref) < 1e-10:
                        all_values.append([coords[i], disp[i], rot[i]])
        if export:
            ifile = open(filename, "w")
            ifile.write("%-*s%-*s%-*s%-*s%-*s%s\n" % (20, "x", 20, "y", 20, "z", 20, "ut1", 20, "ut2", "uncheck"))
            for v in all_values:
                ifile.write("%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%.12e\n" % (20, v[0][0], 20, v[0][1], 20, v[0][2], 20, v[1][0], 20, v[1][1], v[1][2]))
            ifile.close()
        return all_values

    @staticmethod
    def ExtractForce(model_part, x_ref=None, export=True, filename="force.txt"):
        all_values = []
        if x_ref == None:
            for elem in model_part.Elements:
                coords = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
                m = elem.CalculateOnIntegrationPoints(LOCAL_BENDING_MOMENT, model_part.ProcessInfo)
                n = elem.CalculateOnIntegrationPoints(LOCAL_MEMBRANE_FORCE, model_part.ProcessInfo)
                q = elem.CalculateOnIntegrationPoints(LOCAL_SHEAR_FORCE, model_part.ProcessInfo)
                mc = elem.CalculateOnIntegrationPoints(LOCAL_BENDING_MOMENT_INTEGRAL_CHARACTERISTIC, model_part.ProcessInfo)
                nc = elem.CalculateOnIntegrationPoints(LOCAL_MEMBRANE_FORCE_INTEGRAL_CHARACTERISTIC, model_part.ProcessInfo)
                qc = elem.CalculateOnIntegrationPoints(LOCAL_SHEAR_FORCE_INTEGRAL_CHARACTERISTIC, model_part.ProcessInfo)
                for i in range(0, len(coords)):
                    all_values.append([coords[i], n[i], m[i], q[i], nc[i], mc[i], qc[i]])
        else:
            for elem in model_part.Elements:
                coords = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model_part.ProcessInfo)
                m = elem.CalculateOnIntegrationPoints(LOCAL_BENDING_MOMENT, model_part.ProcessInfo)
                n = elem.CalculateOnIntegrationPoints(LOCAL_MEMBRANE_FORCE, model_part.ProcessInfo)
                q = elem.CalculateOnIntegrationPoints(LOCAL_SHEAR_FORCE, model_part.ProcessInfo)
                mc = elem.CalculateOnIntegrationPoints(LOCAL_BENDING_MOMENT_INTEGRAL_CHARACTERISTIC, model_part.ProcessInfo)
                nc = elem.CalculateOnIntegrationPoints(LOCAL_MEMBRANE_FORCE_INTEGRAL_CHARACTERISTIC, model_part.ProcessInfo)
                qc = elem.CalculateOnIntegrationPoints(LOCAL_SHEAR_FORCE_INTEGRAL_CHARACTERISTIC, model_part.ProcessInfo)
                for i in range(0, len(coords)):
                    if abs(coords[i][0] - x_ref) < 1e-10:
                        all_values.append([coords[i], n[i], m[i], q[i], nc[i], mc[i], qc[i]])
        if export:
            ifile = open(filename, "w")
            ifile.write("%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%s\n" % (20, "x", 20, "y", 20, "z", 20, "n00", 20, "n01", 20, "n10", 20, "n11", 20, "m00", 20, "m01", 20, "m10", 20, "m11", 20, "q0", 20, "q1", 20, "nic00", 20, "nic01", 20, "nic10", 20, "nic11", 20, "mic00", 20, "mic01", 20, "mic10", 20, "mic11", 20, "qic0", "qic1"))
            for v in all_values:
                # print(v[1])
                ifile.write("%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%.12e\n" % (20, v[0][0], 20, v[0][1], 20, v[0][2], 20, v[1][0], 20, v[1][1], 20, v[1][2], 20, v[1][3], 20, v[2][0], 20, v[2][1], 20, v[2][2], 20, v[2][3], 20, v[3][0], 20, v[3][1], 20, v[4][0], 20, v[4][1], 20, v[4][2], 20, v[4][3], 20, v[5][0], 20, v[5][1], 20, v[5][2], 20, v[5][3], 20, v[6][0], v[6][1]))
            ifile.close()
        return all_values

    @staticmethod
    def ConvertUrCheck(values, export=True, filename='ur_check_polar.txt'):
        values_with_ang = []
        for v in values:
            y = v[0][1]
            z = v[0][2]
            a = math.atan2(y, z)
            values_with_ang.append([a, v[1][0], v[1][1], v[1][2], v[2][0], v[2][1], v[2][2]])
        sorted_values_with_ang = sorted(values_with_ang, key = lambda x: x[0])
        if export:
            ifile = open(filename, "w")
            ifile.write("%-*s%-*s%-*s%-*s%-*s%-*s%s\n" % (20, "a", 20, "ut1", 20, "ut2", 20, "uncheck", 20, "psit1", 20, "psit2", "psin"))
            for v in sorted_values_with_ang:
                ifile.write("%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%.12e\n" % (20, v[0], 20, v[1], 20, v[2], 20, v[3], 20, v[4], 20, v[5], v[6]))
            ifile.close()
        return sorted_values_with_ang

    @staticmethod
    def ConvertForce(values, export=True, filename='force_polar.txt'):
        values_with_ang = []
        for v in values:
            y = v[0][1]
            z = v[0][2]
            a = math.atan2(y, z)
            values_with_ang.append([a, v[1][0], v[1][1], v[1][2], v[1][3], v[2][0], v[2][1], v[2][2], v[2][3], v[3][0], v[3][1], v[4][0], v[4][1], v[4][2], v[4][3], v[5][0], v[5][1], v[5][2], v[5][3], v[6][0], v[6][1]])
        sorted_values_with_ang = sorted(values_with_ang, key = lambda x: x[0])
        if export:
            ifile = open(filename, "w")
            ifile.write("%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%s\n" % (20, "a", 20, "n00", 20, "n01", 20, "n10", 20, "n11", 20, "m00", 20, "m01", 20, "m10", 20, "m11", 20, "q0", 20, "q1", 20, "nic00", 20, "nic01", 20, "nic10", 20, "nic11", 20, "mic00", 20, "mic01", 20, "mic10", 20, "mic11", 20, "qic0", "qic1"))
            for v in sorted_values_with_ang:
                ifile.write("%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%-*.12e%.12e\n" % (20, v[0], 20, v[1], 20, v[2], 20, v[3], 20, v[4], 20, v[5], 20, v[6], 20, v[7], 20, v[8], 20, v[9], 20, v[10], 20, v[11], 20, v[12], 20, v[13], 20, v[14], 20, v[15], 20, v[16], 20, v[17], 20, v[18], 20, v[19], v[20]))
            ifile.close()
        return sorted_values_with_ang

    @staticmethod
    def ComputeMeshSize(mpatch):

        patch = mpatch[1].GetReference()
        ku = patch.FESpace().KnotU
        kus = list(set(ku))
        kus.sort()
        cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
        max_dist = -1e99
        for i in range(0, len(kus)-1):
            xi1 = [kus[i], 0.5]
            xi2 = [kus[i+1], 0.5]
            p1 = cgf.GetValue(xi1)
            p2 = cgf.GetValue(xi2)
            dist = math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
            if dist > max_dist:
                max_dist = dist
        return max_dist

    def Run(self, output=True, logging=True, disp=None, rot=None):
        length_bar = self.length / self.thickness
        radius_bar = self.radius / self.thickness
        p = self.pressure
        # create the scaled multipatch. On the scaled multipatch, the VARIABLE and TRUE_VARIABLE are not the same.
        mpatch = self.CreateMultiPatch(length_bar, radius_bar, self.order, self.mode)
        mpatch = self.Refine(mpatch, self.nsampling)
        # create the unscaled multipatch. On the unscaled multipatch, the VARIABLE and TRUE_VARIABLE are the same. They can be used interchangeably.
        # mpatch_orig = self.CreateMultiPatch(self.length, self.radius, self.order, self.mode)
        # mpatch_orig = self.Refine(mpatch_orig, self.nsampling)
        mpatch_orig = mpatch
        print("########################################")
        mpatch.Enumerate()
        print(mpatch)
        # sys.exit(0)

        self.h = self.ComputeMeshSize(mpatch)

        if output:
            mpatch_export.Export(mpatch, "one_patch.m")

        patch_elems, mpatch_mp = self.CreateModel(mpatch)
        model_part = mpatch_mp.GetModelPart()

        #############ANALYSIS MODEL#######################################
        sim_params = model_iga_include.StaticParameters()
        sim_params["linear_solver"] = self.plinear_solver
        sim_params["log_residuum"] = logging
        sim_params["convergence_criteria"] = "displacement"
        sim_params["rel_tol"] = 1e-10
        sim_params["abs_tol"] = 1e-10
        model = model_iga_include.Model('one_patch', os.getcwd()+"/", model_part, sim_params)
        model.InitializeModel()
        time = 0.0
        # model.Solve(time, 0, 0, 0, 0)

        # fix displacement on the boundary
        tol = 1e-6
        for node in model.model_part.Nodes:
            if (abs(node.X0) < tol) or (abs(node.X0 - length_bar) < tol):
                node.Fix(DISPLACEMENT_X)
                node.Fix(ROTATION_X)

            if (abs(node.Z0) < tol):
                node.Fix(DISPLACEMENT_Z)
                node.Fix(DISPLACEMENT_Y)
                node.Fix(DISPLACEMENT_X)
                node.Fix(ROTATION_X)
                node.Fix(ROTATION_Y)
                node.Fix(ROTATION_Z)

        # apply initial displacement and rotation
        if disp != None:
            for node in model.model_part.Nodes:
                node.SetSolutionStepValue(DISPLACEMENT_X, disp[node.Id][0])
                node.SetSolutionStepValue(DISPLACEMENT_Y, disp[node.Id][1])
                node.SetSolutionStepValue(DISPLACEMENT_Z, disp[node.Id][2])

        if rot != None:
            for node in model.model_part.Nodes:
                node.SetSolutionStepValue(ROTATION_X, rot[node.Id][0])
                node.SetSolutionStepValue(ROTATION_Y, rot[node.Id][1])
                node.SetSolutionStepValue(ROTATION_Z, rot[node.Id][2])

        ## pressure load
        for element in model.model_part.Elements:
            element.SetValue(POSITIVE_FACE_PRESSURE, 0.0)
            element.SetValue(NEGATIVE_FACE_PRESSURE, p)

        # analysis
        time = 1.0
        model.Solve(time, 0, 0, 0, 0)

        self.ndofs = model.solver.solver.builder_and_solver.GetEquationSystemSize()

        # # preparation for post-processing
        # transfer_util = BezierPostUtility()
        # # transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())
        # transfer_util.TransferVariablesToNodes(TRUE_SHEAR_ROTATION, model.model_part, SuperLUSolver())
        # # true_rotation_control_values = transfer_util.TransferVariablesToNodalArray(TRUE_SHEAR_ROTATION, model.model_part, SuperLUSolver())

        # true_rotation_patch_values = {}
        # for patch_id, elems in patch_elems.items():
        #     true_rotation_control_values = transfer_util.TransferVariablesToNodalArray(TRUE_SHEAR_ROTATION, model.model_part, elems, SuperLUSolver())
        #     true_rotation_patch_values[patch_id] = true_rotation_control_values
        # # print(true_rotation_patch_values)

        # # ######Synchronize back the results ttransfer_util.TransferVariablesToNodeso multipatch
        # # mpatch_mp.SynchronizeBackward(DISPLACEMENT)
        # # mpatch_mp.SynchronizeBackward(ROTATION)
        # # mpatch_mp.SynchronizeBackward(TRUE_SHEAR_ROTATION)
        # # # mpatch_mp.SynchronizeBackward(THREED_STRESSES)
        # ##################################################################

        if output:
            ## post processing

            patch_orig_elems, mpatch_mp_orig = self.CreateModel(mpatch_orig)
            model_part_orig = mpatch_mp_orig.GetModelPart()
            model_orig = model_iga_include.Model('one_patch', os.getcwd()+"/", model_part_orig, model_iga_include.StaticParameters())
            model_orig.InitializeModel()

            for node_orig in model_part_orig.Nodes:
                node = model_part.Nodes[node_orig.Id]
                node_orig.SetSolutionStepValue(DISPLACEMENT, node.GetSolutionStepValue(DISPLACEMENT))
                rot_orig = Vector(3)
                rot = node.GetSolutionStepValue(ROTATION)
                rot_orig[0] = rot[0] / self.thickness
                rot_orig[1] = rot[1] / self.thickness
                rot_orig[2] = rot[2] / self.thickness
                node_orig.SetSolutionStepValue(ROTATION, rot_orig) # obtain the true rotation

            mpatch_mp_orig.SynchronizeBackward(DISPLACEMENT)
            mpatch_mp_orig.SynchronizeBackward(ROTATION)

            params_post = {}
            params_post['name'] = "one_patch"
            params_post['base condition name'] = "DummyLineCondition"
            params_post['division mode'] = "uniform"
            # params_post['uniform division number'] = 40
            params_post['division mode'] = "non-uniform"
            params_post['division number u'] = 20
            params_post['division number v'] = 5
            params_post['division number w'] = 1
            params_post['variables list'] = [DISPLACEMENT, ROTATION]
            dim = 2
            model_iga_include.PostMultiPatch(mpatch_orig, dim, time, params_post)

            ##################################################################

            # E = self.young_modulus
            # nu = self.poisson_ratio
            # h = self.thickness
            # G = E/(2*(1+nu))
            # R = self.radius
            # q = self.face_load

            # analytical_solution = AnalyticalSolution(q, G, nu, h, R)

            # wck = analytical_solution.get_kh_displacement(0.0, 0.0)
            # wc = analytical_solution.get_displacement(0.0, 0.0)
            # print("Analytical deflection at center: %.10e" % (wc))
            # print("Kirchhoff solution: %.10e" % (wck))

            # patch5 = mpatch_orig[5].GetReference()
            # # cgf = patch5.GridFunction(CONTROL_POINT_COORDINATES)
            # [stat, xi] = patch5.LocalCoordinates([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
            # # print("stat: %d" % stat)
            # dgf = patch5.GridFunction(DISPLACEMENT)
            # disp = dgf.GetValue(xi)
            # # print(disp)
            # disp_error = abs(disp[2] - wc) / abs(wc) * 100
            # print("Computed deflection at center: %.10e, error = %.10e %%" % (disp[2], disp_error))

            # pgf = patch5.GridFunction(ROTATION)
            # psi = pgf.GetValue(xi)
            # print("Computed psi at center: %.10e, %.10e" % (psi[0], psi[1]))

            # u_error = self.ComputeUError(model_part_orig, analytical_solution)
            # print("U error: %.10e" % u_error)
            # psi_error = self.ComputePsiError(model_part_orig, analytical_solution)
            # print("Psi error: %.10e" % psi_error)

        if logging and output:
            self.ExportUr("ur_computed.txt", mpatch_orig, self.radius / self.thickness, npoints=20)
            # self.ExportUrCheck("ur_check_computed.txt", model.model_part)
            middle_values = self.ExtractUrCheck(model.model_part, x_ref=0.5*self.length, export=False)
            self.ConvertUrCheck(middle_values, export=True, filename="ur_check_polar.txt")
            # analytical_solution.ExportUr("ur_analytical.txt", npoints=100)
            # # self.CheckDur(mpatch_orig, self.radius, 1e-8)

            middle_force_values = self.ExtractForce(model.model_part, x_ref=0.5*self.length, export=False)
            self.ConvertForce(middle_force_values, export=True, filename="force_polar.txt")

        return model
