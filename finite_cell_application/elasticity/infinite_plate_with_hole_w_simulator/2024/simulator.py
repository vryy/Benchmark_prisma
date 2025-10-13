import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.FiniteCellApplication import *
from KratosMultiphysics.FiniteCellStructuralApplication import *

import finite_cell_simulator
from finite_cell_simulator import *

class Simulator(FiniteCellSimulator, object):

    def __init__(self, params):
        super(Simulator, self).__init__(params)
        ########Create the B-Rep representing the boundary
        self.brep = InverseLevelSet(CircularLevelSet(0.0, 0.0, 1.0))
        #########Some default parameters##################
        if "write_output_per_each_step" not in self.params:
            self.params["write_output_per_each_step"] = False
        #############variable_transfer_utility############
        self.variable_transfer_utility = VariableTransferUtility(MKLPardisoSolver())

    ###SIMULATION DRIVER#############
    def Initialize(self, model_solid):
        ######## Extract the solid elements
        self.aux_util = FiniteCellAuxiliaryUtility()
        self.solid_elements = ElementsArray()
        for elem in model_solid.model_part.Elements:
            self.aux_util.AddElement(self.solid_elements, elem)
#        self.aux_util.GetElements(solid_elements, model_solid.model_part, model_solid.layer_sets['solid'])
        print("len(solid_elements):", len(self.solid_elements))

        if (self.params["quadrature_method"] == "quadtree") or (self.params["quadrature_method"] == "moment-fit quadtree"):
            self.all_elements = self.solid_elements
        elif self.params["quadrature_method"] == "moment-fit subcell":
            self.all_elements = ElementsArray()
#            print("len(solid_elements):", len(self.solid_elements))
#            for elem in self.proper_cut_elems:
#                print elem.Id
            for elem in self.solid_elements:
                if elem.Id not in self.proper_cut_elems:
                    self.aux_util.AddElement(self.all_elements, elem)
            print("len(all_subcell_elems):", len(self.all_subcell_elems))
            for elem in self.all_subcell_elems:
                self.aux_util.AddElement(self.all_elements, elem)
#            print("len(fict_elems):", len(self.fict_elems))
#            for elem in self.fict_elems:
#                self.aux_util.AddElement(self.all_elements, elem)
            print("len(all_elements):", len(self.all_elements))

        super(Simulator, self).Initialize(model_solid, self.solid_elements)
        model_solid.solver.solver.MoveMeshFlag = False

        self.variable_projection_utility = VariableProjectionUtility(self.solid_elements, MKLPardisoSolver())
        # self.variable_transfer_utility = VariableAdvancedTransferUtility(self.solid_elements, 0.31, 0.31, 0.31)

        E = model_solid.model_part.Properties[1].GetValue(YOUNG_MODULUS)
        nu = model_solid.model_part.Properties[1].GetValue(POISSON_RATIO)
        self.G = E/(2*(1+nu))
        self.kappa = (3-nu)/(1+nu) # plane stress

    def Run(self, model_solid):
        start_time = time_module.time()

        self.Initialize(model_solid)
        ##################################################################
        ###  SIMULATION  #################################################
        ##################################################################

        self.P = 10.0
        self.ri = 1.0
        FunctionR3Rn().Assign(LOAD_FUNCTION, LoadFunctionR3RnPlateWithTheHoleX(self.P, self.ri), model_solid.model_part.Properties[2])
        FunctionR3Rn().Assign(LOAD_FUNCTION, LoadFunctionR3RnPlateWithTheHoleY(self.P, self.ri), model_solid.model_part.Properties[3])

        time = 1.0
        model_solid.Solve(time, 0, 0, 0, 0)

        end_time = time_module.time()
        print("analysis time: " + str(end_time-start_time) + " s")

        #self.variable_transfer_utility.TransferVariablesToNodes(model_solid.model_part, VON_MISES_STRESS)
#        self.variable_transfer_utility.TransferVariablesToNodes(VON_MISES_STRESS, model_solid.model_part.ProcessInfo)

        if self.params['write_output_per_each_step']:
            model_solid.WriteOutput(time)

    def ComputeError(self, model_solid):

        ###COMPUTE THE STRAIN ENERGY USING ANALYTICAL SOLUTION ###
        analytical_solution = self.params["analytical_solution"]

        a = 4.0
        f_se = analytical_solution.get_strain_energy_function(self.P, self.ri, self.G, self.kappa)

        #print(f_se.GetValue(0.0, 4.0))
        #print(f_se.GetValue(2.0, 2.0))
        #print(f_se.GetValue(0.3, 1.2))

        strain_energy_exact = 0.0
        for element in model_solid.model_part.Elements:
            if element.GetValue(IS_INACTIVE) == False:
                strain_energy_exact = strain_energy_exact + f_se.Integrate(element, 2)
        print("strain_energy_exact: %.15e" % strain_energy_exact)

        ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
        nom = 0.0
        denom = 0.0
        for element in self.all_elements:
            if element.GetValue(IS_INACTIVE) == False:
                u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model_solid.model_part.ProcessInfo)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_solid.model_part.ProcessInfo)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_solid.model_part.ProcessInfo)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_solid.model_part.ProcessInfo)
                for i in range(0, len(u)):
                    ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1], self.P, self.ri, self.G, self.kappa)
                    nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2)) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
        l2_error = math.sqrt(nom / denom)
        print("Global displacement (L2) error: %.15e" % l2_error)

        ###COMPUTE GLOBAL STRESS (H1) ERROR###
        nom = 0.0
        denom = 0.0
        for element in self.all_elements:
            if element.GetValue(IS_INACTIVE) == False:
                o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model_solid.model_part.ProcessInfo)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_solid.model_part.ProcessInfo)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_solid.model_part.ProcessInfo)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_solid.model_part.ProcessInfo)
#                print("element " + str(element.Id) + " is counted " + str(len(W)) + " integration_points")
#                print("W:", W)
#                print("o:", o)
                for i in range(0, len(o)):
                    ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], self.P, self.ri, self.G, self.kappa)
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        h1_error = math.sqrt(nom / denom)
        print("Global stress (H1) error: %.15e" % h1_error)

        return strain_energy_exact, l2_error, h1_error

    def WriteVMS(self, ang, vms, filename):
        # sort the angles
        ang_dict = {}
        vms_dict = {}
        for i in range(0, len(ang)):
            ang_dict[i] = ang[i]
            vms_dict[i] = vms[i]
        sorted_ang = []
        sorted_vms = []
        # for key, value in sorted(ang_dict.iteritems(), key=lambda (k,v): (v,k)):
        for key, value in sorted(ang_dict.items(), key=lambda kv: (kv[1], kv[0])):
            sorted_ang.append(value)
            sorted_vms.append(vms[key])
        fid = open(filename, "w")
        fid.write("angle\tvon_mises_stress\n")
        for i in range(0, len(sorted_ang)):
            fid.write(str(sorted_ang[i]) + "\t" + str(sorted_vms[i]) + "\n")
        fid.close()

    def GetAnalyticalVMSOnBoundary(self, n, filename):
        analytical_solution = self.params["analytical_solution"]
        ang = []
        vms = []
        for i in range(0, n+1):
            a = (i*math.pi/2) / n
            ang.append(a)
            x = self.ri * math.cos(a)
            y = self.ri * math.sin(a)
            stress_3d = analytical_solution.get_stress_3d(x, y, self.P, self.ri, self.G, self.kappa)
            p = (stress_3d[0] + stress_3d[1] + stress_3d[2]) / 3
            o_xx = stress_3d[0] - p
            o_yy = stress_3d[1] - p
            o_zz = stress_3d[2] - p
            o_xy = stress_3d[3]
            o_yz = stress_3d[4]
            o_xz = stress_3d[5]
            q = math.sqrt(1.5*(o_xx*o_xx + o_yy*o_yy + o_zz*o_zz + 2*o_xy*o_xy + 2*o_yz*o_yz + 2*o_xz*o_xz))
            vms.append(q)
        self.WriteVMS(ang, vms, filename)

    def GetVMSOnBoundary(self, model_part, tol, filename):
        ###EXTRACT THE VON MISES STRESS ALONG THE BOUNDARY for the post model_part###
        ang = []
        vms = []
        for node in model_part.Nodes:
            r = math.sqrt(math.pow(node.X0, 2) + math.pow(node.Y0, 2))
            if abs(r - self.ri) < tol:
#                print("node ", node.Id, node.X0, node .Y0)
                ang.append(math.acos(node.X0/r))
                vms.append(node.GetSolutionStepValue(VON_MISES_STRESS))
        self.WriteVMS(ang, vms, filename)

    def GetVMSOnBoundaryByInterpolation(self, model_part, npoint, filename):
        util = ImmersedBoundaryUtility()
        util.InitializeBinning(model_part.Elements, 0.11, 0.11, 0.11)
        ###EXTRACT THE VON MISES STRESS ALONG THE BOUNDARY###
        ang = []
        vms = []
        P = Vector(3)
        P[2] = 0.0
        for i in range(0, npoint+1):
            a = (0.5*math.pi)*i/npoint
            P[0] = self.ri * math.cos(a)
            P[1] = self.ri * math.sin(a)
            ang.append(a)
            vms.append(util.GetValueOnPoint(VON_MISES_STRESS, P, model_part.Elements))
        self.WriteVMS(ang, vms, filename)
