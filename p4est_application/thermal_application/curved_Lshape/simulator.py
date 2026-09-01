import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

import p4est_simulator

###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
def compute_L2_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            u = element.GetValuesOnIntegrationPoints(TEMPERATURE, process_info)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
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

###COMPUTE GLOBAL DISPLACEMENT (H1) ERROR###
def compute_H1_error(elements, solution, process_info):
    nom = 0.0
    denom = 0.0
    for element in elements:
        if element.Is(ACTIVE):
            du = element.GetValuesOnIntegrationPoints(TEMPERATURE_GRADIENT, process_info)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
            for i in range(0, len(du)):
                ana_du = solution.GetTemperatureGradientAt(Q[i][0], Q[i][1], Q[i][2])
                dist = 0.0
                norm_ana_du = 0.0
                for j in range(0, 3):
                    dist += pow(du[i][j] - ana_du[j], 2)
                    norm_ana_du += pow(ana_du[j], 2)
                nom = nom + dist * W[i][0] * J0[i][0]
                denom = denom + norm_ana_du * W[i][0] * J0[i][0]
    print("nom:", nom)
    print("denom:", denom)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

class Simulator:

    def __init__(self, params):
        self.params = params
        #########Some default parameters##################
        if "write_output_per_each_step" not in self.params:
            self.params["write_output_per_each_step"] = False
        if "compute_condition_number" not in self.params:
            self.params["compute_condition_number"] = False

        if "mesh_smoothing" not in self.params:
            self.params['mesh_smoothing'] = "none"

    ###SIMULATION DRIVER#############
    def Initialize(self, model):
        # register the boundary
        self.cx = 2.0
        self.cy = 0.0
        self.r1 = 1.0
        self.r2 = 2.0
        self.r3 = 3.0
        tol = 1.0e-6
        ########Create the B-Rep representing the boundary
        brep1 = CircularLevelSet(self.cx, self.cy, self.r1)
        brep2 = CircularLevelSet(self.cx, self.cy, self.r2)
        brep3 = CircularLevelSet(self.cx, self.cy, self.r3)
        self.params['boundaries'] = {}
        self.params['boundaries'][1] = brep1
        self.params['boundaries'][2] = brep2
        self.params['boundaries'][3] = brep3

        # create p4est model
        p4est_simulator.SetLayerIndex(model)
        p4est_simulator.SetBoundaryFlags(model)
        self.p4est_model = p4est_simulator.CreateP4estModel(model.model_part, self.params)

    def Update(self, model):
        # create the model_part
        [mp, layer_nodes_sets, layer_sets, layer_cond_sets, node_groups, element_assignments] = p4est_simulator.ConstructSystemModelPart(self.p4est_model, self.params)
        model.SetModelPart(mp)
        model.layer_nodes_sets = layer_nodes_sets
        model.layer_sets = layer_sets
        model.layer_cond_sets = layer_cond_sets
        model.node_groups = node_groups
        model.InitializeModel()

        model.solver.solver.MoveMeshFlag = False

    def Run(self, model, time=1.0, logging=True):
        ##################################################################
        ###  SIMULATION  #################################################
        ##################################################################
        compute_condition_number = self.params["compute_condition_number"]
        output = self.params['write_output_per_each_step']

        # if output:
        #     model.WriteOutput(0.0)

        if compute_condition_number:
            model.solver.solver.ReformDofSetAtEachStep = False
            eigen_solver = SpectraEigenvaluesSolver()

        for cond in model.model_part.Conditions:
            cond.SetValue(ACTIVATION_LEVEL, -1)

        # print(model.layer_nodes_sets['boundary'])
        for node_id in model.layer_nodes_sets['boundary']:
            node = model.model_part.Nodes[node_id]
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0.0)
            # print(f"Node {node.Id} TEMPERATURE is fixed")

        util = HeatStdProblem1Solution()
        for node_id in model.layer_nodes_sets['boundary']:
            node = model.model_part.Nodes[node_id]
            ana_temp = util.CalculateTemperatureOnNode(node)
            node.SetSolutionStepValue(TEMPERATURE, ana_temp)

        model.Solve(time, 0, 0, 0, 0)

        for node in model.model_part.Nodes:
            ana_temp = util.CalculateTemperatureOnNode(node)
            temp = node.GetSolutionStepValue(TEMPERATURE)
        #    if abs(ana_temp) < 1.0e-10:
        #        error = temp-ana_temp
        #    else:
        #        error = (temp-ana_temp)/abs(ana_temp)
            error = temp - ana_temp
            node.SetSolutionStepValue(TEMPERATURE_ERROR, error)
            node.SetSolutionStepValue(REFERENCE_TEMPERATURE, ana_temp)

        if output:
            for cond in model.model_part.Conditions:
                cond.SetValue(ACTIVATION_LEVEL, 0)
                cond.Set(ACTIVE, True)
            model.WriteOutput(time)

        l2_error = compute_L2_error(model.model_part.Elements, util, model.model_part.ProcessInfo)
        h1_error = compute_H1_error(model.model_part.Elements, util, model.model_part.ProcessInfo)
        model.l2_error = l2_error
        model.h1_error = h1_error
        print("Global L2 error: %.16e" % l2_error)
        print("Global H1 error: %.16e" % h1_error)

        if compute_condition_number:
            eigen_solver = SpectraEigenvaluesSolver()
            values = eigen_solver.SolveLargestSym(model.solver.solver.A, 3)
            print("largest eigenvalues:", values)
            values = eigen_solver.SolveSmallestSPD(model.solver.solver.A, SuperLUSolver(), 2)
            print("smallest eigenvalues:", values)
