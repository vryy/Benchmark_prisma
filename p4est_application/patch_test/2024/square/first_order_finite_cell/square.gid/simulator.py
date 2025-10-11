import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.FiniteCellApplication import *
from KratosMultiphysics.FiniteCellStructuralApplication import *

import finite_cell_simulator
from finite_cell_simulator import *

class Simulator(FiniteCellSimulator, object):

    def __init__(self, params):
        super(Simulator, self).__init__(params)
        ########Create the B-Rep representing the boundary
        self.brep = InverseLevelSet(CircularLevelSet(0.0, 0.0, 0.5))
        #########Some default parameters##################
        if "write_output_per_each_step" not in self.params:
            self.params["write_output_per_each_step"] = False
        #############variable_transfer_utility############
        #self.variable_transfer_utility = VariableTransferUtility(MKLPardisoSolver())

    ###SIMULATION DRIVER#############
    def Initialize(self, model_solid):
        ######## Extract the solid elements
        self.aux_util = FiniteCellAuxiliaryUtility()
        self.solid_elements = ElementsArray()
        for elem in model_solid.model_part.Elements:
            self.aux_util.AddElement(self.solid_elements, elem)
#        self.aux_util.GetElements(solid_elements, model_solid.model_part, model_solid.layer_sets['solid'])
        print("len(solid_elements):", len(self.solid_elements))

        super(Simulator, self).Initialize(model_solid, self.solid_elements)
        model_solid.solver.solver.MoveMeshFlag = False

        # self.variable_projection_utility = VariableProjectionUtility(MKLPardisoSolver())
        # self.variable_transfer_utility = VariableAdvancedTransferUtility(self.solid_elements, 0.31, 0.31, 0.31)

        # reset the displacements
        # this is to prevent weird problem of seldomly failed test
        for node in model_solid.model_part.Nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    def Run(self, model_solid):
        start_time = time_module.time()

        self.Initialize(model_solid)
        ##################################################################
        ###  SIMULATION  #################################################
        ##################################################################

        tol = 1.0e-6
        prescribed_nodes = []
        for node in model_solid.model_part.Nodes:
            if abs(node.X0) < tol:
                node.Fix(DISPLACEMENT_X)
            if abs(node.Y0) < tol:
                node.Fix(DISPLACEMENT_Y)
            if abs(node.X0 - 1.0) < tol:
                node.Fix(DISPLACEMENT_X)
                prescribed_nodes.append(node)
                print(node.Id)

        time = 100.0
        model_solid.Solve(time, 0, 0, 0, 0)
        if self.params['write_output_per_each_step']:
            model_solid.WriteOutput(time)

        time = 101.0
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.1)
        model_solid.Solve(time, 0, 0, 0, 0)

        end_time = time_module.time()
        print("analysis time: " + str(end_time-start_time) + " s")

        #self.variable_transfer_utility.TransferVariablesToNodes(model_solid.model_part, VON_MISES_STRESS)
#        self.variable_transfer_utility.TransferVariablesToNodes(VON_MISES_STRESS, model_solid.model_part.ProcessInfo)
        if self.params['write_output_per_each_step']:
            model_solid.WriteOutput(time)

        timer = Timer()
        print(timer)
