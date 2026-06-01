import os
import math
import time as time_module
import pprint
import six
from KratosMultiphysics import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.MmgApplication import *

## refinement process that refine everything
class MmgRefinementProcessAll():
    def __init__(self, mesh_size):
        self.h = mesh_size

    def Initialize(self, mmg_model):
        mmg_model.SetValue(MMG_HAUSDORFF_DISTANCE, self.h)

    def Execute(self, model_part):
        for node in model_part.Nodes:
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, self.h)

## refinement process that refines only the nodes specified by the refine vector
class MmgRefinementProcessBasedOnRefineVector():
    def __init__(self, refine_vector):
        self.refine_vector = refine_vector

    def Initialize(self, mmg_model):
        pass # DO NOTHING

    def Execute(self, model_part):
        for node_id in self.refine_vector:
            node = model_part.Nodes[node_id]
            h = node.GetSolutionStepValue(NODAL_MMG_SCALAR_METRIC)
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 0.5*h)

## create the Mmg model from the model_part
def CreateMmgModel(model_part, params):
    mmg_model = Mmg2DModel()

    # refinement
    if 'initial_refinement_process' in params:
        if params['initial_refinement_process'] != None:
            params['initial_refinement_process'].Initialize(mmg_model)
            params['initial_refinement_process'].Execute(model_part)

    # mmg_model.SetValue(MMG_GRADATION, 1.1)
    if 'hausdorff_distance' in params:
        mmg_model.SetValue(MMG_HAUSDORFF_DISTANCE, params['hausdorff_distance'])
    mmg_model.SetIParam(MMG2D_Param.IPARAM_verbose, 10)
    mmg_model.Initialize(model_part)
    mmg_model.Remesh()

    return mmg_model

## create the model_part out from mmg_model
def ConstructSystemModelPart(mmg_model, params):
    # start to construct the Kratos model_part
    mmg_model.BeginModelPart()

    mp = mmg_model.GetModelPart()
    mp.SetBufferSize(2)
    mp.AddNodalSolutionStepVariable(TEMPERATURE)
    mp.AddNodalSolutionStepVariable(LAGRANGE_TEMPERATURE)
    mp.AddNodalSolutionStepVariable(TEMPERATURE_ERROR)
    mp.AddNodalSolutionStepVariable(REFERENCE_TEMPERATURE)
    mp.AddNodalSolutionStepVariable(NODAL_MMG_SCALAR_METRIC)

    # create nodes
    mmg_model.CreateNodes()
    import structural_solver_advanced
    structural_solver_advanced.AddDofs( mp )

    element_postfix = "2D3N"
    condition_postfix = "2D2N"

    ######## create elements ########
    timer = Timer()
    timer.Start("Create elements")
    # create elements for layer
    layer_name = 'Layer0'
    layer_index = 1
    prop_index = 1
    element_name = params['element_name'] + element_postfix
    prop = mp.Properties[prop_index]
    # prop.SetValue(THERMAL_CONDUCTIVITY, 1.0 )
    # prop.SetValue(THICKNESS,            1.0 )
    # util = ThermalAuxiliaryUtility()
    # util.SetValueForProperties(prop, HEAT_SOURCE, HeatSourceStdProblem1())
    # util.SetValueForProperties(prop, TEMPERATURE_FUNCTION, HeatStdProblem1Solution())
    if params['order'] == 1:
        prop.SetValue(INTEGRATION_ORDER,     1 )
    elif params['order'] == 2:
        prop.SetValue(INTEGRATION_ORDER,     2 )
    # prop.SetValue(LAYER_NAME, layer_name)
    # print(str(mpi.rank) + ": Create element " + element_name + " for layer " + layer_name)
    bulk_elems = mmg_model.AddElements(layer_index, element_name, prop)
    print("Set the element type to " + element_name + " for the layer index " + str(layer_index))
    timer.Stop("Create elements")
    ######## create elements completed ########

    ######## create conditions ########
    if 'condition_name' in params:
        timer.Start("Create conditions")
        for boundary_id in [2, 3, 4, 5]:
            prop_index = boundary_id
            prop = mp.Properties[prop_index]
            prop.SetValue(THICKNESS, 1.0 )
            condition_name = params['condition_name'] + condition_postfix
            mmg_model.AddConditions(boundary_id, condition_name, prop)
            print("Set the condition type to " + condition_name + " for the boundary_index " + str(boundary_id))
        timer.Stop("Create conditions")
    ######## create conditions completed ########

    mmg_model.EndModelPart()

    ######## increase the order ########
    if params['order'] == 2:
        model_part_utils = ModelPartUtilities()
        new_mp = ModelPart("tri6")
        model_part_utils.CopyAndIncreaseOrder(mp, new_mp)
        # print(new_model_part)
    elif params['order'] == 1:
        new_mp = mp

    ######## boundary correction ########

    r_boundary_nodes = {}
    r_boundary_nodes[1] = []
    r_boundary_nodes[2] = []
    r_boundary_nodes[3] = []
    r_boundary_nodes[4] = []
    for cond in new_mp.Conditions:
        if cond.Properties.Id == 2:
            for node in cond.GetNodes():
                r_boundary_nodes[1].append(node.Id)
        if cond.Properties.Id == 3:
            for node in cond.GetNodes():
                r_boundary_nodes[2].append(node.Id)
        if cond.Properties.Id == 4:
            for node in cond.GetNodes():
                r_boundary_nodes[3].append(node.Id)
        if cond.Properties.Id == 5:
            for node in cond.GetNodes():
                r_boundary_nodes[4].append(node.Id)
    r_boundary_nodes[1] = list(set(r_boundary_nodes[1]))
    r_boundary_nodes[2] = list(set(r_boundary_nodes[2]))
    r_boundary_nodes[3] = list(set(r_boundary_nodes[3]))
    r_boundary_nodes[4] = list(set(r_boundary_nodes[4]))

    for boundary_id, brep in params['boundaries'].items():
        for node_id in r_boundary_nodes[boundary_id]:
            node = new_mp.Nodes[node_id]
            node.Set(FREEZE, False)
            brep.ProjectOnSurface(node)
            node.Set(FREEZE, True)

    # # for quadratic edge we shall place the mid-node in the middle of the chord
    # if params['order'] == 2:
    #     Ri = brep1.Radius()
    #     Ro = brep2.Radius()
    #     for cond in new_mp.Conditions:
    #         if (cond.Properties.Id == 2) or (cond.Properties.Id == 3):
    #             if cond.Properties.Id == 2:
    #                 R = Ri
    #             elif cond.Properties.Id == 3:
    #                 R = Ro
    #             n1x = cond.GetNodes()[0].X0
    #             n1y = cond.GetNodes()[0].Y0
    #             n2x = cond.GetNodes()[1].X0
    #             n2y = cond.GetNodes()[1].Y0
    #             a1 = math.atan2(n1y, n1x)
    #             a2 = math.atan2(n2y, n2x)
    #             am = 0.5*(a1 + a2)
    #             nmx = R*math.cos(am)
    #             nmy = R*math.sin(am)
    #             mnode = cond.GetNodes()[2]
    #             mnode.Set(FREEZE, False)
    #             mnode.SetPosition(nmx, nmy, 0.0)
    #             mnode.SetInitialPosition(nmx, nmy, 0.0)
    #             mnode.Set(FREEZE, True)

    #############################################
    ####### assign the layer information ########
    #############################################

    layer_nodes_sets = {}
    layer_sets = {}
    layer_cond_sets = {}
    element_assignments = {}

    ## assign node group
    node_groups = {}
    node_groups['r1'] = r_boundary_nodes[1]
    node_groups['r2'] = r_boundary_nodes[2]
    node_groups['r3'] = r_boundary_nodes[3]
    node_groups['rl'] = r_boundary_nodes[4]

    ## assign node sets
    layer_nodes_sets['boundary'] = []
    layer_nodes_sets['boundary'].extend(r_boundary_nodes[1])
    layer_nodes_sets['boundary'].extend(r_boundary_nodes[2])
    layer_nodes_sets['boundary'].extend(r_boundary_nodes[3])
    layer_nodes_sets['boundary'].extend(r_boundary_nodes[4])

    ## assign elements for layer
    layer_name = 'Layer0'
    layer_sets[layer_name] = []
    for elem in bulk_elems:
        layer_sets[layer_name].append(elem.Id)

    print("ModelPart is constructed successfully")

    return [new_mp, layer_nodes_sets, layer_sets, layer_cond_sets, node_groups, element_assignments]
