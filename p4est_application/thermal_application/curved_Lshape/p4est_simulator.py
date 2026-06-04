import os
import math
import time as time_module
import pprint
import six
from KratosMultiphysics import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.P4estApplication import *

## refinement process that refine everything
class P4RefinementProcessAll():
    def __init__(self, nrefine):
        self.number_of_refinement = nrefine

    def Execute(self, p4est_model):
        nsampling = 5
        for step in range(0, self.number_of_refinement):
            ncell = p4est_model.NumberOfCells()
            refine_vector = ncell*[1]
            balance_option = 0
            p4est_model.Mark(refine_vector)
            p4est_model.Refine()
            p4est_model.Repartition(balance_option)

## refinement process that refine using an indicator vector
class P4RefinementProcessBasedOnRefineVector():
    def __init__(self, refine_vector):
        self.refine_vector = refine_vector
        self.recursive = 0
        self.max_level = 99
        self.balance_option = 0

    def Execute(self, p4est_model):
        # create the refinement vector
        ncell = p4est_model.NumberOfCells()
        print(f"ncell: {ncell}")
        refine_vector = ncell*[0]
        for cell_id in self.refine_vector:
            refine_vector[cell_id - 1] = 1
        #
        p4est_model.Mark(refine_vector)
        p4est_model.Refine(self.recursive, self.max_level)
        p4est_model.Repartition(self.balance_option)

## set the layer index for the elements in a model
def SetLayerIndex(model):
    for elem in model.model_part.Elements:
        elem.SetValue(LAYER_INDEX, 1)

## set the boundary indices for the elements in a model
def SetBoundaryFlags(model):
    p4est_util = P4estUtilities()
    for cond_id, elem_id in model.element_assignments.items():
        cond = model.model_part.Conditions[cond_id]
        elem = model.model_part.Elements[elem_id]
        fid = p4est_util.FindFaceIndex(cond, elem)
        # print(f"fid: {fid}, prop id: {cond.Properties.Id}")
        # print(" cond nodes:")
        # for node in cond.GetNodes():
        #     print(f"  {node.Id}")
        # print(" elem nodes:")
        # for node in elem.GetNodes():
        #     print(f"  {node.Id}")
        if elem.Has(BOUNDARY_INDICES):
            bindices = elem.GetValue(BOUNDARY_INDICES)
            bindices[fid] = cond.Properties.Id
            elem.SetValue(BOUNDARY_INDICES, bindices)
        else:
            bindices = [-1]*4
            bindices[fid] = cond.Properties.Id
            elem.SetValue(BOUNDARY_INDICES, bindices)

    for elem in model.model_part.Elements:
        if elem.Has(BOUNDARY_INDICES):
            print(elem.GetValue(BOUNDARY_INDICES))

## create the P4est model from the model_part
def CreateP4estModel(model_part, params):
    order = params['order']

    # specify the order of the integration
    if order == 1:
        p4est_order = P4estFirstOrder()
    elif order == 2:
        p4est_order = P4estSecondOrder()

    # create a sample quad for the p4est tree
    p4est_quad_nodal = P4estQuadData(p4est_order.Value)
    p4est_quad_nodal.Initialize()
    p4est_quad_nodal.Register(DISPLACEMENT)
    p4est_quad_nodal.Finalize()

    p4est_quad_int = P4estQuadData(p4est_order.Value)
    p4est_quad_int.Initialize()
    p4est_quad_int.Finalize()

    p4est_quad = P4estQuad(p4est_quad_nodal, p4est_quad_int)

    p4est_util = P4estUtilities()
    p4est_model = p4est_util.CreateModel(p4est_order, mpi.world, p4est_quad)
    p4est_model.Initialize(model_part)
    p4est_model.Repartition(0)
    mpi.world.barrier()

    # ## test seg fault
    # balance_option = 0
    # p4est_model.Mark([1], P4estRefineFlag.TO_REFINE)
    # p4est_model.Refine()
    # p4est_model.Repartition(balance_option)
    # ####

    # register the boundary
    if 'boundaries' in params:
        if params['boundaries'] != None:
            for boundary_index, brep in six.iteritems(params['boundaries']):
                # p4est_model.RegisterBoundaryAndFlagFaces(boundary_index, brep)
                p4est_model.RegisterBoundary(boundary_index, brep)

    # refinement
    if 'initial_refinement_process' in params:
        if params['initial_refinement_process'] != None:
            params['initial_refinement_process'].Execute(p4est_model)

    return p4est_model

## create the model_part out from p4est_model
def ConstructSystemModelPart(p4est_model, params):
    # start to construct the Kratos model_part
    p4est_model.BeginModelPart()

    mp = p4est_model.GetModelPart()
    mp.SetBufferSize(2)
    import structural_solver_advanced
    mp.AddNodalSolutionStepVariable(TEMPERATURE)
    mp.AddNodalSolutionStepVariable(LAGRANGE_TEMPERATURE)
    mp.AddNodalSolutionStepVariable(TEMPERATURE_ERROR)
    mp.AddNodalSolutionStepVariable(REFERENCE_TEMPERATURE)
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX)

    # create nodes
    hanging_nodes = p4est_model.CreateNodes()
    for node in mp.Nodes:
        node.AddDof(TEMPERATURE)

    # print("Hanging nodes:")
    # for node in hanging_nodes:
    #     print(node.Id)

    # create constraints
    timer = Timer()
    timer.Start("Create constraints")

    last_constraint_id = 0
    constraint_name = "LinearMasterSlaveConstraint"
    last_constraint_id = p4est_model.CreateConstraints(hanging_nodes, last_constraint_id, constraint_name, TEMPERATURE)

    timer.Stop("Create constraints")

    if params['order'] == 1:
        element_postfix = "2D4N"
        condition_postfix = "2D2N"
    elif params['order'] == 2:
        element_postfix = "2D9N"
        condition_postfix = "2D3N"

    ######## create elements ########
    timer.Start("Create elements")
    # create elements for layer
    layer_name = 'Layer0'
    layer_index = 1
    prop_index = 1
    element_name = params['element_name'] + element_postfix
    prop = mp.Properties[prop_index]
    if params['order'] == 1:
        prop.SetValue(INTEGRATION_ORDER,     2 )
    elif params['order'] == 2:
        prop.SetValue(INTEGRATION_ORDER,     3 )
    # prop.SetValue(LAYER_NAME, layer_name)
    # print(str(mpi.rank) + ": Create element " + element_name + " for layer " + layer_name)
    bulk_elems = p4est_model.AddElements(layer_index, element_name, prop)
    print("Set the element type to " + element_name + " for the layer index " + str(layer_index))
    timer.Stop("Create elements")
    ######## create elements completed ########

    ######## create conditions ########
    if 'condition_name' in params:
        timer.Start("Create conditions")
        for boundary_id in [1, 2, 3, 4]:
            prop_index = boundary_id
            prop = mp.Properties[prop_index]
            prop.SetValue(THICKNESS, 1.0 )
            condition_name = params['condition_name'] + condition_postfix
            p4est_model.AddConditions(boundary_id, condition_name, prop)
            print("Set the condition type to " + condition_name + " for the boundary_index " + str(boundary_id))
        timer.Stop("Create conditions")
    ######## create conditions completed ########

    p4est_model.EndModelPart()
    # sys.exit(0)

    for elem in mp.Elements:
        print(f"elem {elem.Id} nodes:")
        for node in elem.GetNodes():
            print(f" {node.Id}")

    ## obtain and mark the boundary nodes
    r_boundary_nodes = {}

    for i in [1, 2, 3, 4]:
        r_boundary_nodes[i] = []
        boundary_nodes = p4est_model.GetBoundaryNodes(i)
        for node in boundary_nodes:
            node.Set(BOUNDARY, True)
            r_boundary_nodes[i].append(node.Id)

    #### mesh smoothing
    if params['mesh_smoothing'] == "laplacian":
        find_nodal_neighbour_process = FindNodalNeighboursProcess(mp, 10, 10)
        find_nodal_neighbour_process.Execute()
        for i in range(0, 1):
            mesh_smoothing_process = LaplacianMeshSmoothingProcess(mp)
            mesh_smoothing_process.Execute()

    #############################################
    ####### assign the layer information ########
    #############################################

    layer_nodes_sets = {}
    layer_sets = {}
    cell_to_element = {}
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
    # print(f"r_boundary_nodes[1]: {r_boundary_nodes[1]}")
    # print(f"r_boundary_nodes[2]: {r_boundary_nodes[2]}")
    # print(f"r_boundary_nodes[3]: {r_boundary_nodes[3]}")
    # print(f"r_boundary_nodes[4]: {r_boundary_nodes[4]}")
    # sys.exit(0)

    ## assign elements for layer
    layer_name = 'Layer0'
    layer_sets[layer_name] = []
    for elem in bulk_elems:
        layer_sets[layer_name].append(elem.Id)
        cell_to_element[elem.GetValue(PARENT_CELL_ID)] = elem.Id

    print("ModelPart is constructed successfully")

    return [mp, layer_nodes_sets, layer_sets, layer_cond_sets, node_groups, element_assignments]
