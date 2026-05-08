import os
import math
import time as time_module
import pprint
import six
from KratosMultiphysics import *
# from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.FiniteCellApplication import *
from KratosMultiphysics.FiniteCellStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
# from KratosMultiphysics.MortarApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.P4estApplication import *

## refinemenet process along the boundary
class P4RefinementProcessAlongBoundary():
    def __init__(self, brep, nrefine):
        self.brep = brep
        self.number_of_refinement = nrefine

    def Execute(self, p4est_model):
        nsampling = 5
        for step in range(0, self.number_of_refinement):
            all_stats = p4est_model.ProbeCutStatus(self.brep, nsampling)
            refine_vector = []
            for stat in all_stats:
                if stat == BRep._CUT:
                    refine_vector.append(1)
                else:
                    refine_vector.append(0)
            balance_option = 0
            p4est_model.Mark(refine_vector)
            p4est_model.Refine()
            p4est_model.Repartition(balance_option)

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

## set the layer index for the elements in a model
def SetLayerIndex(model):
    for elem in model.model_part.Elements:
        elem.SetValue(LAYER_INDEX, 1)

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
    p4est_quad_int.Register(STRESSES)
    p4est_quad_int.Register(PRESTRESS)
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
                p4est_model.RegisterBoundaryAndFlagFaces(boundary_index, brep)

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
    structural_solver_advanced.AddVariables( mp )
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX)

    # create nodes
    hanging_nodes = p4est_model.CreateNodes()
    import structural_solver_advanced
    structural_solver_advanced.AddDofs( mp )

    # create constraints
    timer = Timer()
    timer.Start("Create constraints")

    last_constraint_id = 1
    constraint_name = "LinearMasterSlaveConstraint"
    last_constraint_id = p4est_model.CreateConstraints(hanging_nodes, last_constraint_id, constraint_name, DISPLACEMENT_X)
    last_constraint_id = p4est_model.CreateConstraints(hanging_nodes, last_constraint_id, constraint_name, DISPLACEMENT_Y)

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
    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0
    gravity[2] = 0.0
    prop = mp.Properties[prop_index]
    prop.SetValue(DENSITY,            0 )
    prop.SetValue(BODY_FORCE,            ZeroVector(3) )
    prop.SetValue(YOUNG_MODULUS,      2.1e+11 )
    prop.SetValue(POISSON_RATIO,          0.3 )
    prop.SetValue(THICKNESS, 1.0 )
    prop.SetValue(TENSILE_STRENGTH,          0.27712812921102037e9 )
    if params['order'] == 1:
        prop.SetValue(INTEGRATION_ORDER,     2 )
    elif params['order'] == 2:
        prop.SetValue(INTEGRATION_ORDER,     3 )
    # prop.SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
    prop.SetValue(CONSTITUTIVE_LAW, ValuesContainerConstitutiveLaw(PlaneStrain()) ) # select this in order to visualize the L^2 and H^1 error
    print("Linear elastic model selected for Layer0, description: Steel")
    # prop.SetValue(LAYER_NAME, layer_name)
    # print(str(mpi.rank) + ": Create element " + element_name + " for layer " + layer_name)
    bulk_elems = p4est_model.AddElements(layer_index, element_name, prop)
    print("Set the element type to " + element_name + " for the layer index " + str(layer_index))
    timer.Stop("Create elements")
    ######## create elements completed ########

    ######## create conditions ########
    if 'condition_name' in params:
        timer.Start("Create conditions")
        boundary_id = 1
        prop_index = 2
        prop = mp.Properties[prop_index]
        prop.SetValue(THICKNESS, 1.0 )
        condition_name = params['condition_name'] + condition_postfix
        p4est_model.AddConditions(boundary_id, condition_name, prop)
        print("Set the condition type to " + condition_name + " for the boundary_index " + str(boundary_id))
        timer.Stop("Create conditions")
    ######## create conditions completed ########

    p4est_model.EndModelPart()
    # sys.exit(0)

    #############################################
    ####### assign the layer information ########
    #############################################

    layer_sets = {}
    cell_to_element = {}
    layer_cond_sets = {}
    element_assignments = {}
    node_groups = {}

    ## assign elements for layer
    layer_name = 'Layer0'
    layer_sets[layer_name] = []
    for elem in bulk_elems:
        layer_sets[layer_name].append(elem.Id)
        cell_to_element[elem.GetValue(PARENT_CELL_ID)] = elem.Id

    print("ModelPart is constructed successfully")

    return [mp, layer_sets, layer_cond_sets, node_groups, element_assignments]
