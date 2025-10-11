import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import square_include
from square_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def ConstructModelPart(p4est_model, order):
    # start to construct the Kratos model_part
    p4est_model.BeginModelPart()

    mp = p4est_model.GetModelPart()
    mp.SetBufferSize(2)
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( mp )
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX)

    prop = mp.Properties[1]
    prop.SetValue(CONSTITUTIVE_LAW, DummyConstitutiveLaw())
    prop.SetValue(BODY_FORCE, ZeroVector(3))
    if order == 1:
        element_name = "KinematicLinearFiniteCell2D4N"
    elif order == 2:
        element_name = "KinematicLinearFiniteCell2D9N"

    constraint_name = "LinearMasterSlaveConstraint"

    hanging_nodes = p4est_model.CreateNodes()
    import structural_solver_advanced
    structural_solver_advanced.AddDofs( mp )

    start_adding_time = time_module.time()
    constraint_id = 1
    # for hnode in hanging_nodes:
    #     print(hnode)
    #     p4est_model.CreateConstraint(hnode, constraint_id, constraint_name, DISPLACEMENT_X)
    #     constraint_id = constraint_id + 1
    #     p4est_model.CreateConstraint(hnode, constraint_id, constraint_name, DISPLACEMENT_Y)
    #     constraint_id = constraint_id + 1
    constraint_id = p4est_model.CreateConstraints(hanging_nodes, constraint_id, constraint_name, DISPLACEMENT_X)
    constraint_id = p4est_model.CreateConstraints(hanging_nodes, constraint_id, constraint_name, DISPLACEMENT_Y)
    end_adding_time = time_module.time()
    print("Time to add constraints: " + str(end_adding_time - start_adding_time) + " s")

    layer_idx = 1
    p4est_model.AddElements(layer_idx, element_name, prop)
    p4est_model.EndModelPart()

    return mp

def main(logging=True, output=True):

    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    # set the layer index. It is important to set the layer index.
    # The elements will be created based on layer
    for elem in model.model_part.Elements:
        elem.SetValue(LAYER_INDEX, 1)

    # specify the order of the integration
    p4est_order = P4estFirstOrder()

    # create a sample quad for the p4est tree
    p4est_quad_nodal = P4estQuadData(p4est_order.Value)
    p4est_quad_nodal.Initialize()
    p4est_quad_nodal.Register(DISPLACEMENT)
    p4est_quad_nodal.Finalize()

    p4est_quad_int = P4estQuadData(p4est_order.Value)
    p4est_quad_int.Initialize()
    p4est_quad_int.Register(STRESSES)
    p4est_quad_int.Finalize()

    p4est_quad = P4estQuad(p4est_quad_nodal, p4est_quad_int)

    # create the P4est model
    p4est_util = P4estUtilities()
    p4est_model = p4est_util.CreateModel(p4est_order, mpi.world, p4est_quad)
    p4est_model.Initialize(model.model_part)

    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # mark and refine
    p4est_model.Mark([1], P4estRefineFlag.TO_REFINE)
    p4est_model.RefineCoarsen()
    p4est_model.Repartition(0)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    time = 1.1
    if output:
        model.WriteOutput(time)

    ###############################

    # mark and refine
    p4est_model.Mark([1, 2, 4], P4estRefineFlag.TO_REFINE)
    p4est_model.RefineCoarsen()
    p4est_model.Repartition(0)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # export the results without solve
    time = 1.2
    if output:
        model.WriteOutput(time)

    # mark and refine
    p4est_model.Mark([1, 2, 7, 11, 13], P4estRefineFlag.TO_REFINE)
    p4est_model.RefineCoarsen()
    p4est_model.Repartition(0)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # export the results without solve
    time = 1.3
    if output:
        model.WriteOutput(time)

    # mark and coarsen
    p4est_model.Mark([15], P4estRefineFlag.TO_COARSEN)
    p4est_model.RefineCoarsen()
    p4est_model.Repartition(0)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()
    model.solver.solver.Initialize()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # export the results without solve
    time = 1.4
    if output:
        model.WriteOutput(time)

    ########################FINITE CELL SIMULATION#################

    import material_properties_utility
    import simulator

    #import quadrature_data
    matfile_path = os.getcwd() + '/matfile.dat'

    params = {}
    params["quadrature_method"] = "quadtree"
    params["number_of_samplings"] = 4
    params["qt_depth"] = 4
    params["cut_cell_quadrature_method"] = 0x03
    params["small_weight"] = 0.0
    params["material_properties_utility"] = material_properties_utility.MaterialPropertiesUtility(matfile_path)
    params["material_properties_utility"].search_type = "by_name"
    params["material_properties_utility"].mat_type = "elastic"
    params["write_output_per_each_step"] = output
    params["export_physical_integration_point"] = True
    params["physical_integration_point_prop_id"] = 88
    params["export_quadtree_cell"] = False
    params["sample_quad_element_name"] = "DummySurfaceElement2D4N"
    params["write_quadrature_to_file"] = False
    params["quad_filename"] = "quadrature_data.py"
    params["quad_filetype"] = "python"
    params["quad_accuracy"] = 20
    params["enable_ghost_penalty"] = False
    params["space_dim"] = 2
    params["estimated_number_of_neighbours"] = 10
    params["sample_ghost_penalty_condition"] = GhostPenaltyDisplacementGradientCondition()
    params["ghost_penalty_properties"] = model.model_part.Properties[4]
    model.model_part.Properties[4].SetValue(GHOST_PENALTY_STABILIZATION_FACTOR, 1.0e-1)
    model.model_part.Properties[4].SetValue(GHOST_PENALTY_STABILIZATION_ORDER, 2)

    sim = simulator.Simulator(params)
    sim.Run(model)

    return model

def test():
    model = main(logging = False, output = False)

    assert(len(model.model_part.Elements) == 40)

    ref_disp = {}
    ref_disp[35] = [4.0563235786e-02, -2.8104024007e-02]
    ref_disp[36] = [1.7785230712e-02, -4.3482942261e-02]
    ref_disp[40] = [0.0000000000e+00, -5.0421396843e-02]
    ref_disp[41] = [1.0000000000e-01, -1.1762176130e-02]

    for node in model.model_part.Nodes:
        if node.Id in ref_disp:
            ux = node.GetSolutionStepValue(DISPLACEMENT_X)
            uy = node.GetSolutionStepValue(DISPLACEMENT_Y)
            assert(abs(ux - ref_disp[node.Id][0]) < 1e-9)
            assert(abs(uy - ref_disp[node.Id][1]) < 1e-9)

    print("Test passed")

def tag():
    return "p4est,finite-cell"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
