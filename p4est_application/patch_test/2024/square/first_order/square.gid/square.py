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
        element_name = "KinematicLinear2D4N"
    elif order == 2:
        element_name = "KinematicLinear2D9N"

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

def ApplyBC(mp):
    tol = 1.0e-6
    prescribed_nodes = []
    for node in mp.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)
    return prescribed_nodes

def Solve(model, p4est_model, time, output=True):
    prescribed_nodes = ApplyBC(model.model_part)

    time = time + 100.0
    model.Solve(time, 0, 0, 0, 0)

    time = time + 1.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.1)

    model.Solve(time, 0, 0, 0, 0)
    p4est_model.SynchronizeBackwardUsingExtrapolation()

    if output:
        model.WriteOutput(time)
        p4est_model.ExportVTK("square", time, True, True)

    for node in model.model_part.Nodes:
        print("displacement at node " + str(node.Id) + ": " + str(node.GetSolutionStepValue(DISPLACEMENT_X)) + " " + str(node.GetSolutionStepValue(DISPLACEMENT_Y)))

    for elem in model.model_part.Elements:
        stresses = elem.GetValuesOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
        stress_xx = []
        for stress in stresses:
            stress_xx.append(round(stress[0], 2))
        print("stress_xx at element " + str(elem.Id) + ": " + str(stress_xx))

    return time

def main(logging=True, output=True):

    if mpi.size == 1:
        model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    else:
        model = square_parallel_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
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
    p4est_model.Repartition(0)
    mpi.world.barrier()

    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    ##### simulation #####

    time = 0.0
    time = Solve(model, p4est_model, time, output=output)

    ###############################

    # mark and refine
    p4est_model.Mark([1], P4estRefineFlag.TO_REFINE)
    p4est_model.Refine()
    p4est_model.Coarsen()
    p4est_model.Repartition(0)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # solve
    time = Solve(model, p4est_model, time, output=output)

    ###############################

    # mark and refine
    p4est_model.Mark([1,4], P4estRefineFlag.TO_REFINE)
    p4est_model.Refine()
    p4est_model.Coarsen()
    p4est_model.Repartition(0)
    # p4est_model.ExportVTK("square", time, True, True)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # solve
    time = Solve(model, p4est_model, time, output=output)

    ###############################

    # mark and refine
    p4est_model.Mark([2,10], P4estRefineFlag.TO_REFINE)
    p4est_model.Refine()
    p4est_model.Coarsen()
    p4est_model.Repartition(0)
    # p4est_model.ExportVTK("square", time, True, True)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # solve
    time = Solve(model, p4est_model, time, output=output)

    ###############################

    # mark and coarsen
    p4est_model.Mark([3,11], P4estRefineFlag.TO_COARSEN)
    p4est_model.Refine()
    p4est_model.Coarsen()
    p4est_model.Repartition(0)
    # p4est_model.ExportVTK("square", time, True, True)

    # create the model_part out from p4est
    mp = ConstructModelPart(p4est_model, p4est_order.Value)
    model.SetModelPart(mp)
    model.InitializeModel()

    # synchonize the data from forest to the model_part
    p4est_model.SynchronizeForward()

    # solve
    time = Solve(model, p4est_model, time, output=output)

    return model

def test():
    model = main(logging = False, output = False)

    for elem in model.model_part.Elements:
        stresses = elem.GetValuesOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
        stress_xx = []
        for stress in stresses:
            stress_xx.append(round(stress[0], 2))
        for stress in stress_xx:
            assert(abs(stress - 0.2) < 1e-10)
    print("Test passed")

def tag():
    return "p4est"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
