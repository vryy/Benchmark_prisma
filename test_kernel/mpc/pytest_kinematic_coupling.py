from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.LayerApplication import *
# from KratosMultiphysics.mpi import *
# from KratosMultiphysics.P4estApplication import *
# from inspect import currentframe, getframeinfo

def main(output=True, enable_constraint=True):
    # Define a Model
    mp = ModelPart("Main")

    mp.AddNodalSolutionStepVariable(DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(REACTION)
    mp.AddNodalSolutionStepVariable(ROTATION)
    mp.AddNodalSolutionStepVariable(TORQUE)
    mp.AddNodalSolutionStepVariable(ELASTIC_BEDDING_STIFFNESS) # for nodal springs

    # Create Nodes
    mp.CreateNewNode(1, 0.00,0.00,0.00)
    mp.CreateNewNode(2, 0.50,0.00,0.00)
    mp.CreateNewNode(3, 1.00,0.00,0.00)
    mp.CreateNewNode(4, 0.00,0.50,0.00)
    mp.CreateNewNode(5, 0.50,0.50,0.00)
    mp.CreateNewNode(6, 1.00,0.50,0.00)
    mp.CreateNewNode(7, 0.00,1.00,0.00)
    mp.CreateNewNode(8, 0.50,1.00,0.00)
    mp.CreateNewNode(9, 1.00,1.00,0.00)

    # Add Dofs
    for node in mp.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
        node.AddDof(ROTATION_Z, TORQUE_Z)

    # define properties
    mp.GetProperties()[1].SetValue(YOUNG_MODULUS, 2)
    mp.GetProperties()[1].SetValue(POISSON_RATIO, 0.3)
    mp.GetProperties()[1].SetValue(THICKNESS, 1.0)
    mp.GetProperties()[1].SetValue(DENSITY, 1.0)
    mp.GetProperties()[1].SetValue(CONSTITUTIVE_LAW, PlaneStress())
    mp.GetProperties()[1].SetValue(BODY_FORCE, ZeroVector(3))

    springs_stiffness = Vector(3)
    springs_stiffness[0] = 1.0
    springs_stiffness[1] = 1.0
    springs_stiffness[2] = 1.0

    # create Elements
    element_name = "KinematicLinear2D4N"
    mp.CreateNewElement(element_name, 1, [1,2,5,4], mp.GetProperties()[1])
    mp.CreateNewElement(element_name, 2, [2,3,6,5], mp.GetProperties()[1])
    mp.CreateNewElement(element_name, 3, [4,5,8,7], mp.GetProperties()[1])
    mp.CreateNewElement(element_name, 4, [5,6,9,8], mp.GetProperties()[1])

    # create Conditions
    mp.CreateNewCondition("ElasticPointConstraint", 1, [7], mp.GetProperties()[2]) # nodal springs to keep the node not aligned
    mp.Nodes[7].SetSolutionStepValue(ELASTIC_BEDDING_STIFFNESS, springs_stiffness)

    # for elem in mp.Elements:
    #     elem.SetValue(LAYER_INDEX, 1)

    # Constraints
    # to keep the surface flat and behave like rigid body
    if enable_constraint:
        constraint_name = "LinearMasterSlaveConstraint"
        slave_vars = [DISPLACEMENT_X, DISPLACEMENT_Y]
        master_vars = [DISPLACEMENT_X, DISPLACEMENT_Y, ROTATION_Z]
        slave_node_ids = [7, 9]
        master_node_id = 8

        constraint_id = 1
        for slave_node_id in slave_node_ids:
            relation_matrix = Matrix(2, 3)
            relation_matrix[0,0] = 1.0
            relation_matrix[0,1] = 0.0
            relation_matrix[0,2] = -(mp.Nodes[slave_node_id].Y0 - mp.Nodes[master_node_id].Y0)
            relation_matrix[1,0] = 0.0
            relation_matrix[1,1] = 1.0
            relation_matrix[1,2] = mp.Nodes[slave_node_id].X0 - mp.Nodes[master_node_id].X0
            constant_vector = ZeroVector(2)
            mp.CreateNewLinearMasterSlaveConstraint(constraint_name, constraint_id, mp.Nodes[slave_node_id], mp.Nodes[master_node_id], slave_vars, master_vars, relation_matrix, constant_vector)
            constraint_id += 1

    # # impose BC 1
    # mp.Nodes[1].Fix(DISPLACEMENT_Y)
    # mp.Nodes[2].Fix(DISPLACEMENT_Y)
    # mp.Nodes[3].Fix(DISPLACEMENT_Y)

    # mp.Nodes[1].Fix(DISPLACEMENT_X)
    # mp.Nodes[4].Fix(DISPLACEMENT_X)
    # mp.Nodes[7].Fix(DISPLACEMENT_X)

    # mp.Nodes[8].Fix(DISPLACEMENT_X)
    # mp.Nodes[8].Fix(DISPLACEMENT_Y)

    #

    # impose BC 2
    mp.Nodes[1].Fix(DISPLACEMENT_X)
    mp.Nodes[1].Fix(DISPLACEMENT_Y)
    mp.Nodes[2].Fix(DISPLACEMENT_X)
    mp.Nodes[2].Fix(DISPLACEMENT_Y)
    mp.Nodes[3].Fix(DISPLACEMENT_X)
    mp.Nodes[3].Fix(DISPLACEMENT_Y)

    mp.Nodes[8].Fix(DISPLACEMENT_X)
    mp.Nodes[8].Fix(DISPLACEMENT_Y)

    for node in mp.Nodes:
        node.Fix(ROTATION_X)
        node.Fix(ROTATION_Y)
        node.Fix(DISPLACEMENT_Z)

    if not enable_constraint:
        for node in mp.Nodes:
            node.Fix(ROTATION_Z)

    linear_solver = SkylineLUFactorizationSolver()
    builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraints(linear_solver)
    # builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(linear_solver)
    scheme = ResidualBasedIncrementalUpdateStaticScheme()
    convergence_criterion = ResidualCriteria(1e-10, 1e-12)

    max_iters = 100
    compute_reactions = True
    reform_step_dofs = False
    move_mesh_flag = False
    strategy = ResidualBasedNewtonRaphsonStrategy(mp, scheme, linear_solver, convergence_criterion, builder_and_solver, max_iters, compute_reactions,reform_step_dofs, move_mesh_flag)
    # strategy = ResidualBasedLinearStrategy(mp, scheme, linear_solver, builder_and_solver,max_iters,compute_reactions,reform_step_dofs,move_mesh_flag)
    strategy.SetEchoLevel(2)
    strategy.Initialize()
    mp.Nodes[8].SetSolutionStepValue(DISPLACEMENT_X,0.1)
    strategy.Check()
    strategy.Solve()

    print("Analysis completed")

    if output:
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        gid_io = SDGidPostIO( "plate", post_mode, multi_file_flag, write_deformed_flag, write_elements )

        time = 1.0
        gid_io.InitializeMesh( time )
        mesh = mp.GetMesh()
        gid_io.WriteMesh( mesh )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults( time, mp.GetMesh() )
        gid_io.WriteNodalResults(DISPLACEMENT, time, 0)
        gid_io.WriteNodalResults(ROTATION, time, 0)
        gid_io.FinalizeResults()
        gid_io.Reset()

    return mp

def test():
    mp = main(output=False, enable_constraint=True)

    tol = 1e-16
    u7x = mp.Nodes[7].GetSolutionStepValue(DISPLACEMENT_X)
    u7y = mp.Nodes[7].GetSolutionStepValue(DISPLACEMENT_Y)
    u8x = mp.Nodes[8].GetSolutionStepValue(DISPLACEMENT_X)
    u8y = mp.Nodes[8].GetSolutionStepValue(DISPLACEMENT_Y)
    u9x = mp.Nodes[9].GetSolutionStepValue(DISPLACEMENT_X)
    u9y = mp.Nodes[9].GetSolutionStepValue(DISPLACEMENT_Y)
    print("u7x: %.16e" % (u7x))
    print("u7y: %.16e" % (u7y))
    print("u8x: %.16e" % (u8x))
    print("u8y: %.16e" % (u8y))
    print("u9x: %.16e" % (u9x))
    print("u9y: %.16e" % (u9y))
    assert(abs(u7x - 1e-1) < tol)
    assert(abs(u7y - 2.6011381912656194e-02) < tol)
    assert(abs(u8x - 1e-1) < tol)
    assert(abs(u8y) < tol)
    assert(abs(u9x - 1e-1) < tol)
    assert(abs(u9y + 2.6011381912656194e-02) < tol)
    print("Test passed")

def tag():
    return "MPC"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=True, enable_constraint=True)
