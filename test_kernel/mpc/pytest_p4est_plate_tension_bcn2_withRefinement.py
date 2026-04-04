from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.LayerApplication import *
# from KratosMultiphysics.mpi import *
# from KratosMultiphysics.P4estApplication import *
# from inspect import currentframe, getframeinfo

def main(output=True):
    # Define a Model
    mp = ModelPart("Main")

    mp.AddNodalSolutionStepVariable(DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(REACTION)
    mp.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)

    # Create Nodes
    mp.CreateNewNode(1, 0.50,0.00,0.00)
    mp.CreateNewNode(2, 0.50,0.50,0.00)
    mp.CreateNewNode(3, 0.00,0.00,0.00)
    mp.CreateNewNode(4, 0.00,0.50,0.00)
    mp.CreateNewNode(5, 1.00,0.00,0.00)
    mp.CreateNewNode(6, 1.00,0.50,0.00)
    mp.CreateNewNode(7, 0.50,1.00,0.00)
    mp.CreateNewNode(8, 0.00,1.00,0.00)
    mp.CreateNewNode(9, 1.00,0.75,0.00)
    mp.CreateNewNode(10, 0.75,0.75,0.00)
    mp.CreateNewNode(11, 1.00,1.00,0.00)
    mp.CreateNewNode(12, 0.75,1.00,0.00)
    mp.CreateNewNode(13, 0.75,0.50,0.00)
    mp.CreateNewNode(14, 0.50,0.75,0.00)

    # Add Dofs
    for node in mp.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)

    # define properties
    mp.GetProperties()[1].SetValue(YOUNG_MODULUS, 2)
    mp.GetProperties()[1].SetValue(POISSON_RATIO, 0.3)
    mp.GetProperties()[1].SetValue(THICKNESS, 1.0)
    mp.GetProperties()[1].SetValue(DENSITY, 1.0)
    mp.GetProperties()[1].SetValue(CONSTITUTIVE_LAW, PlaneStress())
    mp.GetProperties()[1].SetValue(BODY_FORCE, ZeroVector(3))

    # create Elements
    mp.CreateNewElement("KinematicLinear2D4N", 1, [1,2,4,3], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 2, [1,5,6,2], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 3, [4,2,7,8], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 4, [6,9,10,13], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 5, [9,11,12,10], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 6, [10,12,7,14], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 7, [13,10,14,2], mp.GetProperties()[1])

    # for elem in mp.Elements:
    #     elem.SetValue(LAYER_INDEX, 1)

    #Constraints
    relation_matrix = Matrix(1, 2)
    relation_matrix[0, 0] = 0.5
    relation_matrix[0, 1] = 0.5
    constant_vector = ZeroVector(1)
    mp.CreateNewLinearMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, mp.Nodes[14], [mp.Nodes[2], mp.Nodes[7]], DISPLACEMENT_X, relation_matrix, constant_vector)
    mp.CreateNewLinearMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, mp.Nodes[14], [mp.Nodes[2], mp.Nodes[7]], DISPLACEMENT_Y, relation_matrix, constant_vector)

    mp.CreateNewLinearMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, mp.Nodes[13], [mp.Nodes[2], mp.Nodes[6]], DISPLACEMENT_X, relation_matrix, constant_vector)
    mp.CreateNewLinearMasterSlaveConstraint("LinearMasterSlaveConstraint", 4, mp.Nodes[13], [mp.Nodes[2], mp.Nodes[6]], DISPLACEMENT_Y, relation_matrix, constant_vector)

    # impose BC
    mp.Nodes[3].Fix(DISPLACEMENT_X)
    mp.Nodes[3].Fix(DISPLACEMENT_Y)
    mp.Nodes[4].Fix(DISPLACEMENT_X)
    mp.Nodes[8].Fix(DISPLACEMENT_X)
    mp.Nodes[1].Fix(DISPLACEMENT_Y)
    mp.Nodes[5].Fix(DISPLACEMENT_Y)

    mp.Nodes[5].Fix(DISPLACEMENT_X)
    mp.Nodes[6].Fix(DISPLACEMENT_X)
    mp.Nodes[9].Fix(DISPLACEMENT_X)
    mp.Nodes[11].Fix(DISPLACEMENT_X)

    # create a sample quad for the p4est tree
    # num_integration_points = 1
    # p4est_quad = P4estQuad(num_integration_points)
    # p4est_quad.Initialize()
    # p4est_quad.Register(DISPLACEMENT)
    # p4est_quad.Register(STRESSES)
    # p4est_quad.Finalize()
    # # create the P4est model
    # p4est_model = P4estModel(mpi.world, p4est_quad)
    # # p4est_model.Initialize(mp)


    linear_solver = SkylineLUFactorizationSolver()
    builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraints(linear_solver)
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
    mp.Nodes[5].SetSolutionStepValue(DISPLACEMENT_X,0.1)
    mp.Nodes[6].SetSolutionStepValue(DISPLACEMENT_X,0.1)
    mp.Nodes[9].SetSolutionStepValue(DISPLACEMENT_X,0.1)
    mp.Nodes[11].SetSolutionStepValue(DISPLACEMENT_X,0.1)
    strategy.Check()
    strategy.Solve()

    #Results
    print("ux Node 1 =" , mp.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 2 =" , mp.Nodes[2].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 3 =" , mp.Nodes[3].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 4 =" , mp.Nodes[4].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 5 =" , mp.Nodes[5].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 6 =" , mp.Nodes[6].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 7 =" , mp.Nodes[7].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 8 =" , mp.Nodes[8].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 9 =" , mp.Nodes[9].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 10 =" , mp.Nodes[10].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 11 =" , mp.Nodes[11].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 12 =" , mp.Nodes[12].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 13 =" , mp.Nodes[13].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 14 =" , mp.Nodes[14].GetSolutionStepValue(DISPLACEMENT_X))
    print("____________________________________________________" )
    print("uy Node 1 =" , mp.Nodes[1].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 2 =" , mp.Nodes[2].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 3 =" , mp.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 4 =" , mp.Nodes[4].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 5 =" , mp.Nodes[5].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 6 =" , mp.Nodes[6].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 7 =" , mp.Nodes[7].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 8 =" , mp.Nodes[8].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 9 =" , mp.Nodes[9].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 10 =" , mp.Nodes[10].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 11 =" , mp.Nodes[11].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 12 =" , mp.Nodes[12].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 13 =" , mp.Nodes[13].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 14 =" , mp.Nodes[14].GetSolutionStepValue(DISPLACEMENT_Y))
    print("____________________________________________________" )
    print("rx Node 1" , mp.Nodes[1].GetSolutionStepValue(REACTION_X))
    print("rx Node 2" , mp.Nodes[2].GetSolutionStepValue(REACTION_X))
    print("rx Node 3" , mp.Nodes[3].GetSolutionStepValue(REACTION_X))
    print("rx Node 4" , mp.Nodes[4].GetSolutionStepValue(REACTION_X))
    print("rx Node 5" , mp.Nodes[5].GetSolutionStepValue(REACTION_X))
    print("rx Node 6" , mp.Nodes[6].GetSolutionStepValue(REACTION_X))
    print("rx Node 7" , mp.Nodes[7].GetSolutionStepValue(REACTION_X))
    print("rx Node 8" , mp.Nodes[8].GetSolutionStepValue(REACTION_X))
    print("rx Node 9" , mp.Nodes[9].GetSolutionStepValue(REACTION_X))
    print("rx Node 10" , mp.Nodes[10].GetSolutionStepValue(REACTION_X))
    print("rx Node 11" , mp.Nodes[11].GetSolutionStepValue(REACTION_X))
    print("rx Node 12" , mp.Nodes[12].GetSolutionStepValue(REACTION_X))
    print("rx Node 13" , mp.Nodes[13].GetSolutionStepValue(REACTION_X))
    print("rx Node 14" , mp.Nodes[14].GetSolutionStepValue(REACTION_X))

    print("____________________________________________________" )
    print("ry Node 1" , mp.Nodes[1].GetSolutionStepValue(REACTION_Y))
    print("ry Node 2" , mp.Nodes[2].GetSolutionStepValue(REACTION_Y))
    print("ry Node 3" , mp.Nodes[3].GetSolutionStepValue(REACTION_Y))
    print("ry Node 4" , mp.Nodes[4].GetSolutionStepValue(REACTION_Y))
    print("ry Node 5" , mp.Nodes[5].GetSolutionStepValue(REACTION_Y))
    print("ry Node 6" , mp.Nodes[6].GetSolutionStepValue(REACTION_Y))
    print("ry Node 7" , mp.Nodes[7].GetSolutionStepValue(REACTION_Y))
    print("ry Node 8" , mp.Nodes[8].GetSolutionStepValue(REACTION_Y))
    print("ry Node 9" , mp.Nodes[9].GetSolutionStepValue(REACTION_Y))
    print("ry Node 10" , mp.Nodes[10].GetSolutionStepValue(REACTION_Y))
    print("ry Node 11" , mp.Nodes[11].GetSolutionStepValue(REACTION_Y))
    print("ry Node 12" , mp.Nodes[12].GetSolutionStepValue(REACTION_Y))
    print("ry Node 13" , mp.Nodes[13].GetSolutionStepValue(REACTION_Y))
    print("ry Node 14" , mp.Nodes[14].GetSolutionStepValue(REACTION_Y))

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
        gid_io.FinalizeResults()
        gid_io.Reset()

    return mp

def test():
    mp = main(output=False)

    tol = 1e-16
    assert(abs(mp.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y) - mp.Nodes[5].GetSolutionStepValue(DISPLACEMENT_Y)) < tol)
    assert(abs(0.5*(mp.Nodes[6].GetSolutionStepValue(DISPLACEMENT_Y) + mp.Nodes[11].GetSolutionStepValue(DISPLACEMENT_Y)) - mp.Nodes[9].GetSolutionStepValue(DISPLACEMENT_Y)) < tol)
    print("Test passed")

def tag():
    return "MPC"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=True)
