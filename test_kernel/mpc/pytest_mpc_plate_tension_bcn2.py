from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.LayerApplication import *
# from inspect import currentframe, getframeinfo

def main(output=True):
    # Define a Model
    mp = ModelPart("Main")

    mp.AddNodalSolutionStepVariable(DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(REACTION)
    mp.AddNodalSolutionStepVariable(FACE_LOAD)
    mp.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)

    # Create Nodes
    mp.CreateNewNode(1, 0.00,0.00,0.00)
    mp.CreateNewNode(2, 1.00,0.00,0.00)
    mp.CreateNewNode(3, 1.00,1.00,0.00)
    mp.CreateNewNode(4, 0.00,1.00,0.00)
    mp.CreateNewNode(5, 0.50,0.00,0.00)
    mp.CreateNewNode(6, 1.00,0.50,0.00)
    mp.CreateNewNode(7, 0.50,1.00,0.00)
    mp.CreateNewNode(8, 0.50,0.50,0.00)
    # define properties
    mp.GetProperties()[1].SetValue(YOUNG_MODULUS, 210e3)
    mp.GetProperties()[1].SetValue(POISSON_RATIO, 0.3)
    mp.GetProperties()[1].SetValue(THICKNESS, 1.0)
    mp.GetProperties()[1].SetValue(CONSTITUTIVE_LAW, PlaneStrain())
    mp.GetProperties()[1].SetValue(BODY_FORCE, ZeroVector(3))
    # Creat Elements
    mp.CreateNewElement("KinematicLinear2D4N", 1, [1,5,7,4], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 2, [5,2,6,8], mp.GetProperties()[1])
    mp.CreateNewElement("KinematicLinear2D4N", 3, [8,6,3,7], mp.GetProperties()[1])
    # Add Dofs
    for node in mp.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)

    # mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, mp.Nodes[8], DISPLACEMENT_X, mp.Nodes[7], DISPLACEMENT_X, 1.0, 0)
    # mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, mp.Nodes[8], DISPLACEMENT_X, mp.Nodes[5], DISPLACEMENT_X, 1.0, 0)

    relation_matrix = Matrix(1, 2)
    relation_matrix[0, 0] = 0.5
    relation_matrix[0, 1] = 0.5
    constant_vector = ZeroVector(1)

    mp.CreateNewLinearMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, mp.Nodes[8], [mp.Nodes[7], mp.Nodes[5]], DISPLACEMENT_X, relation_matrix, constant_vector)
    mp.CreateNewLinearMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, mp.Nodes[8], [mp.Nodes[7], mp.Nodes[5]], DISPLACEMENT_Y, relation_matrix, constant_vector)

    # impose BC
    mp.Nodes[1].Fix(DISPLACEMENT_X)
    mp.Nodes[1].Fix(DISPLACEMENT_Y)
    mp.Nodes[2].Fix(DISPLACEMENT_Y)
    mp.Nodes[4].Fix(DISPLACEMENT_X)
    # Load
    load = Vector(3)
    load[0] = -1000
    load[1] = 0.0
    load[2] = 0.0

    condition1 = mp.CreateNewCondition("LinePressure2D2N", 1, [2,6], mp.GetProperties()[1])
    condition2 = mp.CreateNewCondition("LinePressure2D2N", 2, [6,3], mp.GetProperties()[1])
    condition1.SetValue(PRESSURE, load[0])
    condition2.SetValue(PRESSURE, load[0])

    linear_solver = SkylineLUFactorizationSolver()
    builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraints(linear_solver)
    scheme = ResidualBasedIncrementalUpdateStaticScheme()
    convergence_criterion = ResidualCriteria(1e-10, 1e-12)

    max_iters = 100
    compute_reactions = True
    reform_step_dofs = False
    move_mesh_flag = False
    strategy = ResidualBasedNewtonRaphsonStrategy(mp, scheme, linear_solver, convergence_criterion, builder_and_solver, max_iters, compute_reactions,reform_step_dofs, move_mesh_flag)
    # strategy = ResidualBasedLinearStrategy(mp, scheme, linear_solver, builder_and_solver,False,False,False,False)
    strategy.SetEchoLevel(2)
    strategy.Initialize()
    strategy.Check()
    strategy.Solve()

    print("ux Node 1 =" , mp.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 2 =" , mp.Nodes[2].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 3 =" , mp.Nodes[3].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 4 =" , mp.Nodes[4].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 5 =" , mp.Nodes[5].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 6 =" , mp.Nodes[6].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 7 =" , mp.Nodes[7].GetSolutionStepValue(DISPLACEMENT_X))
    print("ux Node 8 =" , mp.Nodes[8].GetSolutionStepValue(DISPLACEMENT_X))
    print("uy Node 1 =" , mp.Nodes[1].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 2 =" , mp.Nodes[2].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 3 =" , mp.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 4 =" , mp.Nodes[4].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 5 =" , mp.Nodes[5].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 6 =" , mp.Nodes[6].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 7 =" , mp.Nodes[7].GetSolutionStepValue(DISPLACEMENT_Y))
    print("uy Node 8 =" , mp.Nodes[8].GetSolutionStepValue(DISPLACEMENT_Y))

    print("rx Node 1" , mp.Nodes[1].GetSolutionStepValue(REACTION_X))
    print("rx Node 2" , mp.Nodes[2].GetSolutionStepValue(REACTION_X))
    print("rx Node 3" , mp.Nodes[3].GetSolutionStepValue(REACTION_X))
    print("rx Node 4" , mp.Nodes[4].GetSolutionStepValue(REACTION_X))
    print("rx Node 5" , mp.Nodes[5].GetSolutionStepValue(REACTION_X))
    print("rx Node 6" , mp.Nodes[6].GetSolutionStepValue(REACTION_X))
    print("rx Node 7" , mp.Nodes[7].GetSolutionStepValue(REACTION_X))
    print("rx Node 8" , mp.Nodes[8].GetSolutionStepValue(REACTION_X))

    print("ry Node 1" , mp.Nodes[1].GetSolutionStepValue(REACTION_Y))
    print("ry Node 2" , mp.Nodes[2].GetSolutionStepValue(REACTION_Y))
    print("ry Node 3" , mp.Nodes[3].GetSolutionStepValue(REACTION_Y))
    print("ry Node 4" , mp.Nodes[4].GetSolutionStepValue(REACTION_Y))
    print("ry Node 5" , mp.Nodes[5].GetSolutionStepValue(REACTION_Y))
    print("ry Node 6" , mp.Nodes[6].GetSolutionStepValue(REACTION_Y))
    print("ry Node 7" , mp.Nodes[7].GetSolutionStepValue(REACTION_Y))
    print("ry Node 8" , mp.Nodes[8].GetSolutionStepValue(REACTION_Y))

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

    print("Analysis completed")

    return mp

def test():
    mp = main(output=False)

    tol = 1e-16
    assert(abs(mp.Nodes[2].GetSolutionStepValue(DISPLACEMENT_Y) - mp.Nodes[1].GetSolutionStepValue(DISPLACEMENT_Y)) < tol)
    assert(abs(0.5*(mp.Nodes[2].GetSolutionStepValue(DISPLACEMENT_Y) + mp.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)) - mp.Nodes[6].GetSolutionStepValue(DISPLACEMENT_Y)) < tol)
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
