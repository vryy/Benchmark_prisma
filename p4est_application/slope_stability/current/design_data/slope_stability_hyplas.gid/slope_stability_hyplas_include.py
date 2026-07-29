##################################################################
import sys
import math
import os
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.P4estApplication import *
kernel = Kernel()   #defining kernel
##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, output_path, logging=True ):
        #setting the domain size for the problem to be solved
        self.domain_size = 3
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("ekate_simulation")
        self.path = path + os.sep
        self.output_path = output_path + os.sep
        self.problem_name = problem_name
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################
        # reading simulation parameters
        number_of_time_steps = 1
        self.analysis_parameters = {}
        # content of analysis_parameters:
        # perform_contact_analysis_flag
        # penalty value for normal contact
        # maximum number of uzawa iterations
        # friction coefficient
        # penalty value for frictional contact
        # contact_double_check_flag
        # contact_ramp_penalties_flag
        # maximum penalty value for normal contact
        # ramp criterion for normal contact
        # ramp factor for normal contact
        # maximum penalty value for frictional contact
        # ramp criterion for frictional contact
        # ramp factor for frictional contact
        # analysis type: static (0), quasi-static (1) or dynamic (2)
        perform_contact_analysis_flag = False
        penalty = 0.0
        maxuzawa = 0.0
        friction = 0.0
        frictionpenalty = 0.0
        contact_double_check_flag = False
        contact_ramp_penalties_flag = False
        maxpenalty = 0.0
        rampcriterion = 0.0
        rampfactor = 0.0
        fricmaxpenalty = 0.0
        fricrampcriterion = 0.0
        fricrampfactor = 0.0
        self.analysis_parameters['perform_contact_analysis_flag'] = perform_contact_analysis_flag
        self.analysis_parameters['penalty'] = penalty
        self.analysis_parameters['maxuzawa'] = maxuzawa
        self.analysis_parameters['friction'] = friction
        self.analysis_parameters['frictionpenalty'] = frictionpenalty
        self.analysis_parameters['contact_double_check_flag'] = contact_double_check_flag
        self.analysis_parameters['contact_ramp_penalties_flag'] = contact_ramp_penalties_flag
        self.analysis_parameters['maxpenalty'] = maxpenalty
        self.analysis_parameters['rampcriterion'] = rampcriterion
        self.analysis_parameters['rampfactor'] = rampfactor
        self.analysis_parameters['fricmaxpenalty'] = fricmaxpenalty
        self.analysis_parameters['fricrampcriterion'] = fricrampcriterion
        self.analysis_parameters['fricrampfactor'] = fricrampfactor
        self.analysis_parameters['print_sparsity_info_flag'] = False
        self.analysis_parameters['analysis_type'] = 0
        self.analysis_parameters['dissipation_radius'] = 0.1
        self.analysis_parameters['decouple_build_and_solve'] = False
        self.analysis_parameters['builder_and_solver_type'] = "residual-based block with constraints"
        self.analysis_parameters['solving_scheme'] = 'monolithic'
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
        self.analysis_parameters['number_of_iterations_for_divergence_check'] = 7
        self.analysis_parameters['log_residuum'] = logging
        self.analysis_parameters['calculate_strain_energy'] = True
        self.analysis_parameters['log_strain_energy'] = logging
        print("builder_and_solver_type: " + str(self.analysis_parameters['builder_and_solver_type']))

        self.abs_tol =        1e-06
        self.rel_tol =        1e-10

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        #import ekate_solver_parallel
        #self.solver = ekate_solver_parallel.EkateSolver( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
        #ekate_solver_parallel.AddVariables( self.model_part )
        #variables used for discontinuities_application
        self.model_part.AddNodalSolutionStepVariable(LINE_LOAD)
        self.model_part.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.output_path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
        (self.solver).CalculateReactionFlag = True
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                cond_id = int(val_set[1]) #hbui added
                elem_id = int(val_set[2]) #hbui added
                if (cond_id in self.model_part.Conditions) and (elem_id in self.model_part.Elements):
                    self.model_part.Conditions[cond_id].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[elem_id].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
        print("input data read OK")
        #print "+++++++++++++++++++++++++++++++++++++++"
        #for node in self.model_part.Nodes:
        #    print node
        #print "+++++++++++++++++++++++++++++++++++++++"

        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        structural_solver_advanced.AddDofs( self.model_part )
        #ekate_solver_parallel.AddDofs( self.model_part )

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        if KratosMKLSolversApplication.Has("MKLPardisoSolver"):
            plinear_solver = MKLPardisoSolver()
        else:
            plinear_solver = SuperLUSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2);

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )

        self.material = "mohr-coulomb"

    def SetUpActivationLevels( self, model_part, activation_list, cond_activation_list ):
        for element in self.model_part.Elements:
            element.SetValue(ACTIVATION_LEVEL, activation_list[element.Id])
        for condition in self.model_part.Conditions:
            if( not (condition.GetValue(IS_TYING_MASTER) or condition.GetValue(IS_CONTACT_MASTER) ) ):
                condition.SetValue(ACTIVATION_LEVEL, activation_list[cond_activation_list[condition.Id-1]])

    def SetModelPart(self, model_part):
        self.model_part = model_part

        import structural_solver_advanced
        number_of_time_steps = 1
        self.solver = structural_solver_advanced.SolverAdvanced(self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol  )

        structural_solver_advanced.AddDofs( self.model_part )

        plinear_solver = SuperLUSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraints(plinear_solver)
        # (self.solver.solver).builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(plinear_solver)

        (self.solver.solver).SetEchoLevel(2)

    def WriteOutput( self, time ):
        mesh_dt = 1e-3 #adjust this parameter to adapt with simulation time step
        self.gid_io.InitializeMesh( time )
        mesh = self.model_part.GetMesh()
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, time, 0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part.Nodes, time, 0)
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(PLASTICITY_INDICATOR, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(INTEGRATION_POINT_STRAIN_VECTOR, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(ELASTIC_STRAIN_VECTOR, self.model_part, time)
        self.gid_io.FinalizeResults()

    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        self.model_part.Properties[1].SetValue(BODY_FORCE, [0, 0, 0])
        self.model_part.Properties[1].SetValue(GRAVITY, [0, -9.81, 0])
        self.model_part.Properties[1].SetValue(DENSITY,       2.038735 )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,        20000.0 )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,         0.49 )
        self.model_part.Properties[1].SetValue(THICKNESS,            1.0 )
        self.model_part.Properties[1].SetValue(INTERNAL_FRICTION_ANGLE,       20.0 )
        self.model_part.Properties[1].SetValue(DILATANCY_ANGLE,            20.0 )
        self.model_part.Properties[1].SetValue(COHESION,            50.0 )
        self.model_part.Properties[1].SetValue(ISOTROPIC_HARDENING_MODULUS,     0.0 )
        if self.material == "mohr-coulomb":
            self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, MohrCoulomb3dImplicit() )
            print("Mohr Coulomb plasticity model selected")
        elif self.material == "drucker-prager":
            self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, DruckerPrager3dImplicit_InnerEdgeMatchWithMohrCoulomb() )
            print("Drucker Prager plasticity model selected")
        elif self.material == "plane-strain":
            self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
            print("Plane strain elasticity model selected")
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        #self.SetUpActivationLevels( self.model_part, self.activation_flags, self.cond_activation_flags )
        self.deac.Initialize( self.model_part )
        self.solver.solver.Initialize()
        self.model_part.Check(self.model_part.ProcessInfo) #important to check the negativeness of Jacobian
        print("activation utility initialized")
        ##################################################################
        ## MESH TYING ####################################################
        ##################################################################
        #self.mesh_tying_utility= MeshTyingUtility()
        ##self.mesh_tying_utility.InitializeMeshTyingUtilityLagrange(self.model_part)
        #self.mesh_tying_utility.InitializeMeshTyingUtility(self.model_part)
        #print "mesh-tying utility successfully initialized"
        #self.model_part.Check( self.model_part.ProcessInfo )
        print("model successfully initialized")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
##################################################################
