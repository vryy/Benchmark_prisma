##################################################################
import sys
import os
import six
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.P4estApplication import *
kernel = Kernel()   #defining kernel
##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, results_path, logging=True ):
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart(problem_name)
        self.path = path+os.sep
        self.results_path = results_path+os.sep
        self.problem_name = problem_name
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################
        # reading simulation parameters
        self.analysis_parameters = {}
        # analysis type: static (0), quasi-static (1) or dynamic (2)
        self.analysis_parameters['print_sparsity_info_flag'] = False
        self.analysis_parameters['analysis_type'] = 0
        self.analysis_parameters['dissipation_radius'] =          0.1
        self.analysis_parameters['decouple_build_and_solve'] = False
        self.analysis_parameters['builder_and_solver_type'] = "residual-based block with constraints"
        # self.analysis_parameters['builder_and_solver_type'] = "residual-based block with constraints element-wise"
        self.analysis_parameters['solving_scheme'] = 'monolithic'
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
        self.analysis_parameters['list_dof'] = True
        self.analysis_parameters['log_residuum'] = logging


        self.abs_tol =        1e-10
        self.rel_tol =        1e-06

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.analysis_parameters, self.abs_tol, self.rel_tol )
        self.model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(LAGRANGE_TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(TEMPERATURE_ERROR)
        self.model_part.AddNodalSolutionStepVariable(REFERENCE_TEMPERATURE)
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = SDGidPostIO( self.results_path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
        (self.solver).CalculateReactionFlag = False
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        self.element_assignments = {}
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
                self.element_assignments[int(val_set[1])] = int(val_set[2])
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
        for node in self.model_part.Nodes:
            node.AddDof(TEMPERATURE)
            node.AddDof(LAGRANGE_TEMPERATURE)

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
        (self.solver.solver).SetEchoLevel(2)
        (self.solver.solver).max_iter = 10 #control the maximum iterations of Newton Raphson loop
        (self.solver.solver).MoveMeshFlag = False
        self.solver.solver.convergence_criteria = DisplacementCriteria(self.rel_tol, self.abs_tol)
        self.solver.solver.ReformDofSetAtEachStep = False

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )

    def SetModelPart(self, model_part):
        print("Calling SetModelPart...")
        self.model_part = model_part
        number_of_time_steps = 1

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.analysis_parameters, self.abs_tol, self.rel_tol )

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        for node in self.model_part.Nodes:
            node.AddDof(TEMPERATURE)
            node.AddDof(LAGRANGE_TEMPERATURE)

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
        (self.solver.solver).SetEchoLevel(2)
        (self.solver.solver).max_iter = 10 #control the maximum iterations of Newton Raphson loop
        (self.solver.solver).MoveMeshFlag = False
        self.solver.solver.convergence_criteria = DisplacementCriteria(self.rel_tol, self.abs_tol)
        self.solver.solver.ReformDofSetAtEachStep = False
        print("Calling SetModelPart...completed")

    def WriteOutput( self, time ):
        self.gid_io.InitializeMesh( time )
        mesh = self.model_part.GetMesh()
        print("mesh at time " + str(time) + " is ready for printing")
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
        self.gid_io.WriteNodalResults(TEMPERATURE, self.model_part.Nodes, time, 0)
        self.gid_io.WriteNodalResults(TEMPERATURE_ERROR, self.model_part.Nodes, time, 0)
        self.gid_io.WriteNodalResults(REFERENCE_TEMPERATURE, self.model_part.Nodes, time, 0)
        self.gid_io.FinalizeResults()
        self.gid_io.Reset()

    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        self.model_part.Properties[1].SetValue(INTEGRATION_ORDER, 2 )
        self.model_part.Properties[1].SetValue(THERMAL_CONDUCTIVITY, 1.0 )
        self.model_part.Properties[1].SetValue(THICKNESS,            1.0 )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW,  DummyConstitutiveLaw() )
        util = ThermalAuxiliaryUtility()
        util.SetValueForProperties(self.model_part.Properties[1], HEAT_SOURCE, HeatSourceStdProblem1())
        for prop in self.model_part.Properties:
            prop.SetValue(THICKNESS, 1.0)
            util.SetValueForProperties(prop, TEMPERATURE_FUNCTION, HeatStdProblem1Solution())
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
        self.solver.solver.Initialize()
        self.model_part.Check( self.model_part.ProcessInfo )
        print("activation utility initialized")
        print("model successfully initialized")

    def WriteRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time)
        serializer = Serializer(fn)
        serializer.Save("ModelPart", self.model_part)
        serializer = 0
        print("Write restart data to " + fn + ".rest completed")

    def LoadRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time)
        serializer = Serializer(fn)
        serializer.Load("ModelPart", self.model_part)
        serializer = 0
        print("Load restart data from " + fn + ".rest completed")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
