##################################################################
from __future__ import absolute_import
import sys
import os
import six # for iteritems
import time as time_module
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.ErsatzAnwendung import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *
kernel = Kernel()   #defining kernel
##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, results_path, logging=True, decouple_build_and_solve=True ):
        #setting the domain size for the problem to be solved
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
        self.analysis_parameters['analysis_type'] = -1 # using custom time scheme
        self.analysis_parameters['decouple_build_and_solve'] = decouple_build_and_solve
        self.analysis_parameters['reform_dofset_at_each_step'] = False
        self.analysis_parameters['solving_scheme'] = 'monolithic'
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
        self.analysis_parameters['list_dof'] = True
        self.analysis_parameters['log_residuum'] = logging
        self.analysis_parameters['time_scheme'] = ResidualBasedTemperatureBasedBackwardEulerScheme()
        self.analysis_parameters['convergence_criteria'] = "displacement"

        self.abs_tol =        1e-06
        self.rel_tol =        1e-10

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.analysis_parameters, self.abs_tol, self.rel_tol )
        self.AddVariables( self.model_part )
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
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
                cond_id = int(val_set[1])
                elem_id = int(val_set[2])
                if elem_id in self.model_part.Elements:
                    elem = self.model_part.Elements[elem_id]
                else:
                    raise Exception(f"Element {elem_id} does not exist")
                self.model_part.Conditions[cond_id].SetValue( ACTIVATION_LEVEL, elem.GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(elem_id) + " to Condition: " + str(cond_id) + " as " + str(elem.GetValue(ACTIVATION_LEVEL)) )
                self.element_assignments[cond_id] = elem_id
        print("input data read OK")
        #print("+++++++++++++++++++++++++++++++++++++++")
        #for node in self.model_part.Nodes:
        #    print(node)
        #print("+++++++++++++++++++++++++++++++++++++++")
        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## POST_PROCESSING ###############################################
        ##################################################################

        self.write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        self.write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        self.post_mode = GiDPostMode.GiD_PostBinary
        self.multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = SDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )

        # define the output variables
        self.output_nodal_variables = []
        self.output_intpt_variables = []

        self.output_nodal_variables.append(TEMPERATURE)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = SuperLUSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)
        (self.solver.solver).max_iter = 10 #control the maximum iterations of Newton Raphson loop
        (self.solver.solver).MoveMeshFlag = False

        # replace the builder and solver by the POD type
        factory = ProjectionBasedPodBuilderAndSolverFactory()
        (self.solver.solver).builder_and_solver = factory.Create((self.solver.solver).builder_and_solver)
        # (self.solver.solver).builder_and_solver = ProjectionBasedPodBuilderAndSolver((self.solver.solver).builder_and_solver)

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )


    def AddDofsForNodes(self, nodes):
        import structural_solver_advanced
        structural_solver_advanced.AddDofsForNodes( nodes )

    def AddDofs(self, model_part):
        import structural_solver_advanced
        structural_solver_advanced.AddDofs( model_part )
        for node in model_part.Nodes:
            node.AddDof(TEMPERATURE, ENTROPY)

    def AddDofsForNode(self, node):
        import structural_solver_advanced
        structural_solver_advanced.AddDofsForNode( node )
        node.AddDof(TEMPERATURE, ENTROPY)

    def AddVariables(self, model_part):
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( model_part )
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_NULL)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_EINS)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_DT)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_NULL_DT)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_EINS_DT)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_NULL_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(TEMPERATURE_EINS_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(ENTROPY)

    def SetModelPart(self, model_part):
        self.model_part = model_part

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.analysis_parameters, self.abs_tol, self.rel_tol )
        (self.solver).CalculateReactionFlag = False
        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = SuperLUSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)

    def SetOutputPath(self, path):
        self.results_path = path
        self.gid_io = SDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )


    def WriteOutput( self, time ):
        self.gid_io.InitializeMesh( time )
        mesh = self.model_part.GetMesh()
        print("mesh at time " + str(time) + " is ready for printing")
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
        for var in self.output_nodal_variables:
            self.gid_io.WriteNodalResults(var, time, 0)
            print(f"nodal {var} written")
        for var in self.output_intpt_variables:
            self.gid_io.PrintOnGaussPoints(var, self.model_part, time)
            print(f"gauss point {var} written")
        self.gid_io.FinalizeResults()
        self.gid_io.Reset()

    def InitializeModel( self ):
        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        model_layers = __import__(self.problem_name+"_layers")
        ## ELEMENTS on layers ############################################
        self.layer_sets = model_layers.ReadLayerSets()
        ## NODES on layers ###############################################
        self.layer_nodes_sets = model_layers.ReadLayerNodesSets()
        ## CONTACT MASTER NODES ##########################################
        #self.contact_master_nodes = model_layers.ReadContactMasterNodes()
        ## CONTACT SLAVE NODES ###########################################
        #self.contact_slave_nodes = model_layers.ReadContactSlaveNodes()
        ##################################################################
        print("layer sets stored")
        ##################################################################
        ## STORE NODES ON GROUND SURFACE #################################
        ##################################################################
        self.top_surface_nodes = model_layers.ReadTopSurfaceNodes()
        print("nodes on ground surface stored")
        ##################################################################
        ## STORE NODES ON SIDE ###########################################
        ##################################################################
        self.boundary_nodes = model_layers.ReadBoundaryNodes()
        print("nodes on side surface stored")
        ##################################################################
        ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        ##################################################################
        self.node_groups = model_layers.ReadNodeGroups()
        print("node groups stored")
        ##################################################################
        ## EXTRACT CONDITIONS FROM NODE GROUPS ###########################
        ##################################################################
        start_time = time_module.time()
        self.layer_cond_sets = {}
        for layer, node_group in six.iteritems(self.node_groups):
            self.layer_cond_sets[layer] = []
        for layer, node_group in six.iteritems(self.node_groups):
            for cond in self.model_part.Conditions:
                in_group = True
                for node in cond.GetNodes():
                    if node.Id not in node_group:
                        in_group = False
                        break
                if in_group:
                    self.layer_cond_sets[layer].append(cond.Id)
        end_time = time_module.time()
        print("conditions in node groups stored, time = " + str(end_time - start_time) + "s")
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        aux_util = ThermalAuxiliaryUtility()
        self.model_part.Properties[1].SetValue(THERMAL_CONDUCTIVITY, 0.01 )
        self.model_part.Properties[1].SetValue(DENSITY, 1.0 )
        self.model_part.Properties[1].SetValue(SPECIFIC_HEAT_CAPACITY, 1.0 )
        self.model_part.Properties[1].SetValue(THICKNESS, 1.0 )
        aux_util.SetValueForProperties(self.model_part.Properties[1], HEAT_SOURCE,    NullHeatSource() )
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

    # solve with deactivation/reactivation
    # element/condition with nonzero ACTIVATION_LEVEL in [from_deac, to_deac] will be deactivated
    # element/condition with negative ACTIVATION_LEVEL will also be deactivated
    # element/condition with ACTIVATION_LEVEL in [from_reac, to_reac] will be re-activated
    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()

    # solve nothing (good for debugging) with deactivation/reactivation
    def DrySolve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)

    # solve without deactivation
    def SolveModel(self, time):
        self.model_part.CloneTimeStep(time)
        success = self.solver.Solve()
        return success

    # solve nothing without deactivation
    def DrySolveModel(self, time):
        self.model_part.CloneTimeStep(time)

##################################################################
