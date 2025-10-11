##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
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
from KratosMultiphysics.OofemApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *
kernel = Kernel()   #defining kernel
##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, results_path, logging=True ):
        #setting the domain size for the problem to be solved
        self.domain_size = 3
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart(problem_name)
        self.path = path
        self.results_path = results_path+os.sep
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
        self.analysis_parameters['dimension'] = self.domain_size # 3
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
        self.analysis_parameters['decouple_build_and_solve'] = True
        self.analysis_parameters['solving_scheme'] = 'monolithic'
        self.analysis_parameters['residuum_tolerance'] = 5e-3
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
        self.analysis_parameters['list_dof'] = True
        self.analysis_parameters['log_residuum'] = logging

        self.abs_tol =        1e-06
        self.rel_tol =        1e-6

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        #import ekate_solver_parallel
        #self.solver = ekate_solver_parallel.EkateSolver( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        self.AddVariables( self.model_part )
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
        (self.solver).CalculateReactionFlag = True
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

        self.output_nodal_variables.append(DISPLACEMENT)
        self.output_nodal_variables.append(REACTION)

        self.output_intpt_variables.append(STRESSES)
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
        (self.solver.solver).max_iter = 20 #control the maximum iterations of Newton Raphson loop
        (self.solver.solver).MoveMeshFlag = False

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

    def AddDofsForNode(self, node):
        import structural_solver_advanced
        structural_solver_advanced.AddDofsForNode( node )

    def AddVariables(self, model_part):
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( model_part )

    def SetModelPart(self, model_part):
        self.model_part = model_part
        number_of_time_steps = 1

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        (self.solver).CalculateReactionFlag = True
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
        #
        # self.model_part.Properties[1].SetValue(DENSITY,         0.0 )
        # self.model_part.Properties[1].SetValue(YOUNG_MODULUS,      30e9 )
        # self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.15 )
        # self.model_part.Properties[1].SetValue(THICKNESS, 1.0 )
        # self.model_part.Properties[1].SetValue(INTEGRATION_ORDER, 1 )
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, OofemIsotropicLinearElasticMaterial3D() )
        # print("Linear elastic model selected for LinearElastic, description: Bulk")
        #
        # self.model_part.Properties[1].SetValue(IFT_Material_density,         0.0 )
        self.model_part.Properties[1].SetValue(IFT_IsotropicLinearElasticMaterial_e,         30e9 )
        self.model_part.Properties[1].SetValue(IFT_IsotropicLinearElasticMaterial_n,         0.15 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_wf,         9.3755e-4 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_fc,         3e6 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_ft,         1e6 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_hp,         0.01 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_yieldtol,         1e-5 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_asoft,         5.0 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_softeningType,         1 )
        self.model_part.Properties[1].SetValue(IFT_ConcreteDPM2_kinit,         0.3 )
        self.model_part.Properties[1].SetValue(INTEGRATION_ORDER, 1 )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, OofemConcreteDPM2() )
        self.model_part.Properties[1].SetValue(MAT_RESPONSE_MODE, 1 ) # 0: tangent; 1: secant; 3: elastic
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, OofemConcreteDPM2DC() )
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, ElasticTangentConstitutiveLaw( OofemConcreteDPM2() ) )
        print("ConcreteDPM2 model selected for bulk, description: Bulk")
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
        self.solver.Solve()

    # solve nothing without deactivation
    def DrySolveModel(self, time):
        self.model_part.CloneTimeStep(time)

##################################################################
