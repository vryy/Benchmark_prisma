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
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.ThermalStructuralApplication import *
from KratosMultiphysics.ThermalSoilMechanicsApplication import *
from KratosMultiphysics.MKLSolversApplication import *
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
        self.analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
        self.analysis_parameters['list_dof'] = True
        self.analysis_parameters['log_residuum'] = logging

        self.abs_tol =        1e-10
        self.rel_tol =        1e-10

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
        self.model_part.AddNodalSolutionStepVariable(ENTROPY)
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.results_path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
        (self.solver).CalculateReactionFlag = True
        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs( self.model_part )

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

    def FixPressureNodes( self, free_node_list_water, free_node_list_air):
        for node in self.model_part.Nodes:
            if (node.IsFixed(WATER_PRESSURE)==0):
                node.Fix(WATER_PRESSURE)
                free_node_list_water.append(node)
            if (node.IsFixed(AIR_PRESSURE)==0):
                node.Fix(AIR_PRESSURE)
                free_node_list_air.append(node)

    def ApplyInsituWaterPressure( self, free_node_list_water, free_node_list_air, z_zero, gravity_z):
        water_density=1000.0;
        for node in self.model_part.Nodes:
            water_pressure= water_density*gravity_z*(z_zero-(node.Z-node.GetSolutionStepValue(DISPLACEMENT_Z,0)))
            if( water_pressure < 1.0 ):
                water_pressure = 1.0
            node.SetSolutionStepValue(WATER_PRESSURE, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_EINS, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_NULL, water_pressure)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

    def SetReferenceWaterPressure( self ):
        for element in self.model_part.Elements:
            water_pressures = element.GetValuesOnIntegrationPoints( WATER_PRESSURE, self.model_part.ProcessInfo )
            pressure_list = []
            for item in water_pressures:
                pressure_list.append( item[0] )
            element.SetValuesOnIntegrationPoints( REFERENCE_WATER_PRESSURE, pressure_list, self.model_part.ProcessInfo )

    def FreePressureNodes(self,free_node_list_water, free_node_list_air):
        for item in free_node_list_water:
            #self.model_part.Nodes[item].Free(WATER_PRESSURE)
            item.Free(WATER_PRESSURE)
        for item in free_node_list_air:
            #self.model_part.Nodes[item].Free(AIR_PRESSURE)
            item.Free(AIR_PRESSURE)

    def WriteMaterialParameters( self, time, indices ):
        self.gid_io.OpenResultFile( self.results_path+self.problem_name, GiDPostMode.GiD_PostBinary)
        #self.gid_io.ChangeOutputName( self.results_path+self.problem_name +str(time), GiDPostMode.GiD_PostBinary )
        for index in indices:
            self.gid_io.SuperPrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, index)
        self.gid_io.CloseResultFile()

    def WriteMonitoringSectionResults( self, time ):
        outfile = open("step_"+str(time)+".dat",'w')
        outfile.write("ekate result file for step "+str(time)+"\n")
        outfile.close()

    def WriteOutput( self, time ):
        self.gid_io.InitializeMesh( time )
        mesh = self.model_part.GetMesh()
        print("mesh at time " + str(time) + " is ready for printing")
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 0)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 1)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 2)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 3)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 4)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 5)
        #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 6)
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, time, 0)
        self.gid_io.WriteNodalResults(TEMPERATURE, self.model_part.Nodes, time, 0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part.Nodes, time, 0)
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        print("gauss point STRESSES written")
        #self.gid_io.PrintOnGaussPoints(CONTACT_PENETRATION, self.model_part, time)
        #self.gid_io.PrintOnGaussPoints(NORMAL, self.model_part, time, 0)
        self.gid_io.FinalizeResults()

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
        self.layer_cond_sets = {}
        for layer, node_group in self.node_groups.items():
            self.layer_cond_sets[layer] = []
        for layer, node_group in self.node_groups.items():
            for cond in self.model_part.Conditions:
                in_group = True
                for node in cond.GetNodes():
                    if node.Id not in node_group:
                        in_group = False
                        break
                if in_group:
                    self.layer_cond_sets[layer].append(cond.Id)
        print("conditions in node groups stored")
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        K       = 801938.0
        G       = 164206.0
        E       = 9*K*G/(3*K+G)
        nu      = (3*K-2*G)/2/(3*K+G)
        y0      = 450.0
        h       = 129.24
        y0inf   = 715.0
        delta   = 16.93
        density = 7.8e-9
        alpha   = 1e-5
        kappa   = 45.0
        c       = 0.46e9
        w0      = 0.002
        wh      = 0.002
        T0      = 0.0
        self.model_part.Properties[1].SetValue(DENSITY,            density )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,      E )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          nu )
        self.model_part.Properties[1].SetValue(THICKNESS, 1.0 )
        self.model_part.Properties[1].SetValue(INTEGRATION_ORDER, 2 )
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, VonMises3dImplicitDC() )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, VonMises3dtImplicit() )
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, VonMises3dtImplicitDC() )
        self.model_part.Properties[1].SetValue(THERMAL_CONDUCTIVITY, kappa )
        self.model_part.Properties[1].SetValue(THERMAL_EXPANSION, alpha )
        self.model_part.Properties[1].SetValue(SPECIFIC_HEAT_CAPACITY, c ) # not important for steady analysis
        self.model_part.Properties[1].SetValue(TAYLOR_QUINNEY_FACTOR, 0.9 ) # to account for heat generated by plasticity
        hlaw = LinearThermalHardeningLaw(y0, y0inf, h, w0, wh, delta, T0)
        hlaw.Assign(ISOTROPIC_HARDENING_LAW, hlaw, self.model_part.Properties[1])
        self.material_type = "j2-thermal-im"
        print("J2 thermal model (implicit) is selected")
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
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
##################################################################
