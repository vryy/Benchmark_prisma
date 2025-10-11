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
import sys
import os
import time as time_module
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FractureApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.LayerApplication import *
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, results_path, logging=True ):
        print("logging:", logging)
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
        #print "+++++++++++++++++++++++++++++++++++++++"
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

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = MKLPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)
        (self.solver.solver).max_iter = 10 #control the maximum iterations of Newton Raphson loop
        (self.solver.solver).MoveMeshFlag = True
        self.solver.solver.convergence_criteria = DisplacementCriteria(self.rel_tol, self.abs_tol)

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )


    def AddDofsForNodes(self, nodes):
        import structural_solver_advanced
        structural_solver_advanced.AddDofsForNodes( nodes )
        for node in nodes:
            node.AddDof(DAMAGE, DAMAGE_REACTION)
            node.AddDof(PLASTICITY, PLASTICITY_REACTION)

    def AddDofs(self, model_part):
        import structural_solver_advanced
        structural_solver_advanced.AddDofs( model_part )
        for node in model_part.Nodes:
            node.AddDof(DAMAGE, DAMAGE_REACTION)
            node.AddDof(PLASTICITY, PLASTICITY_REACTION)

    def AddDofsForNode(self, node):
        import structural_solver_advanced
        structural_solver_advanced.AddDofsForNode( node )
        node.AddDof(DAMAGE, DAMAGE_REACTION)
        node.AddDof(PLASTICITY, PLASTICITY_REACTION)

    def AddVariables(self, model_part):
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( model_part )
        model_part.AddNodalSolutionStepVariable(DAMAGE)
        model_part.AddNodalSolutionStepVariable(DAMAGE_REACTION)
        model_part.AddNodalSolutionStepVariable(PLASTICITY)
        model_part.AddNodalSolutionStepVariable(PLASTICITY_REACTION)

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
        plinear_solver = MKLPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)

    def SetOutputPath(self, path):
        self.results_path = path
        self.gid_io = SDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )

    def FixPressureNodes( self, free_node_list_water, free_node_list_air):
        for node in self.model_part.Nodes:
            if (node.IsFixed(WATER_PRESSURE)==0):
                node.Fix(WATER_PRESSURE)
                free_node_list_water.append(node)
            if (node.IsFixed(AIR_PRESSURE)==0):
                node.Fix(AIR_PRESSURE)
                free_node_list_air.append(node)

    def ApplyInsituWaterPressure( self, free_node_list_water, free_node_list_air, z_zero, gravity_z, water_density=1000.0 ):
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

    def ApplyInsituWaterPressure2D( self, free_node_list_water, free_node_list_air, y_zero, gravity_y, water_density=1000.0 ):
        for node in self.model_part.Nodes:
            water_pressure= water_density*gravity_y*(y_zero-(node.Y-node.GetSolutionStepValue(DISPLACEMENT_Y,0)))
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
        self.SetReferenceWaterPressureForElements(self.model_part.Elements)

    def SetReferenceWaterPressureForElements( self, elements ):
        for element in elements:
            self.SetReferenceWaterPressureForElement(element)

    def SetReferenceWaterPressureForElement( self, element ):
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
        self.gid_io.WriteNodalResults(DISPLACEMENT, time, 0)
        print("nodal DISPLACEMENT written")
        self.gid_io.WriteNodalResults(REACTION, time, 0)
        print("nodal REACTION written")
        self.gid_io.WriteNodalResults(DAMAGE, time, 0)
        print("nodal DAMAGE written")
        self.gid_io.WriteNodalResults(PLASTICITY, time, 0)
        print("nodal PLASTICITY written")
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        print("gauss point STRESSES written")
        self.gid_io.PrintOnGaussPoints(LOCAL_DAMAGE, self.model_part, time)
        print("gauss point LOCAL_DAMAGE written")
        self.gid_io.PrintOnGaussPoints(PLASTICITY_INDICATOR, self.model_part, time)
        print("gauss point PLASTICITY_INDICATOR written")
        self.gid_io.PrintOnGaussPoints(YIELD_STATE, self.model_part, time)
        print("gauss point YIELD_STATE written")
        self.gid_io.FinalizeResults()
        self.gid_io.Reset()

    def InitializeModel( self ):
        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        # model_layers = __import__(self.problem_name+"_layers")
        # ## ELEMENTS on layers ############################################
        # self.layer_sets = model_layers.ReadLayerSets()
        # ## NODES on layers ###############################################
        # self.layer_nodes_sets = model_layers.ReadLayerNodesSets()
        # ## CONTACT MASTER NODES ##########################################
        # #self.contact_master_nodes = model_layers.ReadContactMasterNodes()
        # ## CONTACT SLAVE NODES ###########################################
        # #self.contact_slave_nodes = model_layers.ReadContactSlaveNodes()
        # ##################################################################
        # print("layer sets stored")
        # ##################################################################
        # ## STORE NODES ON GROUND SURFACE #################################
        # ##################################################################
        # self.top_surface_nodes = model_layers.ReadTopSurfaceNodes()
        # print("nodes on ground surface stored")
        # ##################################################################
        # ## STORE NODES ON SIDE ###########################################
        # ##################################################################
        # self.boundary_nodes = model_layers.ReadBoundaryNodes()
        # print("nodes on side surface stored")
        # ##################################################################
        # ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        # ##################################################################
        # self.node_groups = model_layers.ReadNodeGroups()
        # print("node groups stored")
        # ##################################################################
        # ## EXTRACT CONDITIONS FROM NODE GROUPS ###########################
        # ##################################################################
        # start_time = time_module.time()
        # self.layer_cond_sets = {}
        # for layer, node_group in self.node_groups.iteritems():
        #     self.layer_cond_sets[layer] = []
        # for layer, node_group in self.node_groups.iteritems():
        #     for cond in self.model_part.Conditions:
        #         in_group = True
        #         for node in cond.GetNodes():
        #             if node.Id not in node_group:
        #                 in_group = False
        #                 break
        #         if in_group:
        #             self.layer_cond_sets[layer].append(cond.Id)
        # end_time = time_module.time()
        # print("conditions in node groups stored, time = " + str(end_time - start_time) + "s")
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        self.model_part.Properties[1].SetValue(DENSITY,         7.8e-9 )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,      2.1e+5 )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.3 )
        self.model_part.Properties[1].SetValue(THICKNESS,            1 )
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
        # print("Linear elastic model selected")
        mat_params = Vector(8)
        mat_params[0] = 70000.
        mat_params[1] = 0.33
        mat_params[2] = 1.5     # k1
        mat_params[3] = 1.0     # k2
        mat_params[4] = 2.25    # k3
        mat_params[5] = 0.      # fN
        mat_params[6] = 0.119      # sN
        mat_params[7] = 0.183115   # eN
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, GradientEnhancedDamagePlasticityMultiplicativeFiniteStrainBridgingConstitutiveLaw_Cauchy(GradientEnhancedDamagePlasticityGTNLaw()) )
        self.model_part.Properties[1].SetValue(GRADIENT_DAMAGE_PENALIZATION, 1.0 ) # beta_f
        self.model_part.Properties[1].SetValue(GRADIENT_PLASTICITY_PENALIZATION, 1.0e1 ) # beta_p
        self.model_part.Properties[1].SetValue(LOCAL_NONLOCAL_DAMAGE_PENALIZATION, 1.0 ) # alpha_f
        self.model_part.Properties[1].SetValue(LOCAL_NONLOCAL_PLASTICITY_PENALIZATION, 1.0e1 ) # alpha_p
        ## enable this to enable local GTN
        # self.model_part.Properties[1].SetValue(GRADIENT_DAMAGE_PENALIZATION, 0.0 ) # beta_f
        # self.model_part.Properties[1].SetValue(GRADIENT_PLASTICITY_PENALIZATION, 0.0 ) # beta_p
        # self.model_part.Properties[1].SetValue(LOCAL_NONLOCAL_DAMAGE_PENALIZATION, 0.0 ) # alpha_f
        # self.model_part.Properties[1].SetValue(LOCAL_NONLOCAL_PLASTICITY_PENALIZATION, 0.0 ) # alpha_p
        # self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, GTN3dImplicit() )
        ### end of enable this
        self.model_part.Properties[1].SetValue(MATERIAL_PARAMETERS, mat_params )
        hpoints = [375.00, 0.00,
            427.42, 0.02,
            458.41, 0.04,
            484.44, 0.06,
            507.70, 0.08,
            529.10, 0.10,
            549.13, 0.12,
            568.07, 0.14,
            586.14, 0.16,
            603.48, 0.18,
            620.19, 0.20,
            636.36, 0.22,
            652.05, 0.24,
            667.31, 0.26,
            682.19, 0.28,
            696.73, 0.30,
            710.94, 0.32,
            724.87, 0.34,
            738.53, 0.36,
            751.94, 0.38,
            765.12, 0.40,
            778.08, 0.42,
            790.84, 0.44,
            803.41, 0.46,
            815.80, 0.48,
            828.03, 0.50,
            840.09, 0.52,
            852.00, 0.54,
            863.77, 0.56,
            875.39, 0.58,
            886.89, 0.60,
            898.26, 0.62,
            909.51, 0.64,
            920.64, 0.66,
            931.67, 0.68,
            942.58, 0.70,
            953.40, 0.72,
            964.11, 0.74,
            974.74, 0.76,
            985.26, 0.78,
            995.70, 0.80,
            1006.06, 0.82,
            1016.33, 0.84,
            1026.52, 0.86,
            1036.63, 0.88,
            1046.67, 0.90,
            1056.64, 0.92,
            1066.53, 0.94,
            1076.35, 0.96,
            1086.11, 0.98,
            1095.80, 1.00,
            1105.43, 1.02,
            1114.99, 1.04,
            1124.50, 1.06,
            1133.94, 1.08,
            1143.33, 1.10,
            1152.66, 1.12,
            1161.94, 1.14,
            1171.16, 1.16,
            1180.33, 1.18,
            1189.45, 1.20,
            1198.52, 1.22,
            1207.54, 1.24,
            1216.52, 1.26,
            1225.44, 1.28,
            1234.32, 1.30,
            1243.16, 1.32,
            1251.95, 1.34,
            1260.70, 1.36,
            1269.40, 1.38,
            1278.07, 1.40,
            1286.69, 1.42,
            1295.28, 1.44,
            1303.82, 1.46,
            1312.32, 1.48,
            1320.79, 1.50,
            1329.22, 1.52,
            1337.62, 1.54,
            1345.98, 1.56,
            1354.30, 1.58,
            1362.59, 1.60,
            1370.84, 1.62,
            1379.06, 1.64,
            1387.25, 1.66,
            1395.40, 1.68,
            1403.53, 1.70,
            1411.62, 1.72,
            1419.68, 1.74,
            1427.71, 1.76,
            1435.71, 1.78,
            1443.68, 1.80,
            1451.62, 1.82,
            1459.53, 1.84,
            1467.42, 1.86,
            1475.27, 1.88,
            1483.10, 1.90,
            1490.90, 1.92,
            1498.68, 1.94,
            1506.43, 1.96,
            3746.00, 10.00]
        npoints = len(hpoints)/2
        hlaw = PiecewiseLinearHardeningLaw()
        for i in range(0, npoints):
            hlaw.AddPoint(hpoints[2*i+1], hpoints[2*i])
        hlaw.Assign(ISOTROPIC_HARDENING_LAW, hlaw, self.model_part.Properties[1])
        print("GTN (perfectly plastic) model selected")
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
