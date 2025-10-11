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
import math
import time as time_module
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
# from KratosMultiphysics.ExternalSolversApplication import *
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
        self.model_part = ComplexModelPart(problem_name)
        self.path = path
        self.results_path = results_path+os.sep
        self.problem_name = problem_name
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################

        freq = 10.0 # Hz

        self.abs_tol =        1e-10
        self.rel_tol =        1e-10

        ## generating solver
        import structural_solver_static
        self.solver = structural_solver_static.StaticStructuralSolver( self.model_part, self.domain_size, abs_tol=self.abs_tol, rel_tol=self.rel_tol )
        self.solver.time_scheme = ComplexSteadyStateDynamicsScheme(2*math.pi*freq)
        self.solver.CalculateReactionFlag = True
        self.solver.ReformDofSetAtEachStep = True
        self.solver.MoveMeshFlag = False
        self.AddVariables( self.model_part )
        #ekate_solver_parallel.AddVariables( self.model_part )
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        self.model_part_io = ComplexModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
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
        self.gid_io = ComplexSDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = MKLComplexPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)
        (self.solver.solver).max_iter = 10 #control the maximum iterations of Newton Raphson loop

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )

    def AddVariables(self, model_part):
        model_part.AddNodalSolutionStepVariable(COMPLEX_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(COMPLEX_DISPLACEMENT_DT)
        model_part.AddNodalSolutionStepVariable(COMPLEX_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(COMPLEX_VELOCITY)
        model_part.AddNodalSolutionStepVariable(COMPLEX_REACTION)

    def AddDofsForNodes(self, nodes):
        for node in nodes:
            node.AddDof(COMPLEX_DISPLACEMENT_X, COMPLEX_REACTION_X)
            node.AddDof(COMPLEX_DISPLACEMENT_Y, COMPLEX_REACTION_Y)
            node.AddDof(COMPLEX_DISPLACEMENT_Z, COMPLEX_REACTION_Z)

    def AddDofs(self, model_part):
        for node in model_part.Nodes:
            node.AddDof(COMPLEX_DISPLACEMENT_X, COMPLEX_REACTION_X)
            node.AddDof(COMPLEX_DISPLACEMENT_Y, COMPLEX_REACTION_Y)
            node.AddDof(COMPLEX_DISPLACEMENT_Z, COMPLEX_REACTION_Z)

    def SetModelPart(self, model_part):
        self.model_part = model_part
        number_of_time_steps = 1

        ## generating solver
        self.solver = structural_solver_static.StaticStructuralSolver( self.model_part, self.domain_size, abs_tol=self.abs_tol, rel_tol=self.rel_tol )
        self.solver.ReformDofSetAtEachStep = True
        self.solver.time_scheme = ComplexSteadyStateDynamicsScheme(2*math.pi*10)
        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        self.AddDofs(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = MKLComplexPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2)

    def SetOutputPath(self, path):
        self.results_path = path
        self.gid_io = ComplexSDGidPostIO( self.results_path+self.problem_name, self.post_mode, self.multi_file_flag, self.write_deformed_flag, self.write_elements )

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

    def ApplyInsituWaterPressure2D( self, free_node_list_water, free_node_list_air, y_zero, gravity_y):
        water_density=1000.0;
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
        self.gid_io.WriteNodalResults(COMPLEX_DISPLACEMENT, time, 0)
        print("nodal COMPLEX_DISPLACEMENT written")
        self.gid_io.WriteNodalResults(COMPLEX_DISPLACEMENT_DT, time, 0)
        print("nodal COMPLEX_DISPLACEMENT_DT written")
        self.gid_io.WriteNodalResults(COMPLEX_ACCELERATION, time, 0)
        print("nodal COMPLEX_ACCELERATION written")
        self.gid_io.WriteNodalResults(COMPLEX_REACTION, time, 0)
        print("nodal COMPLEX_REACTION written")
        # self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        # print("gauss point STRESSES written")
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
        end_time = time_module.time()
        print("conditions in node groups stored, time = " + str(end_time - start_time) + "s")
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        self.model_part.Properties[1].SetValue(DENSITY,         2000 )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,     1.91e+08 )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.3 )
        self.model_part.Properties[1].SetValue(THICKNESS, 1.0 )
        self.model_part.Properties[1].SetValue(COMPLEX_CONSTITUTIVE_LAW, ComplexIsotropic3D() )
        print("Linear elastic model selected, description: ComplexIsotropic3D")
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
