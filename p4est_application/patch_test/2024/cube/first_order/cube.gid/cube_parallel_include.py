##################################################################
################### distributed_include.py   #####################
##################################################################
##### KRATOS Multiphysics                                    #####
##### include file for distributed-memory simulation         #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.DistributedBuildersApplication import *
# from KratosMultiphysics.TrilinosSolversApplication import *
from KratosMultiphysics.PetscSolversApplication import *
from KratosMultiphysics.P4estApplication import *
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path, results_path ):
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("Kratos_simulation")
        self.path = path
        self.results_path = results_path
        self.problem_name = problem_name
        ##################################################################
        ## ADD VARIABLES #################################################
        ##################################################################
        ## generating solver
        import structural_solver_advanced
        structural_solver_advanced.AddVariables( self.model_part )
        self.model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        # post_mode = GiDPostMode.GiD_PostAscii
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.results_path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
        if( mpi.rank == 0 ):
            self.mergefile = open( self.path+self.problem_name+"_merge_results.bch", 'w' )
            self.mergefile.write("Postprocess\n")
            self.mergefile.write("mescape\n \n")
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
        print "input data read OK"
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

        ##################################################################
        ## PARALLEL SPACE AND COMMUNICATOR ###############################
        ##################################################################
        # define space and parallel communicator
        self.parallel_space = PETScSpace()
        self.comm = self.parallel_space.CreateCommunicator(mpi.world)

        # defining linear solver
        self.structure_linear_solver = PETScSolver()

        # initialize first time the solver
        self.SetModelPart(self.model_part)

    def SetModelPart(self, model_part):
        print("SetModelPart is called")
        self.model_part = model_part

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        import structural_solver_advanced
        structural_solver_advanced.AddDofs( self.model_part )
        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################

        # defining builder_and_solver
        # self.builder_and_solver = PETScResidualBasedEliminationBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
        # self.builder_and_solver = PETScResidualBasedBlockBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
        self.builder_and_solver = PETScResidualBasedBlockBuilderAndSolverDeactivationWithConstraints(self.comm, self.structure_linear_solver)
        #self.builder_and_solver = PETScResidualBasedBlockBuilderAndSolverDeactivationWithConstraintsElementWise(self.comm, self.structure_linear_solver)

        # defining time scheme
        self.time_scheme = PETScResidualBasedIncrementalUpdateStaticScheme()

        # defining convergence criteria
        self.conv_criteria = PETScMultiphaseFlowCriteria(       1e-06,        1e-09)

        # defining solving strategy flags
        self.MaxNewtonRapshonIterations = 30
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True
#        self.EchoLevel = 0b0000000000011100
        self.EchoLevel = 0x0001 + 0x0020 + 0x0040 + 0x0080 + 0x0008 + 0x0010 + 0x0004

        # defining solving strategy
        import distributed_strategies
        self.solver = distributed_strategies.ResidualBasedNewtonRaphsonStrategy(self.parallel_space, \
                        self.comm, \
                        self.model_part, \
                        self.time_scheme, \
                        self.structure_linear_solver, \
                        self.conv_criteria, \
                        self.builder_and_solver, \
                        self.MaxNewtonRapshonIterations, \
                        self.CalculateReactionFlag, \
                        self.ReformDofSetAtEachStep, \
                        self.MoveMeshFlag)
        self.solver.SetEchoLevel(self.EchoLevel)

    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        self.model_part.Properties[1].SetValue(DENSITY,         0.0 )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,      2.0 )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.3 )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
        print "Linear elastic model selected"

        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        model_layers = __import__(self.problem_name+"_layers")
        ## ELEMENTS on layers ############################################
        self.layer_sets = model_layers.ReadLayerSets()
        ## NODES on layers ###############################################
        self.layer_nodes_sets = model_layers.ReadLayerNodesSets()
        ##################################################################
        print "layer sets stored"
        ##################################################################
        ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        ##################################################################
        self.node_groups = model_layers.ReadNodeGroups()
        print "node groups stored"
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
        self.solver.Initialize()
        if(mpi.rank == 0):
            print "activation utility initialized"
        self.model_part.Check( self.model_part.ProcessInfo )
        if(mpi.rank == 0):
            print "model successfully initialized"

    def WriteOutput( self, time ):
        mpi_fn_step = 0.0001
        meshname = time+mpi_fn_step*(mpi.rank+1)
        self.gid_io.InitializeMesh( meshname )
        mesh = self.model_part.GetMesh()
        # mesh = self.model_part.GetCommunicator().LocalMesh()
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        if(mpi.rank == 0):
            print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( meshname, mesh )
        if(mpi.rank == 0):
            print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
        self.gid_io.WriteNodalResults(PARTITION_INDEX, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        self.gid_io.FinalizeResults()
        if(mpi.rank == 0):
            meshname = time+mpi_fn_step
            self.mergefile.write("Files Read "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
            self.mergefile.write("mescape\n")
            for rank in range(2, mpi.size+1 ):
                meshname = time+mpi_fn_step*(rank)
                self.mergefile.write("Files Add "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
                self.mergefile.write("mescape\n")
            self.mergefile.write("Files SaveAll BinMeshesSets\n")
            self.mergefile.write(self.path+"merged_"+self.problem_name+"_"+ str(time)+".bin\n")
            self.mergefile.write("mescape\n")
            self.mergefile.write("\n")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()
        if(mpi.rank == 0):
            self.mergefile.close()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
##################################################################
