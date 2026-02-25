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
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.DistributedBuildersApplication import *
from KratosMultiphysics.TrilinosSolversApplication import *
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path ):
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("Kratos_simulation")
        self.path = path
        self.problem_name = problem_name
        ##################################################################
        ## ADD VARIABLES #################################################
        ##################################################################
        import structural_solver_static as this_solver_strategy
        this_solver_strategy.AddVariables(self.model_part)
        set_communicator_process = SetMPICommunicatorProcess(self.model_part)
        set_communicator_process.Execute() #set the communicator of the model_part to MPICommunicator. This is important for the model_part to read the partition correctly
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.read_partition_file = True
        number_of_partitions = mpi.size
        if number_of_partitions > 1:
            if self.read_partition_file == True:
                if mpi.rank == 0
                    dimension = 3
                    verbose = 0
                    synchronize_conds = True
                    model_part_io = ModelPartIO(self.path+self.problem_name)
                    partitioning_process = MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions, dimension, verbose, synchronize_conds)
                    partitioning_process.Execute()
                mpi.world.barrier()
            self.model_part_io = ModelPartIO(self.path+self.problem_name+'_'+str(mpi.rank))
            self.model_part_io.ReadModelPart(self.model_part)
        else:
            self.model_part_io = ModelPartIO(self.path+self.problem_name)
            self.model_part_io.ReadModelPart(self.model_part)
        mpi.world.barrier()
        self.meshWritten = False
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
        self.cond_file.close()

        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        this_solver_strategy.AddDofs(self.model_part)
        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        # define space and parallel communicator
        self.parallel_space = TrilinosEpetraSpace()
        self.comm = self.parallel_space.CreateCommunicator(mpi.world)

        # defining linear solver
        self.structure_linear_solver = TrilinosEpetraSolver()
        # Please add your preferred solver for Trilinos Epetra backend

        # defining builder_and_solver
        self.builder_and_solver = TrilinosEpetraResidualBasedEliminationBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)

        # defining time scheme
        self.time_scheme = TrilinosEpetraResidualBasedIncrementalUpdateStaticScheme()

        # defining convergence criteria
        self.conv_criteria = TrilinosEpetraDisplacementCriteria(       1e-06,        1e-09)

        # defining solving strategy flags
        self.MaxNewtonRapshonIterations = 30
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True
        self.EchoLevel = 0b0000000000011100

        # defining solving strategy
        self.solver.SetEchoLevel(self.EchoLevel)

    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False

        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,            1 )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.3 )
        self.model_part.Properties[1].SetValue(THICKNESS,            1 )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStress() )

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
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        if(mpi.rank == 0):
            print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( meshname, mesh )
        if(mpi.rank == 0):
            print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
        self.gid_io.WriteNodalResults(PARTITION_INDEX, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
        self.gid_io.FinalizeResults()
        if(mpi.rank == 0):
            mergefile = open( self.path+self.problem_name+"_"+str(time)+"_merge_results.bch", 'w' )
            mergefile.write("Postprocess\n")
            mergefile.write("mescape\n \n")
            meshname = time+mpi_fn_step
            mergefile.write("Files Read "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
            mergefile.write("mescape\n")
            for rank in range(2, mpi.size+1 ):
                meshname = time+mpi_fn_step*(rank)
                mergefile.write("Files Add "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
                mergefile.write("mescape\n")
            mergefile.write("Files SaveAll BinMeshesSets\n")
            mergefile.write(self.path+"merged_"+self.problem_name+"_"+ str(time)+".bin\n")
            mergefile.write("mescape\n")
            mergefile.write("\n")
            mergefile.close()
        mpi.world.barrier()

    def WriteRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time) + "_" + str(mpi.rank)
        serializer = Serializer(fn)
        serializer.Save("ModelPart", self.model_part)
        serializer = 0
        if mpi.rank == 0:
            print("Write restart data to " + self.problem_name + "_" + str(time) + " completed")

    def LoadRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time) + "_" + str(mpi.rank)
        serializer = Serializer(fn)
        serializer.Load("ModelPart", self.model_part)
        serializer = 0 
        if mpi.rank == 0:
            print("Load restart data from " + self.problem_name + "_" + str(time) + " completed")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
##################################################################
