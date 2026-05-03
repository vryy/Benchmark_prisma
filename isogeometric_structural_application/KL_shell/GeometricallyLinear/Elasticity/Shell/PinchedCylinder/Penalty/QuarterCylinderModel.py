from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
import numpy
import matplotlib
from pylab import *

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
CheckForPreviousImport()

class SolutionScheme:
    def __init__(self, problem_name, path, logging=True, nsampling=40, order=2):
        self.path = path
        self.problem_name = problem_name

        ##
        self.params = {}
        self.params['name'] = self.problem_name
        self.params['division mode'] = "non-uniform"
        self.params['division number u'] = 20
        self.params['division number v'] = 20

        if 'name' not in self.params:
            self.params['name'] = "iga model"
        if 'base element name' not in self.params:
            self.params['base element name'] = "DummyElement"
        if 'base condition name' not in self.params:
            self.params['base condition name'] = "DummyLineCondition"
        if 'last node id' not in self.params:
            self.params['last node id'] = 1
        if 'last element id' not in self.params:
            self.params['last element id'] = 1
        if 'last condition id' not in self.params:
            self.params['last condition id'] = 1
        if 'division mode' not in self.params:
            self.params['division mode'] = "uniform"
        # default variables list
        if 'variables list' not in self.params:
            self.params['variables list'] = []
            self.params['variables list'].append(DISPLACEMENT)
            self.params['variables list'].append(REACTION)

        # create geometry
        # note that model_part and model_part_post are declared in CreateGeometry()
        self.CreateFEModel(nsampling=nsampling, order=order)

        self.solver = self.CreateSolver2(logging=logging, abs_tol=1e-6, rel_tol=1e-10)

    def CreateSolver(self, abs_tol = 1e-06, rel_tol = 1e-10):
        # defining linear solver
        if KratosMKLSolversApplication.Has("MKLPardisoSolver"):
            structure_linear_solver = MKLPardisoSolver()
        else:
            structure_linear_solver = SuperLUSolver()

        # defining builder_and_solver
        builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(structure_linear_solver)

        # defining time scheme
        time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # defining convergence criteria
        conv_criteria = DisplacementCriteria(rel_tol, abs_tol)

        # defining solving strategy flags
        MaxNewtonRapshonIterations = 30
        CalculateReactionFlag = True
        ReformDofSetAtEachStep = True
        MoveMeshFlag = False

        # defining solving strategy
        parallel_space = UblasSparseSpace()

        import shell_analysis_strategy
        solver = shell_analysis_strategy.ShellAnalysisSolver( parallel_space, \
                                                              None, \
                                                              self.model_part, \
                                                              time_scheme, \
                                                              structure_linear_solver, \
                                                              conv_criteria, \
                                                              builder_and_solver, \
                                                              MaxNewtonRapshonIterations, \
                                                              CalculateReactionFlag, \
                                                              ReformDofSetAtEachStep, \
                                                              MoveMeshFlag, \
                                                              None )

        return solver

    def CreateSolver2(self, logging=True, abs_tol = 1e-06, rel_tol = 1e-10):
        analysis_parameters = {}
        analysis_parameters['dimension'] = 2
        analysis_parameters['perform_contact_analysis_flag'] = False
        analysis_parameters['print_sparsity_info_flag'] = False
        analysis_parameters['analysis_type'] = 0
        analysis_parameters['dissipation_radius'] = 0.1
        analysis_parameters['decouple_build_and_solve'] = True
        analysis_parameters['solving_scheme'] = 'monolithic'
        analysis_parameters['residuum_tolerance'] = 5e-3
        analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
        analysis_parameters['convergence_criteria'] = "displacement"
        analysis_parameters['list_dof'] = True
        analysis_parameters['log_residuum'] = logging

        import structural_solver_advanced
        solver = structural_solver_advanced.SolverAdvanced( self.model_part, 2, 1, analysis_parameters, abs_tol, rel_tol )
        self.AddVariables( self.model_part )
        self.AddDofs(self.model_part)
        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        if KratosMKLSolversApplication.Has("MKLPardisoSolver"):
            plinear_solver = MKLPardisoSolver()
        else:
            plinear_solver = SuperLUSolver()
        solver.structure_linear_solver = plinear_solver
        solver.Initialize()
        solver.solver.SetEchoLevel(2)
        solver.solver.max_iter = 20 #control the maximum iterations of Newton Raphson loop
        solver.solver.MoveMeshFlag = False

        return solver

    def Solve( self, time ):
        self.model_part.CloneTimeStep(time)
        return self.solver.Solve()

    def InitializeModel( self, thickness=3.0 ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStress() )
        self.model_part.Properties[1].SetValue(THICKNESS, thickness)
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS, 3.0e6)
        self.model_part.Properties[1].SetValue(POISSON_RATIO, 0.3)

        #self.model_part.Check(self.model_part.ProcessInfo)

    def AddVariables(self,model_part):
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part.AddNodalSolutionStepVariable(FORCE)

    def AddDofs(self,model_part):
        for node in model_part.Nodes:
             node.AddDof(DISPLACEMENT_X, REACTION_X)
             node.AddDof(DISPLACEMENT_Y, REACTION_Y)
             node.AddDof(DISPLACEMENT_Z, REACTION_Z)
             node.AddDof(ROTATION_X)
             node.AddDof(ROTATION_Y)
             node.AddDof(ROTATION_Z)

    ################################
    ######################################################
    def CreateControlNet(self):

        nurbs_fespace_library = BSplinesFESpaceLibrary()
        grid_lib = ControlGridLibrary()
        multipatch_util = MultiPatchUtility()
        multipatch_refine_util = MultiPatchRefinementUtility()
        bsplines_patch_util = BSplinesPatchUtility()
        mpatch_export = MultiNURBSPatchMatlabExporter()

        import geometry_factory

        #############################################
        ############### patch 1 #####################
        ####### create arc 1
        R = 300.0
        L = 600.0

        arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'y', R, 180.0, 90.0)
        arc1 = arc1_ptr.GetReference()

        ####### create arc 2
        arc2_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R, 180.0, 90.0)
        arc2 = arc2_ptr.GetReference()

        # create patch 1
        patch1_ptr = bsplines_patch_util.CreateLoftPatch(arc1, arc2)
        patch1 = patch1_ptr.GetReference()
        patch1.LayerIndex = 1
        patch1.Id = 1

        ######create multipatch
        mpatch = MultiPatch2D()
        mpatch.AddPatch(patch1_ptr)

        return mpatch

    def CreateFEModel(self, nsampling = 40, order = 2):

        multipatch_util = MultiPatchUtility()
        multipatch_refine_util = MultiPatchRefinementUtility()
        bending_strip_utility = BendingStripUtility()

        # initial control net
        self.mpatch = self.CreateControlNet()

        print("############REFINEMENT###############")
        multipatch_refine_util = MultiPatchRefinementUtility()
        knot_u = []
        for i in range(1, nsampling):
            knot_u.append(float(i) / nsampling)
        knot_v = knot_u
        if order >= 2:
            order_u = order - 2
            order_v = order - 1
        else:
            raise Exception(f"Invalid order {order}")
        multipatch_refine_util.DegreeElevate(self.mpatch[1], [order_u, order_v])
        multipatch_refine_util.InsertKnots(self.mpatch[1], [knot_u, knot_v])

        # print(self.mpatch)
        # sys.exit()

        # create model part
        self.mpatch_mp= MultiPatchModelPart2D(self.mpatch)

        # begin model part
        self.mpatch_mp.BeginModelPart()

        self.model_part = self.mpatch_mp.GetModelPart()
        self.post_model_part = ModelPart("iga-fem mesh ")

        # define variables
        self.AddVariables(self.model_part)
        self.AddVariables(self.post_model_part)
        self.model_part.Properties[1].SetValue(NUM_IGA_INTEGRATION_METHOD, 3)

        # create nodes
        self.mpatch_mp.CreateNodes()

        # add elements
        element_name = "KirchhoffLoveLinearShellBezier2D3"
        last_elem_id = multipatch_util.GetLastElementId(self.model_part)
        prop = self.model_part.Properties[1]
        elems1 = self.mpatch_mp.AddElements(self.mpatch[1], element_name, last_elem_id+1, prop)

        #sys.exit()
        #
        # end model part
        self.mpatch_mp.EndModelPart()

        # add dofs
        self.AddDofs( self.model_part )
        self.AddDofs( self.post_model_part)

    def WriteGiD(self, time):
        #######WRITE TO GID
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        gid_io = StructuralGidIO(self.params['name'], post_mode, multi_file_flag, write_deformed_flag, write_elements)
        gid_io.InitializeMesh( time )
        post_mesh = self.post_model_part.GetMesh()
        gid_io.WriteMesh( post_mesh )
        print("mesh written...")
        gid_io.FinalizeMesh()
        gid_io.InitializeResults( time, post_mesh )
        print("write nodal results")
        for var in self.params['variables list']:
            gid_io.WriteNodalResults(var, self.post_model_part.Nodes, time, 0)
        gid_io.FinalizeResults()

    def WriteOutput( self,  time ):

        fem_mesh = NonConformingMultipatchLagrangeMesh2D(self.mpatch)

        fem_mesh.SetBaseElementName(self.params['base element name'])
        fem_mesh.SetBaseConditionName(self.params['base condition name'])
        fem_mesh.SetLastNodeId(self.params['last node id'])
        fem_mesh.SetLastElemId(self.params['last element id'])
        fem_mesh.SetLastCondId(self.params['last condition id'])
        if self.params['division mode'] == "uniform":
            fem_mesh.SetUniformDivision(self.params['uniform division number'])
        elif self.params['division mode'] == "non-uniform":
            for patch_ptr in self.mpatch.Patches():
                patch = patch_ptr.GetReference()
                fem_mesh.SetDivision(patch.Id, 0, self.params['division number u'])
                fem_mesh.SetDivision(patch.Id, 1, self.params['division number v'])

        self.mpatch_mp.SynchronizeBackward(DISPLACEMENT)
        self.mpatch_mp.SynchronizeBackward(REACTION)
        fem_mesh.WriteModelPart(self.post_model_part)
        self.WriteGiD( time)

    def FindInnerCoordinate(self, side= "u1"):
        patch1_ptr = self.mpatch[1]
        patch1 = patch1_ptr.GetReference()

        if side == "u1":
            bpatch_ptr = patch1.ConstructBoundaryPatch(BoundarySide2D.U1)
            bpatch = bpatch_ptr.GetReference()
            cg = bpatch.GridFunction(CONTROL_POINT_COORDINATES).ControlGrid
            n = cg.size()
            return cg[n-2][1]
        elif side == "v1":
            bpatch_ptr = patch1.ConstructBoundaryPatch(BoundarySide2D.V1)
            bpatch = bpatch_ptr.GetReference()
            cg = bpatch.GridFunction(CONTROL_POINT_COORDINATES).ControlGrid
            return cg[1][2]
