import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ThermalApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.IsogeometricThermalApplication import *
# from KratosMultiphysics.ExternalSolversApplication import *
# from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.MultigridSolversApplication import *

kernel = Kernel()   #defining kernel

import model_iga_include
from model_iga_include import *

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

import geometry_factory

def CreateMultiPatch():

    #### create rectagle
    a = 1.0
    rec_ptr = geometry_factory.CreateRectangle([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
    rec = rec_ptr.GetReference()
    rec.Id = 1
    rec.LayerIndex = 1

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(rec_ptr)

    multipatch_refine_util.DegreeElevate(mpatch[1], [2, 2])

    return mpatch

def Refine(mpatch, ins_knot):
    print("###############REFINEMENT###############")

    trans_list = {}
    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        trans_list[patch.Id] = []

    trans = multipatch_refine_util.InsertKnotsGetTrans(mpatch[1], [ins_knot, ins_knot])
    for pid, mat in trans.items():
        trans_list[pid].append(mat)
#    for pid, mat in trans1.items():
#        print("pid1:", pid)
#        print("mat1:", str(mat))

    trans_tot = {}
    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        is_empty = True
        for mat in trans_list[patch.Id]:
            if is_empty:
                trans_tot[patch.Id] = mat
                is_empty = False
            else:
                trans_tot[patch.Id] = mat*trans_tot[patch.Id]

#    for pid, mat in trans_tot.items():
#        print("pid:", pid)
#        print("mat_tot:", str(mat))

    return [mpatch, trans_tot]

def CreateModelPart(mpatch):
    mpatch_util = MultiPatchUtility()
    element_name = "LinearPoissonBezier2D"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)

    mpatch_mp.CreateNodes()

    #problem data
    prop = model_part.Properties[1]
    prop.SetValue(THERMAL_CONDUCTIVITY, 1.0 )
    prop.SetValue(THICKNESS, 1)
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 1)

    patch_ids = [1]
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)

    mpatch_mp.EndModelPart()
    #    print(mpatch_mp)

    # fix temperature on left
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0.0)

    # fix temperature on right
    for node in model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0.0)

    return mpatch_mp

def CreateModels(logging=True):

    #refine_list = [[0.5], [0.25, 0.75], [0.125, 0.375, 0.625, 0.875]]
    refine_list = [[0.5], [0.25, 0.75]]
#    refine_list = [[0.5]]
    model_list = []
    mpatch_mp_list = []
    trans_mat_list = []

    for i in range(0, len(refine_list)+1):
        mpatch = CreateMultiPatch()

        for j in range(0, i):
            [mpatch, trans_tot] = Refine(mpatch, refine_list[j])
#            print(trans_tot)
            if j == i-1:
#                print(trans_tot)
                trans_mat_list.append(trans_tot)

        mpatch.Enumerate()
#        trans_mat_list.append(trans_tot)

    #    print(mpatch)
        mpatch_export.Export(mpatch, "square_" + str(i) + ".m")

        mpatch_mp = CreateModelPart(mpatch)
        mpatch_mp_list.append(mpatch_mp)

        #############ANALYSIS MODEL#######################################
        model_part = mpatch_mp.GetModelPart()
        params = model_iga_include.StaticParameters()
        params["builder_and_solver_type"] = "residual-based block"
        params["log_residuum"] = logging
        model = model_iga_include.Model('linear_poisson_square', os.getcwd()+"/", model_part, params)
        model.solver.solver.convergence_criteria = DisplacementCriteria(model.rel_tol, model.abs_tol)
        for node in model.model_part.Nodes:
            node.AddDof(TEMPERATURE)
        model.InitializeModel()

        # assign temperature on right
        tol = 1.0e-6
        for node in model_part.Nodes:
            if abs(node.X0 - 1.0) < tol:
                node.SetSolutionStepValue(TEMPERATURE, 10.0)

        model_list.append(model)

    # here we make the first model as the finest one
    mpatch_mp_list.reverse()
    model_list.reverse()
    trans_mat_list.reverse()
    return [mpatch_mp_list, model_list, trans_mat_list]

def main(logging=True, output=True):
    #############ANALYSIS MODEL#######################################
    [mpatch_mp_list, model_list, trans_mat_list] = CreateModels(logging=logging)

    model = model_list[0]
    mpatch_mp = mpatch_mp_list[0]
    mpatch = mpatch_mp.GetMultiPatch()

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)

    if output:
        ######Synchronize back the results to multipatch
        mpatch_mp.SynchronizeBackward(TEMPERATURE)
        ##################################################################

        ## post processing
        params_post = {}
        params_post['name'] = "linear_poisson_square"
        params_post['division mode'] = "uniform"
        params_post['uniform division number'] = 20
    #    params_post['division mode'] = "non-uniform"
    #    params_post['division number u'] = 10
    #    params_post['division number v'] = 10
    #    params_post['division number w'] = 1
        params_post['variables list'] = [TEMPERATURE]
        dim = 2
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    return model

def main2(logging=True, output=True):
    #############ANALYSIS MODEL#######################################
    [mpatch_mp_list, model_list, trans_mat_list] = CreateModels(logging=logging)

    model0 = model_list[0]
    mpatch_mp0 = mpatch_mp_list[0]

    print("len(model_list):", len(model_list))
    for i in range(0, len(model_list)):
        print(" nnodes model " + str(i) + ": " + str(model_list[i].model_part.NumberOfNodes()))
#    print("len(trans_mat_list):", len(trans_mat_list))
#    for trans_mat in trans_mat_list:
# #       print(trans_mat)
#        for pid, mat in trans_mat.items():
#            print("pid:", pid)
#            print("trans_mat:", str(mat))
    # sys.exit(0)

    ##################################################################
    ####SETUP GMG SOLVER##############################################
    ##################################################################

    nlevels = len(model_list)
    block_size = 1

    ## construct multigrid solver for first level
    #defining linear solver
    plinear_solver = MultilevelSolver()

    plinear_solver.SetCycle("V")
    plinear_solver.AddPreSmoother(JacobiIterativeSolver(3, 0.66))
    plinear_solver.AddPostSmoother(JacobiIterativeSolver(3, 0.66))
    plinear_solver.SetCoarseSolver(SkylineLUFactorizationSolver())
    plinear_solver.SetEchoLevel(2)

    solver_plist = MGParameterList()
    solver_factory = MultilevelSolverFactory(solver_plist) # this is an empty factory
    solver_factory.SetMute(True)
    plinear_solver.SetFactory(solver_factory)

    for lvl in range(0, nlevels):
        plevel = MatrixBasedMGLevel(lvl)

        if lvl < nlevels-1:
            nnodes_fine = model_list[lvl].model_part.NumberOfNodes()
            nnodes_coarse = model_list[lvl+1].model_part.NumberOfNodes()
            print("At level", lvl, "nnodes_coarse:", nnodes_coarse, ", nnodes_fine:", nnodes_fine)
#            Prolongator = IndexBasedMGProjector(nnodes_fine, nnodes_coarse)
#            Prolongator.SetStride(block_size)
#            Prolongator.Initialize()
#            for pid, mat in trans_mat_list[lvl].items():
#                fine_patch = mpatch_mp_list[lvl].GetMultiPatch().Patches()[pid].GetReference()
#                coarse_patch = mpatch_mp_list[lvl+1].GetMultiPatch().Patches()[pid].GetReference()
#                coarse_indices = coarse_patch.FESpace().FunctionIndices()
#                fine_indices = fine_patch.FESpace().FunctionIndices()
#                print("coarse_indices:", coarse_indices)
#                print("fine_indices:", fine_indices)
#                print("mat:", str(mat))
#                Prolongator.AssembleOperator(fine_indices, coarse_indices, mat)
#            print("assembled Prolongator for level", lvl)

            Prolongator = MatrixBasedMGProjector(nnodes_fine*block_size, nnodes_coarse*block_size)
            for pid, mat in trans_mat_list[lvl].items():
                fine_patch = mpatch_mp_list[lvl].GetMultiPatch().Patches()[pid].GetReference()
                coarse_patch = mpatch_mp_list[lvl+1].GetMultiPatch().Patches()[pid].GetReference()
                coarse_indices = coarse_patch.FESpace().FunctionIndices()
                fine_indices = fine_patch.FESpace().FunctionIndices()
                print("coarse_indices:", coarse_indices)
                print("fine_indices:", fine_indices)
                print("mat:", str(mat))
                Prolongator.AssembleOperator(fine_indices, coarse_indices, mat, block_size)
            print("assembled Prolongator for level", lvl)

            #################################
#            print("make a simple test for the prolongator")
#            pX = Vector(nnodes_coarse*block_size)
#            pX_0 = Vector(nnodes_coarse)
#            pX_1 = Vector(nnodes_coarse)
#            for node in model_list[lvl+1].model_part.Nodes:
#                pX[2*(node.Id-1)] = node.X0
#                pX_0[node.Id-1] = node.X0
#                pX[2*(node.Id-1)+1] = node.Y0
#                pX_1[node.Id-1] = node.Y0

#            pYref = Vector(nnodes_fine*block_size)
#            for node in model_list[lvl].model_part.Nodes:
#                pYref[2*(node.Id-1)] = node.X0
#                pYref[2*(node.Id-1)+1] = node.Y0

#            pY = Vector(nnodes_fine*block_size)
#            Prolongator.Apply(pX, pY)

##            pY_0 = Vector(nnodes_fine)
##            pY_0 = trans_mat_list[lvl][1]*pX_0
##            pY_1 = trans_mat_list[lvl][1]*pX_1
##            pY2 = Vector(nnodes_fine*block_size)
##            for i in range(0, nnodes_fine):
##                pY2[2*i] = pY_0[i]
##                pY2[2*i+1] = pY_1[i]

#            print("pX:", str(pX))
#            print("pY:", str(pY))
##            print("pY2:", str(pY2))
#            print("pYref:", str(pYref))
#            sys.exit(0)
            ##################################

            Restrictor = MGTransposeProjector(Prolongator)
            print("Prolongator.GetBaseSize():", Prolongator.GetBaseSize())
            print("Prolongator.GetProjectedSize():", Prolongator.GetProjectedSize())
            print("Restrictor.GetBaseSize():", Restrictor.GetBaseSize())
            print("Restrictor.GetProjectedSize():", Restrictor.GetProjectedSize())
        else:
            Prolongator = MGNullProjector()
            Restrictor = MGNullProjector()

        plevel.SetProlongationOperator(Prolongator)
        plevel.SetRestrictionOperator(Restrictor)
        plinear_solver.AddLevel(plevel)

#    sys.exit(0)

    model0.solver.structure_linear_solver = plinear_solver
    model0.solver.Initialize()
    model0.solver.solver.SetEchoLevel(2)
    model0.solver.solver.max_iter = 10 #control the maximum iterations of Newton Raphson loop
    model0.solver.solver.MoveMeshFlag = False
    model0.solver.solver.convergence_criteria = DisplacementCriteria(model0.rel_tol, model0.abs_tol)
    model0.solver.solver.builder_and_solver = ResidualBasedBlockBuilderAndSolver(plinear_solver)

    ##################################################################
    ###################SOLVE##########################################
    ##################################################################

    level_list = []
    for i in range(0, nlevels):
        level_list.append(plinear_solver.GetLevel(i))

    params = {}
    params['print_sparsity_info_flag'] = False
    params['stop_Newton_Raphson_if_not_converge'] = True

    import gmg_newton_raphson_strategy
    gmg_solve_strategy = gmg_newton_raphson_strategy.SolvingStrategyPython(model_list, level_list, params)

    ## solve first step
    time = 1.0
    gmg_solve_strategy.Solve(time, 0, 0, 0, 0)

    if output:
        ######Synchronize back the results to multipatch
        for mpatch_mp in mpatch_mp_list:
            mpatch_mp.SynchronizeBackward(TEMPERATURE)
        ##################################################################

        ## post processing
        params_post = {}
        params_post['division mode'] = "uniform"
        params_post['uniform division number'] = 20
    #    params_post['division mode'] = "non-uniform"
    #    params_post['division number u'] = 10
    #    params_post['division number v'] = 10
    #    params_post['division number w'] = 1
        params_post['variables list'] = [TEMPERATURE]
        dim = 2

        for i in range(0, len(mpatch_mp_list)):
            params_post['name'] = "linear_poisson_square_" + str(i)
            model_iga_include.PostMultiPatch(mpatch_mp_list[i].GetMultiPatch(), dim, time, params_post)

    return model0

def test():
    model = main2(logging=False, output=False)

    # print(len(model.model_part.Nodes))

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.5) < tol:
            temp = node.GetSolutionStepValue(TEMPERATURE)
            # print(temp)
            assert(abs(temp - 5.0) < 1e-12)

    print("Test passed")

def tag():
    return "thermal,gmg,iga"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main2(logging=True, output=True)
