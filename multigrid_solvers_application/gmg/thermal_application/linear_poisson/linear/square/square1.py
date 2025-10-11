##################################################################
##### Solution of linear Poisson problem using geometric     #####
#####                              multigrid method          #####
##### copyright Hoang-Giang Bui, 2018                        #####
#####          Ruhr University Bochum                        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
import math
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
##################################################################
##################################################################
import poisson_include
from poisson_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

mesh_util = FiniteCellMeshUtility()

##################################################################
##############PARAMETERS##########################################
##################################################################

nlevels = 3
coarse_m = 1
coarse_n = 1
fine_m = coarse_m*(2**(nlevels-1))
fine_n = coarse_n*(2**(nlevels-1))

sample_element_name = "LinearPoisson2D4N"
element_type = 1
results_path = os.getcwd()+"/"

start_point = Vector(3)
start_point[0] = 0.0
start_point[1] = 0.0
start_point[2] = 0.0

end_point = Vector(3)
end_point[0] = 1.0
end_point[1] = 1.0
end_point[2] = 0.0

##################################################################
####CONSTRUCT MULTIGRID HIERARCHY#################################
##################################################################

## generate the model_part for each level
model_list = []
for i in range(0, nlevels):
    name = "level_" + str(i)
    size_m = fine_m/(2**i)
    size_n = fine_n/(2**i)
    model_part = ModelPart(name)
    poisson_include.AddVariables(model_part)
    mesh_util.GenerateStructuredModelPart(model_part, start_point, end_point, size_m, size_n, sample_element_name, element_type, model_part.Properties[1])
    model = poisson_include.Model(model_part, results_path)
    print("Generate model_part for level " + str(i) + " completed")
    model.InitializeModel()
#    model.WriteOutput(0.0)
    model_list.append(model)
print(model_list[0].solver.solver.InitializeWasPerformed)
## boundary conditions
tol = 1.0e-6
temp = 10.0
for model in model_list:
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 0)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, temp)
    #    if abs(node.Y0 - 0.0) < tol:
    #        node.Fix(TEMPERATURE)
    #        node.SetSolutionStepValue(TEMPERATURE, 0)
    #    if abs(node.Y0 - 1.0) < tol:
    #        node.Fix(TEMPERATURE)
    #        node.SetSolutionStepValue(TEMPERATURE, temp)
    model.model_part.CloneTimeStep(0.0)

##################################################################
####SETUP GMG SOLVER##############################################
##################################################################

## construct multigrid solver for first level
#defining linear solver
plinear_solver = MultilevelSolver()
solver_plist = MGParameterList()
solver_plist.SetIntValue("num_levels", nlevels)
solver_plist.SetIntValue("block_size", 1)
solver_plist.SetIntValue("coarse_div_1", coarse_m)
solver_plist.SetIntValue("coarse_div_2", coarse_n)
plinear_solver.AddPreSmoother(JacobiIterativeSolver(3, 0.66))
plinear_solver.AddPostSmoother(JacobiIterativeSolver(3, 0.66))
plinear_solver.SetCoarseSolver(SkylineLUFactorizationSolver())

solver_factory = GMGStructuredSolverFactory2D(solver_plist)
for i in range(0, len(model_list)):
    solver_factory.SetModelPart(i, model_list[i].model_part)
solver_factory.InitializeMultilevelSolver(plinear_solver)
plinear_solver.SetFactory(solver_factory)

model_list[0].solver.structure_linear_solver = plinear_solver
model_list[0].solver.Initialize()
model_list[0].solver.solver.SetEchoLevel(2)
model_list[0].solver.solver.max_iter = 10 #control the maximum iterations of Newton Raphson loop
model_list[0].solver.solver.MoveMeshFlag = False
model_list[0].solver.solver.convergence_criteria = DisplacementCriteria(model_list[0].rel_tol, model_list[0].abs_tol)
model_list[0].solver.solver.builder_and_solver = ResidualBasedBlockBuilderAndSolver(plinear_solver)

##################################################################
###################SOLVE##########################################
##################################################################

level_list = []
for i in range(0, len(model_list)):
    level_list.append(plinear_solver.GetLevel(i))

params = {}
params['print_sparsity_info_flag'] = False
params['stop_Newton_Raphson_if_not_converge'] = True

import gmg_newton_raphson_strategy
gmg_solve_strategy = gmg_newton_raphson_strategy.SolvingStrategyPython(model_list, level_list, params)

time = 1.0

gmg_solve_strategy.Solve(time, 0, 0, 0, 0)

for model in model_list:
    model.WriteOutput(time)

