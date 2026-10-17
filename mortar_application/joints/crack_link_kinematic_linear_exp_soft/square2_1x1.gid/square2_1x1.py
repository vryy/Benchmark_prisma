##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025 by Hoang-Giang Bui (UoB)                 #####
#####          2026 by Hoang-Giang Bui (DU)                  #####
##### all rights reserved                                    #####
##################################################################
## This file is generated on (date)
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import square2_1x1_include as simulation_include
try:
    from square2_1x1_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'square2_1x1'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def GetReaction(nodes):
    reac = 0.0
    for node in nodes:
        reac += node.GetSolutionStepValue(REACTION_X)
    return reac

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    ## setup mortar tying links
    util = MortarTyingUtility2D()
    util.SetDefaultSearchParameters()
    util.InitializeSearchTree(model.model_part, 10)
    util.SetValue(TYING_INTEGRATION_ORDER, 3)
    util.SetEchoLevel(1)
    mortar_links = util.SetupTyingLinkElementBased(model.model_part, 10, "crack_link_kinematic_linear_exp_soft", model.model_part.Properties[2])
    joint_stiff = Matrix(2, 2)
    joint_stiff[0,0] = 222
    joint_stiff[1,1] = 99
    mat_params = Vector(2)
    mat_params[0] = 0.3
    mat_params[1] = 0.012
    model.model_part.Properties[2].SetValue(JOINT_STIFFNESS, joint_stiff)
    model.model_part.Properties[2].SetValue(MATERIAL_PARAMETERS, mat_params)

    ## boundary conditions
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 2.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("reaction.log", "w")
        ifile.write("disp\treaction\n")

    time = 0.0
    disp = 0.0
    model.SolveModel(time)
    if logging:
        reac = GetReaction(prescribed_nodes)
        ifile.write("%.6e\t%.10e\n" % (disp, reac))

    ## loading
    delta_disp = 1e-3/2
    delta_time = 1.0
    nsteps = 150*2
    for step in range(0, nsteps):
        time += delta_time
        disp += delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, disp)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            reac = GetReaction(prescribed_nodes)
            ifile.write("%.6e\t%.10e\n" % (disp, reac))

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 2.0) < tol:
            prescribed_nodes.append(node)

    reac = GetReaction(prescribed_nodes)
    # print("reac: %.16e" % (reac))
    ref_reac = 7.2973064910384321e-03
    assert(abs(reac - ref_reac) < 1e-10)
    print("Test passed")

def tag():
    tags = "joint"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
