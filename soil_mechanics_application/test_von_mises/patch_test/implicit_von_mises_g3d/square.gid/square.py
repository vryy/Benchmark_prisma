##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
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
sys.path.append('./square.gid')
import square_include
from square_include import *

##################################################################
###  SIMULATION  #################################################
##################################################################
# =====================
# | USER SCRIPT FOR CALCULATION OF EKATE.GID |
# vvvvvvvvvvvvvvvvvvvvv

def main(logging=True, output=True, nsteps=150):
    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)
        if abs(node.X0 - 0.0) < tol and abs(node.Y0 - 0.0) < tol:
            reaction_node = node

    if logging:
        ifile = open("reaction.txt", "w")
        ifile.write("u\tf\n")

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)

    disp = 0.0
    delta_disp = 1e-6
    delta_time = delta_disp
    for i in range(0, nsteps):
        disp = disp + delta_disp
        for node in prescribed_nodes:
            #node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)

        print("************Step %f started************" % (disp))

        time = time + delta_time
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)

        if logging:
            ifile.write(str(disp) + "\t" + str(reaction_node.GetSolutionStepValue(REACTION_X)) + "\n")

    return model

def test():
    model = main(logging=False, output=False, nsteps = 150)

    test_node = model.model_part.Nodes[3]

    ######### pytesting results #########
    ref_disp_x = 1.5e-4
    ref_disp_y = -1.4414886466e-04
    disp_x = test_node.GetSolutionStepValue(DISPLACEMENT_X)
    disp_y = test_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("disp_x: %.10e" % disp_x)
    print("disp_y: %.10e" % disp_y)
    assert(abs(disp_x - ref_disp_x) / abs(ref_disp_x) < 1e-10)
    assert(abs(disp_y - ref_disp_y) / abs(ref_disp_y) < 1e-10)
    #####################################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, nsteps = 150)
