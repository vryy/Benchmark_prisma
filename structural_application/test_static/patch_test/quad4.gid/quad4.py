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
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import quad4_include
from quad4_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = quad4_include.Model('quad4',current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)

    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.1)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.Id)
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    return model

def test():
    model = main(logging=False, output=False)

    ######### pytesting results #########
    tol = 1.0e-6
    ref_disp_y = -8.5714285714285701e-03
    mon_node = model.model_part.Nodes[35]
    disp_y = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
    # print("%.16e" % (disp_y))
    assert(abs(disp_y - ref_disp_y) < 1e-12)
    print("Test passed")
    #####################################

def tag():
    return "kinematic_linear,linear_elastic"

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)
