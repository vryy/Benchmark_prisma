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
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) or (abs(node.X0 - 1.0) < tol):
            node.Fix(DISPLACEMENT_Z)
        if (abs(node.Y0) < tol) or (abs(node.Y0 - 1.0) < tol):
            node.Fix(DISPLACEMENT_Z)
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)

    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(FACE_LOAD_Z, -1e3)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    ######### pytesting results #########
    mon_node = model.model_part.Nodes[15]
    ref_disp_z = -8.7969924812030128e-04
    disp_z = mon_node.GetSolutionStepValue(DISPLACEMENT_Z)
    print("%.16e" % (disp_z))
    assert(abs(disp_z - ref_disp_z) < 1e-12)
    print("Test passed")
    #####################################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
