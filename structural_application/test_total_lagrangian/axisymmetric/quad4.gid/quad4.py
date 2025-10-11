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
sys.path.append('./quad4.gid')
import quad4_include
from quad4_include import *
# calculate insitu-stress for geology_virgin.gid

##################################################################
###  SIMULATION  #################################################
##################################################################

def main(logging=True, output=True):
    model = quad4_include.Model('quad4',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)

    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            node.SetSolutionStepValue(FACE_LOAD_X, 1.0e3)
            node.SetSolutionStepValue(FACE_LOAD_Y, 0.0)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    ######### pytesting results #########
    ref_disp_x = 0.000530094901587
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            disp_x = node.GetSolutionStepValue(DISPLACEMENT_X)
            assert(abs(disp_x - ref_disp_x) < 1e-12)
    #####################################

def test():
    main(output=False, logging=False)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output = False)
