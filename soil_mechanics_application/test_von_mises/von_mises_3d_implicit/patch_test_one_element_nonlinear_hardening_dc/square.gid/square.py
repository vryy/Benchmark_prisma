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

def main(output=True, logging=True):
    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
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
    if logging:
        ifile.write("0.0\t0.0\n")
        ifile.flush()

    disp = 0.0
    delta_disp = 1e-6
    delta_time = delta_disp
    for i in range(0, 200):
        disp = disp + delta_disp
        for node in prescribed_nodes:
            #node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)

        time = time + delta_time
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        if logging:
            ifile.write(str(disp) + "\t" + str(reaction_node.GetSolutionStepValue(REACTION_X)) + "\n")
            ifile.flush()

    if logging:
        ifile.close()

    return model

def test():
    model = main(output=False, logging=False)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol and abs(node.Y0 - 1.0) < tol:
            test_node = node

    ######pytesting######
    dy = test_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("dy: %.10e" % (dy))
    ref_disp = -1.9309803974e-04
    test = abs(dy - ref_disp) / abs(ref_disp)
    print(test)
    assert(test < 1e-10)
    #####################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=False, logging=True)
