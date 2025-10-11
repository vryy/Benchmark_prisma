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
import time
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
sys.path.append('./mesh_40x80.gid')
import mesh_40x80_include
from mesh_40x80_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################
# =====================
# | USER SCRIPT FOR CALCULATION OF EKATE.GID |
# vvvvvvvvvvvvvvvvvvvvv

def write_log(ifile, du, prescribed_nodes):
    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ifile.write("%-*.10e%.10e\n" % (20, du, reac_force_y))

def main(output=True, logging=True, output_last=False, total_disp=2.0, delta_disp=0.01):
    start = time.time()

    model = mesh_40x80_include.Model('mesh_40x80',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

        if abs(node.Y0 - 2.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("%-*s%s\n" % (20, "uy", "reaction-y"))

    mytime = 0.0
    model.Solve(mytime, 0, 0, 0, 0)

    disp = 0.0
    nsteps = int(total_disp / delta_disp)
    for i in range(0, nsteps):
        disp = disp + delta_disp
        print("disp:", disp)
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

        mytime = mytime + delta_disp
        model.Solve(mytime, 0, 0, 0, 0)
        if output:
            model.WriteOutput(mytime)
        if logging:
            write_log(ifile, disp, prescribed_nodes)

    if output_last and not output:
        model.WriteOutput(mytime)

    end = time.time()

    print("Total time: " + str(end-start) + " s")

    return model

def test():
    model = main(logging=False, output=False, output_last=False, total_disp=0.1, delta_disp=0.01)

    ### pytesting results
    tol = 1e-06
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 2.0) < tol:
            prescribed_nodes.append(node)

    reac_force_y = 0.0
    for node in prescribed_nodes:
        reac_force_y += node.GetSolutionStepValue(REACTION_Y)
    ref_reac = 1.081421966106692e-01
    print("%.15e" % reac_force_y)
    assert(abs(reac_force_y - ref_reac) / ref_reac < 1e-10)
    print("Test passed")

def tag():
    return "Neo-Hookean,Total-Lagrangian"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, output_last=True, total_disp=2.0, delta_disp=0.01)
