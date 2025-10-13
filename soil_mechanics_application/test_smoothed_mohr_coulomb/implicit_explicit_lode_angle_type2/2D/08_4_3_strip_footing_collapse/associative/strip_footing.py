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
sys.path.append('./strip_footing.gid')
import strip_footing_include
from strip_footing_include import *

##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = strip_footing_include.Model('strip_footing',os.getcwd()+"/",logging=logging,material_type="smc-2")
    model.InitializeModel()

    ####boundary condition
    prescribed_nodes = []
    reaction_nodes = []
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol or abs(node.X0 - 500.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            reaction_nodes.append(node)
        if abs(node.Y0 - 500.0) < tol and node.X0 > -tol and node.X0 < 50.0+tol:
            prescribed_nodes.append(node)
            node.Fix(DISPLACEMENT_Y)
        node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, 0.0)
        node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z, 0.0)
    ##################################################################
    if logging:
        ## prepare reporting
        ifile = open("report_smc.grf", "w")
        c = model.model_part.Properties[1].GetValue(COHESION)

        B = 100.0
        ifile.write("Normalized-pressure\t\tNormalized-settlement\n")
    ##################################################################
    disp = 0.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
    model.Solve(disp, 0, 0, 0, 0)
    if output:
        model.WriteOutput(disp)
    if logging:
        P = 0.0
        ifile.write(str(2.0*P/B/c) + "\t" + str(disp/B) + "\n")
    ##################################################################
    disp_max = 100;
    delta_disp_list = []
    delta_disp_list = [0.0005, 0.0005, 0.00025, 0.00025, 0.0005, 0.0005]

    # elem_id = 133
    # point_id = 3
    # plot_util = StressPathPlotUtility("stress_path_elem_"+str(elem_id)+"_"+str(point_id)+".txt")

    for delta_disp in  delta_disp_list:
        disp = disp + delta_disp * disp_max
        for node in prescribed_nodes:
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, -delta_disp * disp_max)
        model.Solve(disp, 0, 0, 0, 0)
        # plot_util.RegisterPoint(disp, model.model_part, elem_id, point_id)
        if output:
            model.WriteOutput(disp)
        if logging:
            P = 0.0
            for node in reaction_nodes:
                P = P + node.GetSolutionStepValue(REACTION_Y)
            ifile.write(str(2.0*P/B/c) + "\t" + str(disp/B) + "\n")
    ##################################################################
    if logging:
        ifile.close()

    return model

def test():
    model = main(output=False, logging=False)

    reaction_nodes = []
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.Y0) < tol:
            reaction_nodes.append(node)

    ######pytesting######
    P = 0.0
    c = model.model_part.Properties[1].GetValue(COHESION)
    B = 100.0
    for node in reaction_nodes:
        P = P + node.GetSolutionStepValue(REACTION_Y)
    print("P/B/c: %.16e" % (2.0*P/B/c))
    ref_value = 1.4975331791180105e+01
    assert(abs(2.0*P/B/c - ref_value) / abs(ref_value) < 1e-10)
    #####################
    print("Test passed")

def tag():
    return "S-MC"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
