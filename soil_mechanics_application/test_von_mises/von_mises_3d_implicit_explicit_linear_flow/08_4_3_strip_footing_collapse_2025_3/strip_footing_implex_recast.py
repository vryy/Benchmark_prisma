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

def main(output=True, logging=True, material_type="j2-implex-recast", ext_option=1, load_option=2):
    model = strip_footing_include.Model('strip_footing',os.getcwd()+"/",logging=logging)
    model.material_type = material_type
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
        if model.material_type == "j2":
            ifile = open("report_j2-" + str(load_option) + ".grf", "w")
            c = model.model_part.Properties[1].GetValue(TENSILE_STRENGTH) / math.sqrt(3)
        elif model.material_type == "j2-implex":
            ifile = open("report_j2_implex-" + str(ext_option) + "-" + str(load_option) + ".grf", "w")
            c = model.model_part.Properties[1].GetValue(TENSILE_STRENGTH) / math.sqrt(3)
        elif model.material_type == "j2-implex-recast":
            ifile = open("report_j2_implex_recast-" + str(ext_option) + "-" + str(load_option) + ".grf", "w")
            c = model.model_part.Properties[1].GetValue(TENSILE_STRENGTH) / math.sqrt(3)
        elif model.material_type == "j2-implex-linearized-flow":
            ifile = open("report_j2_implex_linearized_flow-" + str(ext_option) + "-" + str(load_option) + ".grf", "w")
            c = model.model_part.Properties[1].GetValue(TENSILE_STRENGTH) / math.sqrt(3)
        elif model.material_type == "j2-implex-elastic":
            ifile = open("report_j2_implex_elastic-" + str(ext_option) + "-" + str(load_option) + ".grf", "w")
            c = model.model_part.Properties[1].GetValue(TENSILE_STRENGTH) / math.sqrt(3)
        elif model.material_type == "tr":
            ifile = open("report_tr-" + str(load_option) + ".grf", "w")
            c = 0.5*model.model_part.Properties[1].GetValue(TENSILE_STRENGTH)
        else:
            raise Exception(f"Invalid material type {material_type}")

    if logging:
        B = 100.0
        ifile.write("Normalized-pressure\t\tNormalized-settlement\n")
    ##################################################################
    disp = 0.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
    model.model_part.ProcessInfo[EXTRAPOLATION_FACTOR_1] = 1.0
    model.model_part.ProcessInfo[EXTRAPOLATION_FACTOR_2] = 0.0
    model.SolveModel(disp)
    if output:
        model.WriteOutput(disp)
    if logging:
        P = 0.0
        ifile.write(str(2.0*P/B/c) + "\t" + str(disp/B) + "\n")
    ##################################################################
    ## for implex
    ##################################################################
    if ext_option == 1:
        ## extrapolation (first order accuracy)
        factor1 = 1.0
        factor2 = 0.0
    elif ext_option == 2:
        # ### extrapolation (second order accuracy)
        factor1 = 2.0
        factor2 = -1.0
    #
    if load_option == 1:
        ### coarse load (irregular)
        delta_disp_list = [0.0001, 0.00005, 0.00005, 0.00005, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.00005, 0.0001, 0.0002, 0.0003, 0.0006]
        # delta_disp_list = [0.0001, 0.00005]
        # delta_disp_list = [0.00015]
    elif load_option == 2:
        ### medium fine load (regular)
        delta_disp_list = []
        for i in range(0, 100):
            delta_disp_list.append(0.00002)
    elif load_option == 3:
        ### fine load (regular)
        delta_disp_list = []
        for i in range(0, 2000):
            delta_disp_list.append(0.000001)
    ######################
    disp_max = 100.0
    old_old_delta_disp = delta_disp_list[0]
    old_delta_disp = delta_disp_list[0]
    for delta_disp in delta_disp_list:
        model.model_part.ProcessInfo[EXTRAPOLATION_FACTOR_1] = factor1*delta_disp/old_delta_disp
        model.model_part.ProcessInfo[EXTRAPOLATION_FACTOR_2] = factor2*delta_disp/old_old_delta_disp
        disp = disp + delta_disp * disp_max
        for node in prescribed_nodes:
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, -delta_disp * disp_max)
        model.SolveModel(disp)
        if output:
            model.WriteOutput(disp)
        if logging:
            P = 0.0
            for node in reaction_nodes:
                P = P + node.GetSolutionStepValue(REACTION_Y)
            ifile.write(str(2.0*P/B/c) + "\t" + str(disp/B) + "\n")
        old_old_delta_disp = old_delta_disp
        old_delta_disp = delta_disp
    ##################################################################
    if logging:
        ifile.close()

    return model

def tag():
    return "j2,implex"

def print_tag():
    print("Tag(s): " + tag())

def test():
    model = main(logging=False, output=False, material_type="j2-implex-recast", ext_option=1, load_option=2)

    reaction_nodes = []
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.Y0) < tol:
            reaction_nodes.append(node)

    reac = 0.0
    for node in reaction_nodes:
        reac = reac + node.GetSolutionStepValue(REACTION_Y)

    B = 100.0
    c = model.model_part.Properties[1].GetValue(TENSILE_STRENGTH) / math.sqrt(3)

    reac_scaled = 2.0*reac/B/c
    print("%.16e" % (reac_scaled))
    ref_reac_scaled = 5.1893717012135374e+00
    assert(abs(reac_scaled - ref_reac_scaled) < 1e-10)

    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        # main(logging=True, output=False, material_type="j2", ext_option=1, load_option=1)
        main(logging=True, output=False, material_type="j2-implex-recast", ext_option=2, load_option=3)
