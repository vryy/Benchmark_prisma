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
sys.path.append('./slope_stability_hyplas.gid')
import slope_stability_hyplas_include
from slope_stability_hyplas_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################
# =====================
# | USER SCRIPT FOR CALCULATION OF AUSBLAS.GID |
# vvvvvvvvvvvvvvvvvvvvv

def main(output=True, logging=True, delta_load_lists=[]):
    model = slope_stability_hyplas_include.Model('slope_stability_hyplas',os.getcwd()+"/", logging)
    model.InitializeModel()

    ##boundary condition
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 75.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.X0 - 35.0) < tol and abs(node.Y0 - 40.0) < tol:
            pointA = node

    ##################################################################

    if logging:
        ifile = open("settlement_a.txt", "w")
        ifile.write("dy\tload_fact\n")
        ifile.write("0.0\t0.0\n")

    load = 0.0
    for delta_load in delta_load_lists:
        load = load + delta_load
        print("Step load %f starts" % (load))
        model.model_part.Properties[1].SetValue(DENSITY, 2.038735 * load)
        time = load
        model.Solve( time, 0, 0, 0, 0 )

        if output:
            model.WriteOutput( time )
        if logging:
            ifile.write(str(pointA.GetSolutionStepValue(DISPLACEMENT_Y)) + "\t" + str(load) + "\n")

    if logging:
        ifile.close()

    print("ANALYSIS COMPLETED")
    return model

def test():
    delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.0625, 0.0078125] # load case for testing

    model = main(output=False, logging=False, delta_load_lists=delta_load_lists)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 35.0) < tol and abs(node.Y0 - 40.0) < tol:
            pointA = node

    ######pytesting######
    dy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    print(dy)
    ref_disp = -1.67885320893
    test = abs(dy - ref_disp) / abs(ref_disp)
    print(test)
    assert(test < 1e-9)
    #####################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        #delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.0625, 0.0625, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125]
        #delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125/16]
        delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.0625, 0.0078125] # load case for testing
        # delta_load_lists = [1.0, 1.0, 1.0]
        #delta_load_lists = [1.0, 1.0]
        # delta_load_lists = [1.0]

        main(output=True, logging=True, delta_load_lists=delta_load_lists)
