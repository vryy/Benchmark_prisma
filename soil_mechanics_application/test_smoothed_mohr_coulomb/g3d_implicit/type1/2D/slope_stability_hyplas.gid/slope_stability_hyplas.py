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

def main(output=True, logging=True):

    model = slope_stability_hyplas_include.Model('slope_stability_hyplas',os.getcwd()+"/",logging)
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
            print("pointA.Id: " + str(pointA.Id))

    ##################################################################
    #delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.0625, 0.0625, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125]
    #delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125/16]
    # delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.0625, 0.0078125]

    delta_load_lists = [1.0, 1.0, 1.0, 1.0, 0.125, 0.03, 0.01, 0.005]
    delta_load_lists.extend([0.0025]*4)
    delta_load_lists.extend([0.00125]*7)
    delta_load_lists.extend([0.000625]*8)
    delta_load_lists.extend([0.0003125]*5)# 4.1953125

    if output:
        ifile = open("settlement_a.txt", "w")
        ifile.write("dy\tload_fact\n")
        ifile.write("0.0\t0.0\n")

    load = 0.0
    for delta_load in delta_load_lists:
        load = load + delta_load
        model.model_part.Properties[1].SetValue(DENSITY, 2.038735 * load)
        time = load
        model.Solve( time, 0, 0, 0, 0 )

        if output:
            model.WriteOutput( time )
        if logging:
            ifile.write(str(pointA.GetSolutionStepValue(DISPLACEMENT_Y)) + "\t" + str(load) + "\n")
            ifile.flush()

    if output:
        ifile.close()

    ######pytesting######
    dy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    print("%.15e" % dy)
    ref_disp = -1.481813082425050e+00
    assert(abs(dy - ref_disp) / abs(ref_disp) < 1e-8)
    #####################

    print("ANALYSIS COMPLETED")

def test():
    main(output=False, logging=False)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True)
