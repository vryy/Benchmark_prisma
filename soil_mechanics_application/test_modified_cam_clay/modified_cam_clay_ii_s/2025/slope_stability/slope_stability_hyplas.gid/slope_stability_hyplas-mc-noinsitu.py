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

def main(output=True, logging=True, visual=False):

    model = slope_stability_hyplas_include.Model('slope_stability_hyplas',os.getcwd()+"/",logging=logging)
    model.material = "mohr-coulomb"
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

    # delta_load_lists = [1.0, 1.0, 1.0/2, 1.0/8, 1.0/16] # sum = 2.6875
    delta_load_lists = [1.0, 1.0, 1.0/2, 1.0/8, 1.0/16, 1.0/32] # sum = 2.71875
    # delta_load_lists = [1.0]

    if logging:
        ifile = open("settlement_a-mc-noinsitu.txt", "w")
        ifile.write("dy\tload_fact\n")

    load = 0.0
    first_disp = {}
    for delta_load in delta_load_lists:
        load = load + delta_load
        model.model_part.Properties[1].SetValue(DENSITY, 2.038735 * load)
        time = load
        model.Solve( time, 0, 0, 0, 0 )

        # record the displacement for first step
        if load == 1.0:
            print(f"Record displacements at load step {load}")
            for node in model.model_part.Nodes:
                first_disp[node.Id] = []
                first_disp[node.Id].append(node.GetSolutionStepValue(DISPLACEMENT_X))
                first_disp[node.Id].append(node.GetSolutionStepValue(DISPLACEMENT_Y))

        if output:
            model.WriteOutput( time )
        if logging:
            u = pointA.GetSolutionStepValue(DISPLACEMENT_Y) - first_disp[pointA.Id][1]
            ifile.write("%.16e\t%.16e\n" % (u, load))
            ifile.flush()

    if output:
        ifile.close()

    print("ANALYSIS COMPLETED")

    if visual:
        gui = PostProcessingGUI(model.model_part)
        var_list = VariablesList()
        var_list.Add(DISPLACEMENT)
        var_list2 = VariablesList()
        var_list2.Add(PLASTICITY_INDICATOR)
        gui.SetNodalSolutionStepVariablesList(var_list)
        gui.SetElementalSolutionStepVariablesList(var_list2)
        gui.Run()

    return model

def test():
    model = main(output=False, logging=False, visual=False)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 35.0) < tol and abs(node.Y0 - 40.0) < tol:
            pointA = node

    ######pytesting######
    dy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    print("dy: %.16e" % (dy))
    # print(dy)
    ref_disp = -1.60421242983
    assert(abs(dy - ref_disp) / abs(ref_disp) < 1e-10)
    #####################
    print("Test passed")

def tag():
    return "mohr-coulomb"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True, visual=False)
