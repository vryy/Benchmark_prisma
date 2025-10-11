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
##################################################################
##################################################################
sys.path.append('./patch_test_2cube.gid')
import patch_test_2cube_include
from patch_test_2cube_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True, penalty=1e8):
    model = patch_test_2cube_include.Model('patch_test_2cube',os.getcwd()+"/", logging=logging)
    model.InitializeModel()

    # boundary condition
    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol or abs(node.Y0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        if abs(node.Z0 - 2.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
            prescribed_nodes.append(node)

    model.Solve(0.0, 0, 0, 0, 0)

    # create the correct mortar master/slave surfaces based on node_groups
    master_index_set = IntegerVector(1)
    master_index_set[0] = 10
    for cond in model.model_part.Conditions:
        is_master = True
        for node in cond.GetNodes():
            if node.Id not in model.node_groups['master']:
                is_master = False
                break
        if is_master == True:
            cond.SetValue(MASTER_INDEX_SET, master_index_set)

    slave_index_set = IntegerVector(1)
    slave_index_set[0] = 10
    for cond in model.model_part.Conditions:
        is_slave = True
        for node in cond.GetNodes():
            if node.Id not in model.node_groups['slave']:
                is_slave = False
                break
        if is_slave == True:
            cond.SetValue(SLAVE_INDEX_SET, slave_index_set)

    # setup mortar tying links
    util = MortarTyingUtility3D()
    util.SetDefaultSearchParameters()
    util.InitializeSearchTree(model.model_part, 10)
    util.SetEchoLevel(1)
    mortar_links = util.SetupTyingLinkElementBased(model.model_part, 10, "tying_link_kinematic_linear_penalty")
    for cond in mortar_links:
       cond.SetValue(INITIAL_PENALTY, penalty)

    # load increment
    disp = 0.0
    delta_disp = 0.1
    for i in range(0, 1):
        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Z, disp)
        model.Solve(disp, 0, 0, 0, 0)
        if output:
            model.WriteOutput(disp)

    return model

def test():
    model = main(logging=False, output=False, penalty=1e10)

    node1 = model.model_part.Nodes[11]
    node2 = model.model_part.Nodes[12]

    ref_disp1 = 4.9992933691102671e-02
    ref_disp2 = 5.0007066308896141e-02
    print("%.16e" % (node1.GetSolutionStepValue(DISPLACEMENT_Z)))
    print("%.16e" % (node2.GetSolutionStepValue(DISPLACEMENT_Z)))
    assert(abs(node1.GetSolutionStepValue(DISPLACEMENT_Z) - ref_disp1) < 1e-10)
    assert(abs(node2.GetSolutionStepValue(DISPLACEMENT_Z) - ref_disp2) < 1e-10)

    print("Test passed")

def tag():
    return "mesh-tying"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
