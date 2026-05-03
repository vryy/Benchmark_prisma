##################################################################
# Pinching cylinder problem, analysis with quarter cylinder (single patch)
# Ref:
# [1] Nguyen et al, Rotation free isogeometric thin shell analysis using PHT-splines
# reference solution: 1.8248e-5 [1]
##################################################################

#including kratos path
import sys
import os
import numpy

from KratosMultiphysics.StructuralApplication import *

# loading the problem model
import QuarterCylinderModel
from QuarterCylinderModel import *

def main(logging=True, output=True,nsampling=40,order=2,thickness=3.0):

    solver = QuarterCylinderModel.SolutionScheme('PinchedCylinder',os.getcwd()+"/",nsampling=nsampling,order=order,logging=logging)

    # assign input data to variables
    solver.InitializeModel(thickness=thickness)

    R= 300.0
    L= 600.0
    loading_node1 = []

    fixed_nodes = []
    tol = 1.0e-4
    for node in solver.model_part.Nodes:

        if abs(node.Y0-0.0) < tol:
            node.Fix( DISPLACEMENT_X )
            node.Fix( DISPLACEMENT_Z )

        if abs(node.X0) < tol and abs(node.Z0-R) < tol:
            node.Fix( DISPLACEMENT_X )

        if abs(node.X0+R) < tol and abs(node.Z0) < tol:
            node.Fix( DISPLACEMENT_Z )

        if abs(node.Y0-L/2.0) < tol:
            node.Fix( DISPLACEMENT_Y )


        if abs(node.X0) < tol and abs(node.Z0-R) < tol and abs(node.Y0-L/2.0) < tol:
            loading_node1.append(node)

    # find the coordinates to restrain the rotation
    u2c = solver.FindInnerCoordinate("u1")
    v2c = solver.FindInnerCoordinate("v1")

    #############################
    penalty = 1.0e12
    joint_stiffness = Matrix(6, 6)
    for i in range(0, 6):
        for j in range(0, 6):
            joint_stiffness[i, j] = 0.0

    for i in range(0, 3):
        joint_stiffness[i, i] = penalty
        joint_stiffness[i, i+3] = -penalty
        joint_stiffness[i+3, i] = -penalty
        joint_stiffness[i+3, i+3] = penalty

    print(joint_stiffness)
    solver.model_part.Properties[1].SetValue(JOINT_STIFFNESS, joint_stiffness)

    RowU1= []
    RowU2= []

    RowV1=[]
    RowV2= []
    RowV3= []
    RowV4= []

    for node in solver.model_part.Nodes:
        # restrain the rotation at y=L/2 (front edge)
        if abs(node.Y0-L/2.0)<tol:
            RowU1.append(node.Id)
        if abs(node.Y0-u2c)<tol:
            RowU2.append(node.Id)

        # restrain the rotation at z=0 (bottom edge)
        if abs(node.X0+R)<tol and abs(node.Z0)<tol:
            RowV1.append(node.Id)
        if abs(node.X0+R)<tol and abs(node.Z0-v2c)<tol:
            RowV2.append(node.Id)

        # restrain the rotation at x=0 (left edge)
        if abs(node.X0-0.0)<tol and abs(node.Z0-R)<tol:
            RowV3.append(node.Id)
        if abs(node.X0+v2c)<tol and abs(node.Z0-R)<tol :
            RowV4.append(node.Id)

    #print(RowU1)
    #print(RowV3)
    #sys.exit()
    cnt_cnd = 1
    for i in range(0,len(RowU1)):
        node1= solver.model_part.Nodes[RowU1[i]]
        node2= solver.model_part.Nodes[RowU2[i]]
        cond1 = PointPointJointCondition(cnt_cnd, node1, node2, solver.model_part.Properties[1])
        cnt_cnd = cnt_cnd + 1
        solver.model_part.AddCondition(cond1)

    for i in range(0,len(RowV1)):
        node1= solver.model_part.Nodes[RowV1[i]]
        node2= solver.model_part.Nodes[RowV2[i]]
        cond3 = PointPointJointCondition(cnt_cnd, node1, node2, solver.model_part.Properties[1])
        cnt_cnd = cnt_cnd + 1
        solver.model_part.AddCondition(cond3)

        node3= solver.model_part.Nodes[RowV3[i]]
        node4= solver.model_part.Nodes[RowV4[i]]
        cond4 = PointPointJointCondition(cnt_cnd, node3, node4, solver.model_part.Properties[1])
        cnt_cnd = cnt_cnd + 1
        solver.model_part.AddCondition(cond4)

    #print load_node
    cnt_cnd = cnt_cnd + 1
    for load_node in loading_node1:
        load_node.SetSolutionStepValue(FORCE_X, 0.0)
        load_node.SetSolutionStepValue(FORCE_Y, 0.0)
        load_node.SetSolutionStepValue(FORCE_Z, -0.25)
        cond5 = PointForce3D(cnt_cnd, load_node, solver.model_part.Properties[1])

        solver.model_part.AddCondition(cond5)
        print("force is created at node " + str(load_node.Id))

    #####################################
    ########### solve the model #########
    time = 1.0
    solver.Solve(time)

    #####################################
    ######## post-processing ############
    if output:
        solver.WriteOutput(time )

    for node in loading_node1:
         print(node.GetSolutionStepValue(DISPLACEMENT_Z))

    return solver

def test():
    solver = main(logging=False, output=False, nsampling=40, order=2)

    R = 300.0
    L = 600.0
    loading_node1 = []
    tol = 1e-4
    for node in solver.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Z0-R) < tol and abs(node.Y0-L/2.0) < tol:
            loading_node1.append(node)

    ref_u = -1.775679963897207e-05
    for node in loading_node1:
        u = node.GetSolutionStepValue(DISPLACEMENT_Z)
        assert(abs(u - ref_u) < 1e-10)

    print("Test passed")

def tag():
    return "KL-shell"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True, nsampling=40, order=2, thickness=3.0)
