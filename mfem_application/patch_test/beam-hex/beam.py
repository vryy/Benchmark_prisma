##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Sa 14. Mar 00:15:32 CET 2020
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import beam_include
from beam_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(order, output=True):
    ###
    mesh_file = "beam-hex.mesh"
    vdim = 3
    fec = "H1_FECollection"

    utils = MfemUtilities()

    md = utils.CreateMfemModel(fec)
    md.Initialize(mesh_file, order, vdim)
    # print(md)

    utils.PrintFESpace(md)
    utils.PrintDofs(md, "Vertex")
    utils.PrintDofs(md, "Edge")
    utils.PrintDofs(md, "Face")
    utils.PrintDofs(md, "Interior")

    md.BeginModelPart()
    mp = md.GetModelPart()
    mp.SetBufferSize(2)
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( mp )
    md.CreateNodes()
    prop = mp.Properties[1]
    prop.SetValue(BODY_FORCE, ZeroVector(3) )
    prop.SetValue(CONSTITUTIVE_LAW, DummyConstitutiveLaw() )
    layer_index = 1
    if order == 1:
        prop.SetValue(INTEGRATION_ORDER, 2)
        md.AddElements(layer_index, "KinematicLinear3D8N", prop)
    elif order == 2:
        prop.SetValue(INTEGRATION_ORDER, 3)
        md.AddElements(layer_index, "KinematicLinear3D27N", prop)
    else:
        raise Exception("Order %d is not supported" % (order))
    md.EndModelPart()
    print(mp)
    ###

    model = beam_include.Model('beam',os.getcwd()+"/",os.getcwd()+"/",output)
    model.SetModelPart(mp)
    model.InitializeModel()

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            # print("node %d is fixed in x-direction" % (node.Id))
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
            # print("node %d is fixed in y-direction" % (node.Id))
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
            # print("node %d is fixed in z-direction" % (node.Id))
        if abs(node.X0 - 8.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    time = 1.0
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.3)

    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print()

    # for elem in model.model_part.Elements:
    #     print(elem.Id)
    #     for node in elem.GetNodes():
    #         print("%f,%f,%f" % (node.X0, node.Y0, node.Z0))

    ######pytesting######
    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))
    for node in model.model_part.Nodes:
        if (abs(node.Y0 - 1.0) < tol):
            dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
            assert(abs(dy + 0.01125) < 1e-10)
        if (abs(node.Z0 - 1.0) < tol):
            dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
            assert(abs(dz + 0.01125) < 1e-10)
    #####################

def test():
    main(1, False)
    main(2, False)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(2)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
