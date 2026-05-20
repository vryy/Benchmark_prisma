##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import square_include
from square_include import *

##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def SetMetric(model_part, level_set):
    for node in model_part.Nodes:
        ls_val = level_set.GetValue(node.X0, node.Y0, 0.0)
        # print("Node " + str(node.Id) + " level set: " + str(ls_val))
        node.SetSolutionStepValue(NODAL_MMG_LEVEL_SET, ls_val)

        if abs(ls_val) < 0.01:
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 0.01)
        else:
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 1.0)

def main(logging=True, output=True, output_gmsh=False):
    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    embed1 = CircularLevelSet(0.5, 0.7, 0.1)
    embed2 = CircularLevelSet(0.5, 0.3, 0.1)
    embed = UnionLevelSet(embed1, embed2)

    SetMetric(model.model_part, embed)

    time = 0.0
    if output:
        model.WriteOutput(time)

    brep_util = BRepUtility()

    ## REMESHING ITERATION

    delta_time = 1.0
    for i in range(0, 2):
        time = time + delta_time

        use_level_set = True
        use_metric = True
        mmg_model = Mmg2DModel(use_level_set, use_metric)
        # mesher.SetValue(MMG_GRADATION, 1.1)
        mmg_model.SetValue(MMG_HAUSDORFF_DISTANCE, 0.01)
        mmg_model.SetIParam(MMG2D_Param.IPARAM_verbose, 10)
        mmg_model.SetDParam(MMG2D_Param.DPARAM_ls, 0.0) # level set value to refine
        # mesher.SetDParam(MMG2D_Param.DPARAM_rmc, 0.01) # Remove small connex components in level-set mode
        mmg_model.Initialize(model.model_part)
        mmg_model.Remesh()

        mmg_model.BeginModelPart("refined_square")
        mp = mmg_model.GetModelPart()
        model.AddVariables(mp)

        prop1 = Properties(model.model_part.Properties[1])
        prop1.Id = 1
        mp.AddProperties(prop1)
        mp.Properties[1].SetValue(LAYER_NAME, "Bulk")

        prop2 = Properties(model.model_part.Properties[1])
        prop2.Id = 2
        mp.AddProperties(prop2)
        mp.Properties[2].SetValue(LAYER_NAME, "P2")

        prop3 = Properties(model.model_part.Properties[1])
        prop3.Id = 3
        mp.AddProperties(prop3)
        mp.Properties[3].SetValue(LAYER_NAME, "P3")

        prop4 = Properties(model.model_part.Properties[1])
        prop4.Id = 4
        mp.AddProperties(prop4)
        mp.Properties[4].SetValue(LAYER_NAME, "Edge")

        mmg_model.CreateNodes()
        # after the remeshing, the internal MMG mesh is spawned with 2 and 3 as marker of the new region
        # and 10 is the interface between domain
        mmg_model.AddElements(2, "KinematicLinear2D3N", prop1)
        mmg_model.AddElements(3, "KinematicLinear2D3N", prop2)
        mmg_model.AddConditions(1, "LineForce2D2N", prop1) # add force on top
        mmg_model.AddConditions(10, "LineForce2D2N", prop4)

        mmg_model.EndModelPart()

        # if i > -1:
        #     for elem in mp.Elements:
        #         center = brep_util.ComputeCenter(elem)

        #         ls_val1 = embed1.GetValue(center)
        #         if ls_val1 < 0.0:
        #             elem.Properties = prop2

        #         ls_val2 = embed2.GetValue(center)
        #         if ls_val2 < 0.0:
        #             elem.Properties = prop3

        #         if (ls_val1 > 0.0) and (ls_val2 > 0.0):
        #             elem.Properties = prop1

        if output_gmsh:
            mesher.SaveMesh("square.mesh")
            mesher.SaveMetric("square.sol")
        model.SetModelPart(mp)
        SetMetric(model.model_part, embed)
        # print(mp)
        if output:
            model.WriteOutput(time)

    # analysis

    for prop in model.model_part.Properties:
        prop.SetValue(THICKNESS, 1.0)
        prop.SetValue(YOUNG_MODULUS, 2e5*prop.Id)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)

    for node in model.model_part.Nodes:
        if (abs(node.Y0 - 1.0) < tol):
            node.SetSolutionStepValue(FACE_LOAD_X, 0.0)
            node.SetSolutionStepValue(FACE_LOAD_Y, 1e3)

    time = 100.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    return model

def test():

    model = main(output=False, logging=False, output_gmsh=False)

    assert(len(model.model_part.Nodes) == 1162)
    assert(len(model.model_part.Elements) == 2234)

    mon_node = model.model_part.Nodes[684]
    ux = mon_node.GetSolutionStepValue(DISPLACEMENT_X)
    uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("ux: %.16e" % (ux))
    print("uy: %.16e" % (uy))
    ux_ref = -1.0797686208854366e-03
    uy_ref = 1.5528509112974270e-03
    assert(abs(ux - ux_ref) < 1e-10)
    assert(abs(uy - uy_ref) < 1e-10)

    print("Test passed")

def tag():
    return "mmg"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
