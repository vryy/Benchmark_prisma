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

def main(logging=True, output=True, output_gmsh=False):
    # read in the model_part
    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    if output:
        model.WriteOutput(0.0)

    for node in model.model_part.Nodes:
        if (node.X0 > 0.45) and (node.X0 < 0.55):
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 0.005)
        else:
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 0.1)

    # remeshing
    mmg_model = Mmg2DModel()
    mmg_model.SetValue(MMG_GRADATION, 1.1)
    mmg_model.SetValue(MMG_HAUSDORFF_DISTANCE, 0.001)
    mmg_model.SetIParam(MMG2D_Param.IPARAM_verbose, 10)
    mmg_model.Initialize(model.model_part)
    mmg_model.Remesh()

    # regenerate new model_part
    mmg_model.BeginModelPart("Square")
    mp = mmg_model.GetModelPart()
    model.AddVariables(mp)

    prop1 = Properties(model.model_part.Properties[1])
    prop1.Id = 1
    prop1.SetValue(LAYER_NAME, "Bulk")
    mp.AddProperties(prop1)

    prop2 = Properties(model.model_part.Properties[2])
    prop2.Id = 2
    prop2.SetValue(LAYER_NAME, "Embeb")
    mp.AddProperties(prop2)

    prop3 = Properties(model.model_part.Properties[3])
    prop3.SetValue(THICKNESS,            1 )
    prop3.Id = 3
    prop3.SetValue(LAYER_NAME, "Load")
    mp.AddProperties(prop3)

    mmg_model.CreateNodes()
    mmg_model.AddElements(1, "KinematicLinear2D3N", prop1)
    mmg_model.AddElements(2, "KinematicLinear2D3N", prop2)
    mmg_model.AddConditions(1, "LineForce2D2N", prop3)
    mmg_model.EndModelPart()

    if output_gmsh:
        mmg_model.SaveMesh("square.mesh")
        mmg_model.SaveGmsh("square.msh")
        #mmg_model.SaveMetric("square.sol")
    model.SetModelPart(mp)
    print(mp)
    if output:
        model.WriteOutput(1.0)

    # analysis

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

    assert(len(model.model_part.Nodes) == 7861)
    assert(len(model.model_part.Elements) == 15436)

    mon_node = model.model_part.Nodes[1188]
    ux = mon_node.GetSolutionStepValue(DISPLACEMENT_X)
    uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print("ux: %.16e" % (ux))
    print("uy: %.16e" % (uy))
    ux_ref = -1.0097568362405652e-03
    uy_ref = 1.3540436472784662e-03
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
