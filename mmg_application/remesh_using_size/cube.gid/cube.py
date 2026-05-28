##################################################################
import sys
import os
import math
import time as time_module
##################################################################
import cube_include
from cube_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True, output_gmsh=False, analysis=True, h_in=0.03, h_out=0.3):
    model = cube_include.Model('cube',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    for cond in model.model_part.Conditions:
        cond.Properties = model.model_part.Properties[2]

    model.WriteOutput(0.0)

    for node in model.model_part.Nodes:
        if (node.X0 > 0.4) and (node.X0 < 0.6):
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, h_in)
        else:
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, h_out)

    mmg_model = Mmg3DModel()
    mmg_model.SetValue(MMG_GRADATION, 1.1)
    mmg_model.SetValue(MMG_HAUSDORFF_DISTANCE, 0.001)
    mmg_model.SetIParam(MMG3D_Param.IPARAM_verbose, 10)
    mmg_model.Initialize(model.model_part)
    mmg_model.Remesh()

    mmg_model.BeginModelPart("Cube")
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
    prop3.Id = 3
    prop3.SetValue(LAYER_NAME, "Load")
    mp.AddProperties(prop3)

    mmg_model.CreateNodes()
    mmg_model.AddElements(1, "KinematicLinear3D4N", prop1)
    mmg_model.AddConditions(2, "FaceForce3D3N", prop3)
    mmg_model.EndModelPart()

    if output_gmsh:
        mmg_model.SaveMesh("cube.mesh")
        mmg_model.SaveMetric("cube.sol")
    model.SetModelPart(mp)
    print(mp)
    if output:
        model.WriteOutput(1.0)

    if analysis:
        print("############## ANALYSIS ##############")

        for prop in model.model_part.Properties:
            prop.SetValue(THICKNESS, 1.0)
            prop.SetValue(YOUNG_MODULUS, 2e5*prop.Id)

        tol = 1.0e-6
        for node in model.model_part.Nodes:
            if abs(node.X0) < tol:
                node.Fix(DISPLACEMENT_X)
            if abs(node.Y0) < tol:
                node.Fix(DISPLACEMENT_Y)
            if abs(node.Z0) < tol:
                node.Fix(DISPLACEMENT_Z)

        for node in model.model_part.Nodes:
            if (abs(node.Z0 - 1.0) < tol):
                node.SetSolutionStepValue(FACE_LOAD_X, 0.0)
                node.SetSolutionStepValue(FACE_LOAD_Y, 0.0)
                node.SetSolutionStepValue(FACE_LOAD_Z, 1e3)

        time = 100.0
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)

    return model

def test():

    model = main(output=False, logging=False, output_gmsh=False, analysis=True, h_in=0.05, h_out=0.5)

    print("num nodes: %d" % (len(model.model_part.Nodes)))
    print("num elements: %d" % (len(model.model_part.Elements)))

    if KratosMmgApplication.GetMmgVersion() == "5.5.2":

        assert(len(model.model_part.Nodes) == 7178)
        assert(len(model.model_part.Elements) == 38777)

        mon_node = model.model_part.Nodes[151]
        ux = mon_node.GetSolutionStepValue(DISPLACEMENT_X)
        uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
        uz = mon_node.GetSolutionStepValue(DISPLACEMENT_Z)
        print("ux: %.16e" % (ux))
        print("uy: %.16e" % (uy))
        print("uz: %.16e" % (uz))
        ux_ref = -6.5226439470000288e-04
        uy_ref = -1.9167209820000205e-04
        uz_ref = 1.2033793050000032e-03
        assert(abs(ux - ux_ref) < 1e-10)
        assert(abs(uy - uy_ref) < 1e-10)
        assert(abs(uz - uz_ref) < 1e-10)

    elif KratosMmgApplication.GetMmgVersion() == "5.8.0":

        assert(len(model.model_part.Nodes) == 7144)
        assert(len(model.model_part.Elements) == 38569)

        mon_node = model.model_part.Nodes[151]
        x = mon_node.X0
        y = mon_node.Y0
        z = mon_node.Z0
        ux = mon_node.GetSolutionStepValue(DISPLACEMENT_X)
        uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
        uz = mon_node.GetSolutionStepValue(DISPLACEMENT_Z)
        print("x: %.16e" % (x))
        print("y: %.16e" % (y))
        print("z: %.16e" % (z))
        print("ux: %.16e" % (ux))
        print("uy: %.16e" % (uy))
        print("uz: %.16e" % (uz))
        x_ref = 4.3632775847908845e-01
        y_ref = 1.2935882555313419e-01
        z_ref = 2.3621297167973243e-01
        ux_ref = -6.5449163771863147e-04
        uy_ref = -1.9403823832970105e-04
        uz_ref = 1.1810648583986602e-03
        assert(abs(x - x_ref) < 1e-10)
        assert(abs(y - y_ref) < 1e-10)
        assert(abs(z - z_ref) < 1e-10)
        assert(abs(ux - ux_ref) < 1e-10)
        assert(abs(uy - uy_ref) < 1e-10)
        assert(abs(uz - uz_ref) < 1e-10)

    else:
        raise Exception(f"This test is not implemented for MMg version {KratosMmgApplication.GetMmgVersion()}")

    print("Test passed")

def tag():
    return "mmg"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True, output_gmsh=False, analysis=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
