##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import square_include
from square_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True, output_gmsh=False):
    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    if output:
        model.WriteOutput(0.0)

    ### refinement indicator

    layer1 = LinearLevelSet(0.0, 1.0, -0.2)
    layer2 = LinearLevelSet(0.0, 1.0, -0.4)
    layer3 = LinearLevelSet(0.0, 1.0, -0.6)

    midp = [0.5, 0.6]
    a1 = -0.1
    b1 = 1.0
    a2 = -0.2
    b2 = 1.0
    layer4 = UnionLevelSet(LinearLevelSet(a1, b1, -(a1*midp[0]+b1*midp[1])), LinearLevelSet(a2, b2, -(a2*midp[0]+b2*midp[1])))

    time = 0.0
    delta_time = 1.0

    ## remeshing on soil layer

    it = 0
    for layer in [layer1, layer2, layer3, layer4]:

        print("###############################")
        print("## remeshing for layer " + str(it+1) + " ##")
        print("###############################")

        ##
        ## assign metric
        ##

        for node in model.model_part.Nodes:
            ls_val = layer.GetValue(node.X0, node.Y0, node.Z0)
            node.SetSolutionStepValue(NODAL_MMG_LEVEL_SET, ls_val)
            node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 0.1)

        if output:
            model.WriteOutput(time + 0.1)

        time = time + delta_time

        mmg_model = Mmg2DModel(True, True)
        mmg_model.SetValue(MMG_GRADATION, 1.3)
        mmg_model.SetValue(MMG_HAUSDORFF_DISTANCE, 1.0e-5)
        mmg_model.SetIParam(MMG2D_Param.IPARAM_verbose, 10)
        # mmg_model.SetIParam(MMG2D_Param.IPARAM_nreg, 1) # normal regularization
        # mmg_model.SetDParam(MMG2D_Param.DPARAM_rmc, 1.0e-1)
        # mmg_model.SetIParam(MMG2D_Param.IPARAM_angle, 1)
        # mmg_model.SetDParam(MMG2D_Param.DPARAM_angleDetection, 30.0)
        nmat = it + 1
        mmg_model.SetIParam(MMG2D_Param.IPARAM_numberOfMat, nmat)
        for i in range(0, nmat-1):
            mmg_model.SetLevelSetNoSplitting(i+1)
        mmg_model.SetLevelSetSplitting(nmat, nmat+1, nmat+2)
        mmg_model.Initialize(model.model_part)
        mmg_model.Remesh()
        if output_gmsh:
            mmg_model.SaveMesh("square_" + str(time+0.1) + ".mesh")
            mmg_model.SaveMetric("square_" + str(time+0.1) + ".sol")
        mmg_model.BeginModelPart("square")
        mp = mmg_model.GetModelPart()
        model.AddVariables(mp)

        mmg_model.CreateNodes()

        for i in range(0, nmat+1):
            prop = Properties(model.model_part.Properties[1])
            prop.Id = i+1
            prop.SetValue(LAYER_NAME, "layer"+str(i+1))
            mp.AddProperties(prop)

        for i in range(0, nmat-1):
            mmg_model.AddElements(i+1, "KinematicLinear2D3N", mp.Properties[i+1])
        mmg_model.AddElements(nmat+1, "KinematicLinear2D3N", mp.Properties[nmat])
        mmg_model.AddElements(nmat+2, "KinematicLinear2D3N", mp.Properties[nmat+1])
        mmg_model.AddConditions(1, "LineForce2D2N", mp.Properties[1]) # for the line force on top
        mmg_model.EndModelPart()

        model.SetModelPart(mp)
        # print(mp)
        if output:
            model.WriteOutput(time)

        it = it + 1

    nmat = nmat+1

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

    print(f"len(model.model_part.Nodes): {len(model.model_part.Nodes)}")
    print(f"len(model.model_part.Elements): {len(model.model_part.Elements)}")

    if KratosMmgApplication.GetMmgVersion() == "5.5.2":

        assert(len(model.model_part.Nodes) == 335)
        assert(len(model.model_part.Elements) == 628)

        mon_node = model.model_part.Nodes[151]
        ux = mon_node.GetSolutionStepValue(DISPLACEMENT_X)
        uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
        print("ux: %.16e" % (ux))
        print("uy: %.16e" % (uy))
        ux_ref = -3.0906117989366655e-04
        uy_ref = 1.6505530938269002e-03
        assert(abs(ux - ux_ref) < 1e-10)
        assert(abs(uy - uy_ref) < 1e-10)

    elif KratosMmgApplication.GetMmgVersion() == "5.8.0":

        assert(len(model.model_part.Nodes) == 331)
        assert(len(model.model_part.Elements) == 620)

        mon_node = model.model_part.Nodes[151]
        x = mon_node.X0
        y = mon_node.Y0
        ux = mon_node.GetSolutionStepValue(DISPLACEMENT_X)
        uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
        print("x: %.16e" % (x))
        print("y: %.16e" % (y))
        print("ux: %.16e" % (ux))
        print("uy: %.16e" % (uy))
        x_ref = 4.8936774329093136e-01
        y_ref = 0.6
        ux_ref = -3.0879357652462055e-04
        uy_ref = 1.6510271508482241e-03
        assert(abs(x - x_ref) < 1e-10)
        assert(abs(y - y_ref) < 1e-10)
        assert(abs(ux - ux_ref) < 1e-10)
        assert(abs(uy - uy_ref) < 1e-10)

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
        main(output=True, logging=True, output_gmsh=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
