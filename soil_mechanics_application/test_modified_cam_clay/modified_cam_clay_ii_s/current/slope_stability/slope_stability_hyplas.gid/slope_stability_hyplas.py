##################################################################
import sys
import os
import math
##################################################################
##################################################################
import slope_stability_hyplas_include
from slope_stability_hyplas_include import *
##################################################################
###  SIMULATION  #################################################
### This test is to make sure MCC-S works with prestress
##################################################################

def main(output=True, last_output=False, logging=True, visual=False):

    ############## INSITU MODEL ##############

    model_insitu = slope_stability_hyplas_include.Model('slope_stability_hyplas',os.getcwd()+"/",logging=False)
    model_insitu.material = "linear-elastic"
    model_insitu.InitializeModel()

    ##boundary condition
    tol = 1.0e-6
    for node in model_insitu.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 75.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

    model_insitu.model_part.Properties[1].SetValue(DENSITY, 2.038735)
    time = 1.0
    model_insitu.Solve( time, 0, 0, 0, 0 )
    if output:
        model_insitu.WriteOutput(time)

    ############## SYSTEM MODEL ##############

    model = slope_stability_hyplas_include.Model('slope_stability_hyplas',os.getcwd()+"/",logging)
    model.material = "mcc-s"
    model.InitializeModel()

    ## set the overconsolidation ratio
    ocr = 2.0
    model.model_part.Properties[1].SetValue(OVERCONSOLIDATION_RATIO, ocr)
    k0 = 0.5

    ## setting the prestress and preconsolidation pressure
    for element in model.model_part.Elements:
        stresses = model_insitu.model_part.Elements[element.Id].GetValuesOnIntegrationPoints(STRESSES, model_insitu.model_part.ProcessInfo)
        prestresses = []
        preconsolidation_pressures = []
        for stress in stresses:
            #print("stress:", stress)
            prestress = []
            # for s in stress:
            #     prestress.append(-s)
            prestress = [-k0*stress[1], -stress[1], 0.0]
            prestresses.append(prestress)
            # preconsolidation_pressures.append(-ocr*(stress[0] + stress[1]) / 3)
            preconsolidation_pressures.append(-ocr*stress[1])
        # print("prestresses: " + str(prestresses))
        # print("preconsolidation_pressures: " + str(preconsolidation_pressures))
        element.SetValuesOnIntegrationPoints(PRESTRESS, prestresses, 3, model.model_part.ProcessInfo)
        element.SetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, preconsolidation_pressures, model.model_part.ProcessInfo)
        element.ResetConstitutiveLaw()
    # sys.exit(0)

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
    delta_load_lists = [1.0]
    for i in range(0, 48):
        delta_load_lists.append(0.02)
    for i in range(0, 4):
        delta_load_lists.append(0.01)
    for i in range(0, 4):
        delta_load_lists.append(0.005)
    for i in range(0, 2):
        delta_load_lists.append(0.0025)

    if logging:
        ifile = open("settlement_a.txt", "w")
        ifile.write("dy\tload_fact\n")
        # ifile.write("0.0\t0.0\n")

    load = 0.0
    first_disp = {}
    for delta_load in delta_load_lists:
        load = load + delta_load
        model.model_part.Properties[1].SetValue(DENSITY, 2.038735 * load)

        print(f"Load step {load} started")

        time = load
        model.Solve( time, 0, 0, 0, 0 )

        # # reset the displacement for first step
        # if load == 1.0:
        #     print(f"Reset displacements at load step {load}")
        #     for node in model.model_part.Nodes:
        #         node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        #         node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

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
            u_offset = pointA.GetSolutionStepValue(DISPLACEMENT_Y) - first_disp[pointA.Id][1]
            ifile.write("%.16e\t%.16e\n" % (u_offset, load))
            ifile.flush()

        print(f"Load step {load} completed")

    if logging:
        ifile.close()

    if last_output:
        model.WriteOutput( time )

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
    ref_disp = -6.3983441420996545e-01
    assert(abs(dy - ref_disp) / abs(ref_disp) < 1e-10)
    #####################
    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, last_output=False, logging=True, visual=False)
