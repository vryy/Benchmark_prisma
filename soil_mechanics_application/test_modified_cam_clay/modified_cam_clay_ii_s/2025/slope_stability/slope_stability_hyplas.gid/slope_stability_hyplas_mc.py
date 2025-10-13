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
### This test is to make sure MC works with prestress
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

    model = slope_stability_hyplas_include.Model('slope_stability_hyplas',os.getcwd()+"/",logging=logging,rel_tol=1e-8,abs_tol=1e-6)
    model.material = "mohr-coulomb"
    model.InitializeModel()

    ## set the overconsolidation ratio
    ocr = 2.0
    model.model_part.Properties[1].SetValue(OVERCONSOLIDATION_RATIO, ocr)
    k0 = 0.5

    ## setting the prestress
    for element in model.model_part.Elements:
        stresses = model_insitu.model_part.Elements[element.Id].GetValuesOnIntegrationPoints(STRESSES, model_insitu.model_part.ProcessInfo)
        prestresses = []
        for stress in stresses:
            #print("stress:", stress)
            prestress = []
            # for s in stress:
            #     prestress.append(-s)
            prestress = [-k0*stress[1], -stress[1], 0.0]
            prestresses.append(prestress)
        # print("prestresses: " + str(prestresses))
        element.SetValuesOnIntegrationPoints(PRESTRESS, prestresses, 3, model.model_part.ProcessInfo)
        element.ResetConstitutiveLaw()
    # model.WriteOutput(0.0)
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
    delta_load_lists = [0.5, 0.5, 1.0, 1.0/2, 1.0/8, 1.0/16, 1.0/32] # sum = 2.71875

    if logging:
        ifile = open("settlement_a-mc.txt", "w")
        ifile.write("dy\tload_fact\n")
        # ifile.write("0.0\t0.0\n")

    load = 0.0
    first_disp = {}
    for delta_load in delta_load_lists:
        load = load + delta_load
        model.model_part.Properties[1].SetValue(DENSITY, 2.038735 * load)

        # here we have to ramp up to full gravity. Because applying load factor 1.0 immediately will not converge.
        if load < 1.0:
            prestress_factor = [load]*4
            for element in model.model_part.Elements:
                element.SetValuesOnIntegrationPoints(PRESTRESS_FACTOR, prestress_factor, model.model_part.ProcessInfo)
                element.ResetConstitutiveLaw()

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
            if load >= 1.0:
                u = pointA.GetSolutionStepValue(DISPLACEMENT_Y) - first_disp[pointA.Id][1]
                ifile.write("%.16e\t%.16e\n" % (u, load))
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
    ref_disp = -1.5565077819123057e+00
    assert(abs(dy - ref_disp) / abs(ref_disp) < 1e-10)
    #####################
    print("Test passed")

def tag():
    return "mohr-coulomb,prestress"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=False, last_output=True, logging=True, visual=False)
