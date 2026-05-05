##################################################################
# This example tests the tensile part of C-DPM2. The stress will
# reach and bounded by ft
# Reference; _MAT_LAW124 (CPDM2).pdf
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import cube1_include
from cube1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True, integration_order=1):
    model = cube1_include.Model('cube1',current_dir_,current_dir_,logging=logging)
    model.InitializeModel(integration_order=integration_order)

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)

    time = 0.0
    model.SolveModel(time)

    if logging:
        ifile = open("reaction.log", "w")
        ifile.write("%-*s%-*s%s\n" % (20, "disp", 20, "x-reaction", "sigma-xx"))

    disp = 0.0
    delta_disp = 1e-6
    delta_time = 1.0

    for i in range(0, 2000):
        print("#######################################")
        print(f"##########LOADING STEP {i+1} STARTED##############")
        print("#######################################")

        time += delta_time
        disp += delta_disp
        for node in prescribed_nodes:
            # node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            react_p = 0.0
            for node in prescribed_nodes:
                react_p += node.GetSolutionStepValue(REACTION_X)
            stress = model.model_part.Elements[1].CalculateOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
            ifile.write("%-*.10e%-*.10e%.10e\n" % (20, disp, 20, react_p, stress[0][0]))
            ifile.flush()

    if logging:
        ifile.close()

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    return model

def test():
    model = main(logging=False, output=False, integration_order=2)

    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            prescribed_nodes.append(node)

    react_p = 0.0
    for node in prescribed_nodes:
        react_p += node.GetSolutionStepValue(REACTION_X)
    print("%.16e" % (react_p))

    ref_reac = 7.1047415546459597e+03
    assert(abs(react_p - ref_reac) / abs(ref_reac) < 1e-10)

    print("Test passed")

def tag():
    return "oofem,concrete,cdpm2"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, integration_order=2)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
