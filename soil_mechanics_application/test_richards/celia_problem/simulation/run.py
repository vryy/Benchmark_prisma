##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025 by Hoang-Giang Bui (UoB)                 #####
#####          2026 by Hoang-Giang Bui (DU)                  #####
##### all rights reserved                                    #####
##################################################################
# Reference:
# + geometry: Sepulveda-Cano et al, Numerical analysis of soil desaturation by an air injection method
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import simulation_include
try:
    from simulation_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
# model_name_ = 'ground1'
model_name_ = 'ground'
# model_name_ = 'ground_h20'
# model_name_ = 'ground_h27'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################
## parameters extracted from Celia paper:
## A General Mass-Conservative Numerical Solution for the Unsaturated Flow Equation, Water Resources Research
alpha = 1.611e6
theta_s = 0.287
theta_r = 0.075
beta = 3.96
Ks = 0.00944
A = 1.175e6
gamma = 4.74
##################################################################

aux_util = SoilsAuxiliaryUtility()
vtu = VariableTransferUtility(SuperLUSolver())

Smax = theta_s
Smin = theta_r
pr = math.pow(alpha, 1.0/beta)
s1 = beta
s2 = 1.0

def SetMaterialProperties(elem):
    elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
    elem.SetValue(PERMEABILITY_WATER,       Ks)
    aux_util.SetValue(SWCC_LAW, VanGenuchtenSWCC(s1, s2, pr, Smin, Smax), elem)
    aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, HaverkampRelativePermeabilityWaterHeadLaw(A, gamma), elem)

def main(output=True, logging=True, total_time=360.0, delta_time=1.0, \
    solution_strategy="implicit_Newton_Raphson", analysis_type=2, dissipation_radius=0.9):

    ## solve the system

    model = simulation_include.Model(model_name_,os.getcwd()+"/",os.getcwd()+"/", \
        logging=logging, \
        solution_strategy=solution_strategy, \
        analysis_type=analysis_type, \
        dissipation_radius=dissipation_radius, \
        stop_Newton_Raphson_if_not_converge=True)
    model.InitializeModel()

    # material parameters
    for e in model.layer_sets['Ground']:
        elem = model.model_part.Elements[e]
        SetMaterialProperties(elem)
        elem.Initialize(model.model_part.ProcessInfo)

    # model.WriteOutput(0.0)
    # sys.exit(0)

    # boundary condition
    tol = 1.0e-6
    top_nodes = []
    bottom_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Z0 - 0.0) < tol:
            bottom_nodes.append(node)
        if abs(node.Z0 - 40.0) < tol:
            top_nodes.append(node)

        node.Fix(WATER_PRESSURE)

    h_bottom = 61.5
    h_top = 20.7

    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(WATER_PRESSURE, h_bottom)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, h_bottom)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, h_bottom)

    for node in top_nodes:
        node.SetSolutionStepValue(WATER_PRESSURE, h_top)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, h_top)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, h_bottom)

    # for node in model.model_part.Nodes:
    #     a = node.Z0 / 40.0
    #     p = (1-a)*p_bottom + a*p_top
    #     node.SetSolutionStepValue(WATER_PRESSURE, p)
    #     node.SetSolutionStepValue(WATER_PRESSURE_EINS, p)
    #     node.SetSolutionStepValue(WATER_PRESSURE_NULL, p)

    ## reset the material one more time to account for new information
    for element in model.model_part.Elements:
        element.ResetConstitutiveLaw()

    # release the water pressure on the model
    for node in model.model_part.Nodes:
        node.Free(WATER_PRESSURE)

    # but fix water pressure on top and bottom
    for node in top_nodes + bottom_nodes:
        node.Fix(WATER_PRESSURE)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    time = 0.0
    nsteps = int(total_time / delta_time)
    for i in range(0, nsteps):

        time += delta_time

        print(f"### Time step {time} s", flush=True)

        model.SolveModel(time)

        vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)

        if output:
            model.WriteOutput(time)

        model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

        # for node in top_nodes:
        #     node.SetSolutionStepValue(WATER_PRESSURE, h_top)
        #     node.SetSolutionStepValue(WATER_PRESSURE_EINS, h_top)
        #     node.SetSolutionStepValue(WATER_PRESSURE_NULL, h_bottom)

    if logging:
        ifile = open("pressure_head.txt", "w")
        ifile.write("depth\tpressure_head\twater_flow\n")

        values = []
        for node in model.model_part.Nodes:
            if abs(node.X0) < tol and abs(node.Y0) < tol:
                h = node.GetSolutionStepValue(WATER_PRESSURE)
                f = node.GetSolutionStepValue(WATER_FLOW)
                values.append([node.Z0, h, f])

        sorted_values = sorted(values, key=lambda item: item[0])
        for v in sorted_values:
            ifile.write("%.16e\t%.16e\t%.16e\n" % (v[0], v[1], v[2][2]))

        ifile.close()

    return model

def test():
    model = main(logging=False, output=False, total_time=360.0, delta_time=1.0, analysis_type=4, dissipation_radius=1.0)

    node = model.model_part.Nodes[33]

    p = node.GetSolutionStepValue(WATER_PRESSURE)
    print("%.16e" % (p))
    ref_p = 3.6619148488141356e+01
    assert(abs(p - ref_p) < 1e-10)

    print("Test passed")

def tag():
    if all_modules_are_imported_successfully:
        return "richards"
    else:
        return ""

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        # main(logging=True, output=False, total_time=360.0, delta_time=1.0, analysis_type=4, dissipation_radius=0.5)
        # main(logging=True, output=False, total_time=360.0, delta_time=1.0, analysis_type=4, dissipation_radius=0.5)
        # main(logging=True, output=False, total_time=360.0, delta_time=0.5, analysis_type=4, dissipation_radius=0.5)
        main(logging=True, output=False, total_time=360.0, delta_time=1.0, analysis_type=4, dissipation_radius=1.0)
        # main(logging=True, output=False, total_time=1.0, delta_time=1.0, analysis_type=4, dissipation_radius=1.0)
        # main(logging=True, output=False, total_time=360.0, delta_time=1.0, analysis_type=2, dissipation_radius=0.9)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
