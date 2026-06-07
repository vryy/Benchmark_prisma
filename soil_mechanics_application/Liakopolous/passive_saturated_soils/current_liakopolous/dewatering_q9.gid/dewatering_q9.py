##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import dewatering_q9_include
from dewatering_q9_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):

    model = dewatering_q9_include.Model('dewatering_q9',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    aux_util = SoilsAuxiliaryUtility()

    model_virgin = dewatering_q9_include.Model('dewatering_q9',os.getcwd()+"/",os.getcwd()+"/virgin_results/",logging=False)
    model_virgin.InitializeModel()

    # print(HARDENING_LAW)
    # print(ISOTROPIC_HARDENING_LAW)
    # print(KINEMATIC_HARDENING_LAW)
    # print(SWCC_LAW)
    # print(RELATIVE_PERMEABILITY_WATER_LAW)

    # material parameters
    for e in model_virgin.layer_sets['Layer0']:
        elem = model_virgin.model_part.Elements[e]
        elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
        elem.SetValue(DENSITY,                  1500.0)
        elem.SetValue(DENSITY_WATER,            1.0e3)
        # elem.SetValue(DENSITY_AIR,              1.285)
        # elem.SetValue(BULK_AIR,                 1.01e5)
        elem.SetValue(POROSITY,                 0.2975)
        elem.SetValue(PERMEABILITY_WATER,       4.41e-6) #0.0000044)
        # elem.SetValue(PERMEABILITY_AIR,         0.00000032)
        # elem.SetValue(FIRST_SATURATION_PARAM,   2.5)
        # elem.SetValue(SECOND_SATURATION_PARAM,  0.4)
        # elem.SetValue(AIR_ENTRY_VALUE,          3000.0)
        # aux_util.SetValue(SWCC_LAW, VanGenuchtenSWCC(2.5, 0.4, 3000.0), elem)
        # aux_util.SetValue(SWCC_LAW, FullySaturatedSWCC(), elem)
        aux_util.SetValue(SWCC_LAW, LiakopolousSWCC(), elem)
        # aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, SimpleRelativePermeabilityLaw(), elem)
        aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, LiakopolousRelativePermeabilityWaterLaw(), elem)
        elem.SetValue(POROSITY_CALCULATION_MODE,  2)

        elem.Initialize(model_virgin.model_part.ProcessInfo)

    # boundary condition
    tol = 1.0e-6
    bottom_nodes = []
    for node in model_virgin.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            bottom_nodes.append(node)

        ## uniform water flow
        node.Fix(WATER_PRESSURE)
        node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)

    nsteps = 1
    time = 0.0
    delta_time = 100.0
    for step in range(0, nsteps):
        gravity = Vector(3)
        gravity[0] = 0.0
        gravity[1] = -9.81*(step + 1) / nsteps
        gravity[2] = 0.0

        model_virgin.model_part.Properties[1].SetValue(GRAVITY, gravity)

        time = time + delta_time
        model_virgin.SolveModel(time)
        model_virgin.WriteOutput(time)

    isu = InSituStressUtility()
    isu.SetPreStressFromCurrentStress(model_virgin.model_part, model_virgin.model_part.ProcessInfo)
    for elem in model_virgin.model_part.Elements:
        elem.ResetConstitutiveLaw()

    time = time + delta_time
    model_virgin.SolveModel(time)
    model_virgin.WriteOutput(time)

    max_disp = 0.0
    for node in model_virgin.model_part.Nodes:
        for direction in range(0,3):
            if( abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction])) > max_disp ):
                max_disp = abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction]))

    print("~~ STEP DONE (INSITU STRESS) --> residual displacement= "+str(max_disp)+"~~")

    ## solve the system

    # material parameters
    for e in model.layer_sets['Layer0']:
        elem = model.model_part.Elements[e]
        elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
        elem.SetValue(DENSITY,                  1500.0)
        elem.SetValue(DENSITY_WATER,            1.0e3)
        # elem.SetValue(DENSITY_AIR,              1.285)
        # elem.SetValue(BULK_AIR,                 1.01e5)
        elem.SetValue(POROSITY,                 0.2975)
        elem.SetValue(PERMEABILITY_WATER,       4.41e-6) #0.0000044)
        # elem.SetValue(PERMEABILITY_AIR,         0.00000032)
        # elem.SetValue(FIRST_SATURATION_PARAM,   2.5)
        # elem.SetValue(SECOND_SATURATION_PARAM,  0.4)
        # elem.SetValue(AIR_ENTRY_VALUE,          3000.0)
        # aux_util.SetValue(SWCC_LAW, VanGenuchtenSWCC(2.5, 0.4, 3000.0), elem)
        # aux_util.SetValue(SWCC_LAW, FullySaturatedSWCC(), elem)
        aux_util.SetValue(SWCC_LAW, LiakopolousSWCC(), elem)
        # aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, SimpleRelativePermeabilityLaw(), elem)
        aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, LiakopolousRelativePermeabilityWaterLaw(), elem)
        elem.SetValue(POROSITY_CALCULATION_MODE,  2)

        elem.Initialize(model.model_part.ProcessInfo)

    # boundary condition
    tol = 1.0e-6
    top_nodes = []
    bottom_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            bottom_nodes.append(node)
        if abs(node.Y0 - 1.0) < tol:
            top_nodes.append(node)

        node.Fix(WATER_PRESSURE)

    # transfer insitu stress
    vtu = VariableTransferUtility(SuperLUSolver())
    vtu.TransferPrestressIdentically( model_virgin.model_part, model.model_part )

    def Record(ifile, time):
        vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)
        wf = bottom_nodes[0].GetSolutionStepValue(WATER_FLOW)
        ifile.write("%.10e\t%.10e\n" % (time, wf[1]))
        ifile.flush()

    if logging:
        ifile = open("bottom_water_flow.log", "w")
        ifile.write("time\tflow\n")

    ## reset the material one more time to account for new information
    for element in model.model_part.Elements:
        element.ResetConstitutiveLaw()

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    time = 0.0
    delta_time = 0.000001
    # delta_time = 100.0

    # first solve to get equilibrium
    nsteps = 1
    for step in range(0, nsteps):
        time = time + delta_time
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            Record(ifile, time)
    max_disp = 0.0
    for node in model.model_part.Nodes:
        for direction in range(0,3):
            if( abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction])) > max_disp ):
                max_disp = abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction]))

    print(f"~~ STEP DONE (APPLICATION OF INSITU STRESS) --> residual displacement= {max_disp} ~~")
    # sys.exit(0)

    # release the water pressure on the model
    for node in model.model_part.Nodes:
        node.Free(WATER_PRESSURE)

    # # but fix water pressure on top
    # for node in top_nodes:
    #     node.Fix(WATER_PRESSURE)
    #     node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
    #     node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
    #     node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)

    # but fix water pressure on bottom
    for node in bottom_nodes:
        node.Fix(WATER_PRESSURE)
        node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    delta_time = 0.000001
    for step in range(0,10):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 0.00001
    for step in range(0,11):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 0.0001
    for step in range(0,11):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 0.001
    for step in range(0,11):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 0.01
    for step in range(0,11):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 0.1
    for step in range(0,31):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 1.0
    for step in range(0,51):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 10.0
    for step in range(0,51):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    delta_time = 100.0
    for step in range(0,71):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
    # print(model.model_part.Nodes[217].GetSolutionStepValue(DISPLACEMENT_Y))

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1.0e-6
    bottom_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Y0 - 0.0) < tol:
            bottom_nodes.append(node)

    vtu = VariableTransferUtility(SuperLUSolver())
    vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)
    wf = bottom_nodes[0].GetSolutionStepValue(WATER_FLOW)
    wf2 = wf[1]*60.0/0.01
    print("%.16e" % (wf2))
    ref_wf2 = -1.7730811053657796e-03
    assert(abs(wf2 - ref_wf2) < 1e-10)

    print("Test passed")

def tag():
    return "partially-saturated"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
