##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import dewatering_q9_include as simulation_include
try:
    from dewatering_q9_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'dewatering_q9'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

aux_util = SoilsAuxiliaryUtility()

def SetMaterialProperties(elem, poro_calc_mode=0):
    elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
    elem.SetValue(DENSITY,                  2000.0)
    elem.SetValue(DENSITY_WATER,            1.0e3)
    elem.SetValue(POROSITY,                 0.2975)
    elem.SetValue(PERMEABILITY_WATER,       4.4e-6)
    elem.SetValue(PERMEABILITY_AIR,         3.2e-7)
    aux_util.SetValue(SWCC_LAW, LiakopolousSWCC(), elem)
    aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, LiakopolousRelativePermeabilityWaterLaw(), elem)
    aux_util.SetValue(RELATIVE_PERMEABILITY_AIR_LAW, LiakopolousRelativePermeabilityAirLaw(), elem)
    aux_util.SetValue(GAS_LAW, IdealGasLaw(1.295, 1.188280000e-05), elem)
    elem.SetValue(POROSITY_CALCULATION_MODE, poro_calc_mode)

def main(output=True, logging=True, fix_lateral=False, poro_calc_mode=0, testing=False, damping_alpha=0.0):

    model_virgin = simulation_include.Model(model_name_,os.getcwd()+"/",os.getcwd()+"/virgin_results/",logging=False)
    model_virgin.InitializeModel()

    # material parameters
    for e in model_virgin.layer_sets['Layer0']:
        elem = model_virgin.model_part.Elements[e]
        SetMaterialProperties(elem, poro_calc_mode=poro_calc_mode)
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

        ## uniform air flow
        node.Fix(AIR_PRESSURE)
        node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

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
        if output:
            model_virgin.WriteOutput(time)

    isu = InSituStressUtility()
    isu.SetPreStressFromCurrentStress(model_virgin.model_part, model_virgin.model_part.ProcessInfo)
    for elem in model_virgin.model_part.Elements:
        elem.ResetConstitutiveLaw()

    time = time + delta_time
    model_virgin.SolveModel(time)
    if output:
        model_virgin.WriteOutput(time)

    max_disp = 0.0
    for node in model_virgin.model_part.Nodes:
        for direction in range(0,3):
            if( abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction])) > max_disp ):
                max_disp = abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction]))

    print("~~ STEP DONE (INSITU STRESS) --> residual displacement= "+str(max_disp)+"~~")
    # sys.exit(0)

    ## solve the system

    model = simulation_include.Model(model_name_,os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    model.model_part.Properties[1].SetValue(RAYLEIGH_DAMPING_ALPHA, damping_alpha )
    model.model_part.Properties[1].SetValue(RAYLEIGH_DAMPING_BETA, 0.0 )

    # material parameters
    for e in model.layer_sets['Layer0']:
        elem = model.model_part.Elements[e]
        SetMaterialProperties(elem, poro_calc_mode=poro_calc_mode)
        elem.Initialize(model.model_part.ProcessInfo)

    # boundary condition
    tol = 1.0e-6
    top_nodes = []
    bottom_nodes = []
    lateral_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            lateral_nodes.append(node)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            bottom_nodes.append(node)
        if abs(node.Y0 - 1.0) < tol:
            top_nodes.append(node)

        node.Fix(WATER_PRESSURE)
        node.Fix(AIR_PRESSURE)

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

    # release the water and air pressure on the model
    for node in model.model_part.Nodes:
        node.Free(WATER_PRESSURE)
        node.Free(AIR_PRESSURE)

    # but fix air pressure at bottom and lateral
    # we do not fix air pressure on top create a smooth distribution of AIR_PRESSURE in the center
    # In terms of the experiment, it's in line with closing the water and air flow in on top
    for node in bottom_nodes:
        node.Fix(AIR_PRESSURE)
        node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)
    if fix_lateral:
        for node in lateral_nodes:
            node.Fix(AIR_PRESSURE)
            node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

    # but fix water pressure on bottom
    for node in bottom_nodes:
        node.Fix(WATER_PRESSURE)
        node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    # delta_time = 0.000001
    # ## TESTING
    # for step in range(0,2):
    #     time = time + delta_time
    #     model.Solve( time, 0, 0, 0, 0 )
    #     # model.WriteOutput( time )
    #     Record(ifile, time)
    #     print("delta_time: " + str(delta_time))
    # sys.exit(0)
    # ## END TESTING
    ## REMARK: When the time step is small, the contribution from mass matrix is very large (multiply by (1-alphaf)/(beta*Dt*Dt)).
    ## Thefore, the small time steps shall be skipped for Generalized-alpha.
    # for step in range(0,10):
    #     time = time + delta_time
    #     model.Solve( time, 0, 0, 0, 0 )
    #     # model.WriteOutput( time )
    #     Record(ifile, time)
    #     print("delta_time: " + str(delta_time))
    # delta_time = 0.00001
    # for step in range(0,11):
    #     time = time + delta_time
    #     model.Solve( time, 0, 0, 0, 0 )
    #     # model.WriteOutput( time )
    #     Record(ifile, time)
    #     print("delta_time: " + str(delta_time))
    # delta_time = 0.0001
    # for step in range(0,11):
    #     time = time + delta_time
    #     model.Solve( time, 0, 0, 0, 0 )
    #     # model.WriteOutput( time )
    #     Record(ifile, time)
    #     print("delta_time: " + str(delta_time))
    # delta_time = 0.001
    # for step in range(0,11):
    #     time = time + delta_time
    #     model.Solve( time, 0, 0, 0, 0 )
    #     # model.WriteOutput( time )
    #     Record(ifile, time)
    #     print("delta_time: " + str(delta_time))
    # delta_time = 0.01
    # for step in range(0,11):
    #     time = time + delta_time
    #     model.Solve( time, 0, 0, 0, 0 )
    #     # model.WriteOutput( time )
    #     Record(ifile, time)
    #     print("delta_time: " + str(delta_time))
    delta_time = 0.1
    for step in range(0,31):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
        print("delta_time: " + str(delta_time))
    if testing:
        return model
    delta_time = 1.0
    for step in range(0,51):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
        print("delta_time: " + str(delta_time))
    delta_time = 10.0
    for step in range(0,51):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
        print("delta_time: " + str(delta_time))
    delta_time = 100.0
    for step in range(0,71):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
        if logging:
            Record(ifile, time)
        print("delta_time: " + str(delta_time))
    print(model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_Z))

    if logging:
        ifile.close()

def test():
    model = main(logging=False, output=False, fix_lateral=True, poro_calc_mode=0, testing=True, damping_alpha=0.01)

    tol = 1.0e-6
    bottom_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Z0 - 0.0) < tol:
            bottom_nodes.append(node)

    vtu = VariableTransferUtility(SuperLUSolver())
    vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)
    wf = bottom_nodes[0].GetSolutionStepValue(WATER_FLOW)
    wf2 = wf[1]*60.0/0.01
    print("%.16e" % (wf2))
    ref_wf2 = -2.6399967335706660e-02
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
        main(logging=True, output=False, fix_lateral=True, poro_calc_mode=0, damping_alpha=0.01)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
