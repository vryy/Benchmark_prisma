##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import core_h27_include
from core_h27_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, time, nodes):
    ifile.write("%.6e" % (time))
    for node in nodes:
        ifile.write("\t%.10e" % (node.GetSolutionStepValue(TEMPERATURE_EINS)))
    ifile.write("\n")
    ifile.flush()

def main(output=True, output_last=False, logging=True, nsteps=3600, delta_time = 1.):
    model = core_h27_include.Model('core_h27',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## boundary conditions
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
        if abs(node.Z0 - 0.06) < tol:
            prescribed_nodes.append(node)

    # set the initial density
    initial_density = 8400.0
    values = [initial_density]*27
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(DENSITY, values, model.model_part.ProcessInfo)

    # set the initial temperature
    ref_temp = 293.15 # 20 oC

    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(TEMPERATURE, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_NULL, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_EINS, ref_temp)

    model.model_part.Properties[1].SetValue(REFERENCE_TEMPERATURE, ref_temp )
    model.model_part.Properties[2].SetValue(REFERENCE_TEMPERATURE, ref_temp )

    # set the reference temperature
    values = [ref_temp]*27
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(REFERENCE_TEMPERATURE, values, model.model_part.ProcessInfo)

    # set the initial pressure
    initial_air_pressure = 1e5
    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(AIR_PRESSURE, initial_air_pressure) # Pa
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, initial_air_pressure) # Pa
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, initial_air_pressure) # Pa

    # set the top pressure
    top_air_pressure = 1e6 # 10 bar
    for node in prescribed_nodes:
        node.Fix(AIR_PRESSURE)
        node.SetSolutionStepValue(AIR_PRESSURE, top_air_pressure)
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, top_air_pressure)
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, top_air_pressure)

    #################################

    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.06) < tol):
            pointT = node
        if (abs(node.X0 - 0.015) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.035) < tol):
            pointA = node
        if (abs(node.X0 - 0.015) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.025) < tol):
            pointB = node
        if (abs(node.X0 - 0.015) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.015) < tol):
            pointC = node

    if logging:
        monitoring_nodes = [pointA, pointB, pointC, pointT]

        ifile = open("monitoring.log", "w")
        ifile.write("time\ttemp_A\ttemp_B\ttemp_C\ttemp_T\n")

    #################################
    time = 0.
    if output:
        model.WriteOutput(time)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = True

    time = delta_time
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)
    if logging:
        WriteLog(ifile, time, monitoring_nodes)
    # sys.exit(0)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = False

    for i in range(0, nsteps):
        time += delta_time
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            WriteLog(ifile, time, monitoring_nodes)

    if output_last:
        model.WriteOutput(time)
    if logging:
        ifile.close()

    return model

def test():

    model = main(output=False, logging=False, nsteps=9)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.06) < tol):
            pointT = node
        if (abs(node.X0 - 0.015) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.035) < tol):
            pointA = node
        if (abs(node.X0 - 0.015) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.025) < tol):
            pointB = node
        if (abs(node.X0 - 0.015) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0 - 0.015) < tol):
            pointC = node

    temp_A = pointA.GetSolutionStepValue(TEMPERATURE_EINS)
    temp_B = pointB.GetSolutionStepValue(TEMPERATURE_EINS)
    temp_C = pointC.GetSolutionStepValue(TEMPERATURE_EINS)

    print("temp_A: %.16e" % temp_A)
    print("temp_B: %.16e" % temp_B)
    print("temp_C: %.16e" % temp_C)

    ref_temp_A = 3.4754751532530639e+02
    ref_temp_B = 3.4754739470578807e+02
    ref_temp_C = 3.4760735573635804e+02

    assert(abs(temp_A - ref_temp_A) < 1e-10)
    assert(abs(temp_B - ref_temp_B) < 1e-10)
    assert(abs(temp_C - ref_temp_C) < 1e-10)

    print("Test passed")

def tag():
    return "MH"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, output_last=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
