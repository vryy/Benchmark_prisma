##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mo 10. Jul 16:31:01 CEST 2023 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh2.gid')
import mesh2_include
from mesh2_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True, output_last=False):
    model = mesh2_include.Model('mesh2',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## set reference temperature
    ref_temp = 0.0

    length = 0.1
    thick = 0.01

    ## logging
    if logging:
        log_file = open("temperature.log", "w")
        log_file.write("%-*s%-*s%-*s%s\n" % (20, "time", 20, "temperature-far", 20, "temperature-near", "tip-displacement"))

    values = [ref_temp]*9
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(REFERENCE_TEMPERATURE, values, model.model_part.ProcessInfo)

    ## boundary condition
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            if abs(node.Y0) < tol:
                node.Fix(DISPLACEMENT_Y)
            node.Fix(TEMPERATURE)
            prescribed_nodes.append(node)
        if abs(node.X0 - length) < tol and abs(node.Y0 - thick/2) < tol:
            monitor_node = node
        if abs(node.X0 - length/20) < tol and abs(node.Y0 - thick/2) < tol:
            monitor_node2 = node
        node.SetSolutionStepValue(TEMPERATURE, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_NULL, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_EINS, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_DT, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_NULL_DT, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_EINS_DT, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION, 0.0)

    ## analysis step 0

    time_init = 1.0e-6
    time = time_init
    for node in prescribed_nodes:
        node.SetSolutionStepValue(TEMPERATURE, ref_temp + 50.0)
        node.SetSolutionStepValue(TEMPERATURE_NULL, ref_temp + 50.0)
        node.SetSolutionStepValue(TEMPERATURE_EINS, ref_temp + 50.0)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = True

    model.SolveModel(time)
    if output:
        model.WriteOutput(time*1e6)
    if logging:
        log_file.write("%-*.6e%-*.10e%-*.10e%.10e\n" % (20, time - time_init, 20, monitor_node.GetSolutionStepValue(TEMPERATURE), 20, monitor_node2.GetSolutionStepValue(TEMPERATURE), monitor_node.GetSolutionStepValue(DISPLACEMENT_Y)))
    # sys.exit(0)
    ## transient analysis

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = False

    delta_time = 1.0e-1
    #for n in range(0, 1):
    for n in range(0, 1000):
        time += delta_time

        model.SolveModel(time)
        if output:
            model.WriteOutput(time*1e6)
        if logging:
            log_file.write("%-*.6e%-*.10e%-*.10e%.10e\n" % (20, time - time_init, 20, monitor_node.GetSolutionStepValue(TEMPERATURE), 20, monitor_node2.GetSolutionStepValue(TEMPERATURE), monitor_node.GetSolutionStepValue(DISPLACEMENT_Y)))

    ## reporting

    if logging:
        log_file.close()

    if output_last and not output:
        model.WriteOutput(time*1e6)

    return model

def test():
    model = main(output=False, logging=False, output_last=False)

    length = 0.1
    thick = 0.01

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - length) < tol and abs(node.Y0 - thick/2) < tol:
            monitor_node = node
        if abs(node.X0 - length/20) < tol and abs(node.Y0 - thick/2) < tol:
            monitor_node2 = node

    near_temp_ref = 4.520028576441864e+01
    far_temp_ref = 4.962341876085018e+01
    tip_disp = 2.392800357563666e-04
    # print("near temp: %.15e" % (monitor_node.GetSolutionStepValue(TEMPERATURE)))
    assert((monitor_node.GetSolutionStepValue(TEMPERATURE) - near_temp_ref) / near_temp_ref < 1e-10)
    # print("far temp: %.15e" % (monitor_node2.GetSolutionStepValue(TEMPERATURE)))
    assert((monitor_node2.GetSolutionStepValue(TEMPERATURE) - 4.9623418761e+01) / far_temp_ref < 1e-10)
    # print("tip disp: %.15e" % (monitor_node.GetSolutionStepValue(DISPLACEMENT_Y)))
    assert((monitor_node2.GetSolutionStepValue(DISPLACEMENT_Y) - tip_disp) < 1e-10)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=False, logging=True, output_last=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
