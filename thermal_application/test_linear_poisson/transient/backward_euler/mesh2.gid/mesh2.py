##################################################################
import sys
import os
import math
import time as time_module
##################################################################
sys.path.append('./mesh2.gid')
import mesh2_include
from mesh2_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

# parameters
ref_temp = 0.0
length = 0.1
thick = 0.01

def main(logging=True, output=True):
    model = mesh2_include.Model('mesh2',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ## logging
    if logging:
        log_file = open("temperature.log", "w")
        log_file.write("time\ttemperature-far\ttemperature-near\n")

    ## boundary condition
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
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
        log_file.write("%.6e\t%.10e\t%.10e\n" % (time - time_init, monitor_node.GetSolutionStepValue(TEMPERATURE), monitor_node2.GetSolutionStepValue(TEMPERATURE)))

    ## transient analysis

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = False

    # delta_time = 1.0e-1
    delta_time = 1.0
    # for n in range(0, 1):
    # for n in range(0, 10):
    for n in range(0, 100):
    # for n in range(0, 1000):
        time += delta_time

        model.SolveModel(time)
        if output:
            model.WriteOutput(time*1e6)
        if logging:
            log_file.write("%.6e\t%.10e\t%.10e\n" % (time - time_init, monitor_node.GetSolutionStepValue(TEMPERATURE), monitor_node2.GetSolutionStepValue(TEMPERATURE)))

    ## reporting
    if logging:
        log_file.close()

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - length) < tol and abs(node.Y0 - thick/2) < tol:
            monitor_node = node
        if abs(node.X0 - length/20) < tol and abs(node.Y0 - thick/2) < tol:
            monitor_node2 = node

    temp_far = monitor_node.GetSolutionStepValue(TEMPERATURE)
    temp_near = monitor_node2.GetSolutionStepValue(TEMPERATURE)
    ref_temp_far = 4.5791836586273384e+01
    ref_temp_near = 4.9669831302281516e+01
    print("%.16e" % (temp_far))
    print("%.16e" % (temp_near))
    assert(abs(temp_far - ref_temp_far) < 1e-10)
    assert(abs(temp_near - ref_temp_near) < 1e-10)
    print("Test passed")

def tag():
    tags = "thermal,transient"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
