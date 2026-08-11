##################################################################
## 1D bar heat transfer example
## Ref:
## + Chinesta, Model Order Reduction Based on Proper Orthogonal Decomposition, Section 2.3
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import bar_include as simulation_include
try:
    from bar_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'bar'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

length = 1.0
thick = 0.01

def main(logging=True, output=True, snapshot=False, bc_type=1, pod="none", decouple_build_and_solve=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging, decouple_build_and_solve=decouple_build_and_solve)
    model.InitializeModel()

    if pod == "load":
        pod_process = PodProjectionProcess("Q.bin")
        model.solver.solver.builder_and_solver.SetPodProcess(pod_process)
        if decouple_build_and_solve:
            model.solver.solver.linear_solver = pod_process
            # this must be set if decouple_build_and_solve == True
            # since pod_process is also a linear solver, this is usable and necessary

    ## initial condition
    tol = 1e-6
    ref_temp = 1.0
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Y0 - thick) < tol:
            monitor_node = node
        node.SetSolutionStepValue(TEMPERATURE, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_NULL, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_EINS, ref_temp)
        node.SetSolutionStepValue(TEMPERATURE_DT, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_NULL_DT, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_EINS_DT, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION, 0.0)
        node.SetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION, 0.0)

    ## boundary condition
    lbd = model.model_part.Properties[1].GetValue(THERMAL_CONDUCTIVITY)
    if bc_type == 1:
        flux = -1.0 # heat flux goes inward
        model.model_part.Conditions[1].SetValue(HEAT_FLUX, lbd*flux)
            # we have to multiply the flux with lambda because the Green formula leads to
            # boundary integral with coefficient lambda
        model.model_part.Conditions[2].SetValue(HEAT_FLUX, 0.0)
    elif bc_type == 2:
        model.model_part.Conditions[1].SetValue(HEAT_FLUX, 0.0)
        model.model_part.Conditions[2].SetValue(HEAT_FLUX, 0.0)

    ## logging
    if logging:
        log_file = open("temperature.log", "w")
        log_file.write("time\tnode_%d\n" % (monitor_node.Id))

    time = 0.0
    if logging:
        log_file.write("%.6e\t%.16e\n" % (time, monitor_node.GetSolutionStepValue(TEMPERATURE)))

    ## analysis step 0

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = True

    delta_time = 1.0e-1
    time += delta_time
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    if logging:
        log_file.write("%.6e\t%.16e\n" % (time, monitor_node.GetSolutionStepValue(TEMPERATURE)))

    ## transient analysis

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = False

    marked_time = 1.0
    for n in range(0, 300):
        time += delta_time

        if bc_type == 1:
            if time > 10.0:
                model.model_part.Conditions[1].SetValue(HEAT_FLUX, 0.0)
        elif bc_type == 2:
            if time <= 20.0:
                flux = -time/20.0
            else:
                flux = -(time-30.0)/5.0
            model.model_part.Conditions[1].SetValue(HEAT_FLUX, lbd*flux)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            log_file.write("%.6e\t%.16e\n" % (time, monitor_node.GetSolutionStepValue(TEMPERATURE)))

        if abs(time - marked_time) < 1e-6:
            if snapshot:
                ifile2 = open("snapshot.%d" % (marked_time), "w")
                ifile2.write("x\ttemp\n")
                for node in model.model_part.Nodes:
                    if abs(node.Y0 - thick) < tol:
                        ifile2.write("%.6e\t%.16e\n" % (node.X0, node.GetSolutionStepValue(TEMPERATURE)))
                ifile2.close()
            marked_time += 1.0

    if pod == "save":
        model.solver.solver.builder_and_solver.GetPodProcess().SavePrincipalComponents("Q.bin", 5)

    ## reporting
    if logging:
        log_file.close()

    return model

def test1():
    model = main(logging=False, output=False, snapshot=False, bc_type=1, pod="save")

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Y0 - thick) < tol:
            monitor_node = node

    temp = monitor_node.GetSolutionStepValue(TEMPERATURE)
    print("%.16e" % (temp))
    ref_temp = 1.1176244923232193e+00
    assert(abs(temp - ref_temp) < 1e-10)
    print("Test passed")

def test2():
    model = main(logging=False, output=False, snapshot=False, bc_type=1, pod="load", decouple_build_and_solve=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Y0 - thick) < tol:
            monitor_node = node

    temp = monitor_node.GetSolutionStepValue(TEMPERATURE)
    print("%.16e" % (temp))
    ref_temp = 1.1175733587465031e+00
    assert(abs(temp - ref_temp) < 1e-10)
    print("Test passed")

def test3(): # solve a problem with different BC using the same POD basis
    model = main(logging=False, output=False, snapshot=False, bc_type=2, pod="load", decouple_build_and_solve=True)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and abs(node.Y0 - thick) < tol:
            monitor_node = node

    temp = monitor_node.GetSolutionStepValue(TEMPERATURE)
    print("%.16e" % (temp))
    ref_temp = 9.0681192567052238e-01
    assert(abs(temp - ref_temp) < 1e-10)
    print("Test passed")

def test():
    test1()
    test2()
    test3()

def gen_snapshot_1():
    # solve the FOM with BC 1
    main(logging=False, output=False, snapshot=True, bc_type=1, pod="none")

def gen_snapshot_1_pod():
    # solve the ROM with BC 1
    main(logging=False, output=False, snapshot=True, bc_type=1, pod="load")

def gen_snapshot_2():
    # solve the FOM with BC 2
    main(logging=False, output=False, snapshot=True, bc_type=2, pod="none")

def gen_snapshot_2_pod():
    # solve the ROM with BC 2 using the POD basis of problem with BC 1
    main(logging=False, output=False, snapshot=True, bc_type=2, pod="load")

def tag():
    tags = "unknown"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, snapshot=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
