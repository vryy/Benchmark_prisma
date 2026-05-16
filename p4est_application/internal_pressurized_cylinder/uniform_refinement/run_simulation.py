##################################################################
## Convergence study of the internal pressurized cylinder problem
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
include_path = os.getcwd() + '/../'
sys.path.append(include_path)
import meshxx_p4est_include
from meshxx_p4est_include import *

import simulator
import p4est_simulator
import analytical_solution

##################################################################
###  SIMULATION  #################################################
##################################################################

start = time_module.time()

def main(name="mesh10x10", logging=True, output=True,nrefine=0,stress_recovery_type=0,neighbour_expansion_level=1):

    mdpa_path = os.getcwd() + '/../design_data/' + name + '.gid/'

    model = meshxx_p4est_include.Model(name,mdpa_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    params = {}
    params['write_output_per_each_step'] = output
    params['compute_condition_number'] = False
    params['analytical_solution'] = analytical_solution
    params["report_convergence_name"] = "report_el.grf"
    params["report_disp_b_name"] = "convergence_el_disp_y_0_b.grf"
    params["stress_recovery_type"] = stress_recovery_type
    params["neighbour_expansion_level"] = neighbour_expansion_level
    #####p4est parameters#####
    params["order"] = 1
    params["element_name"] = "KinematicLinear"
    params["condition_name"] = "LinePressure"
    params["initial_refinement_process"] = p4est_simulator.P4RefinementProcessAll(nrefine)
    ##########################
    sim = simulator.Simulator(params)
    sim.Initialize(model)
    sim.Update(model)
    sim.Run(model, logging=logging)

    return model

def test():
    model = main(logging=False, output=False, nrefine=1, stress_recovery_type=1, neighbour_expansion_level=1)
    l2_error_ref = 1.4467271166261504e-03
    h1_error_ref = 3.5049873600182184e-02
    h1_error_rs_ref = 6.6422792656567146e-03
    print("%.16e" % model.l2_error)
    print("%.16e" % model.h1_error)
    print("%.16e" % model.h1_error_rs)
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    assert(abs(model.h1_error_rs - h1_error_rs_ref) < 1e-10)
    print("Test passed")

def study():
    ifile = open("convergence.log", "w")
    ifile.write("%-*s%-*s%-*s%s\n" % (10, "ndofs", 30, "l2_error", 30, "h1_error", "h1_error_rs"))
    for i in range(0, 5):
        model = main(logging=False, output=False, nrefine=i, stress_recovery_type=1, neighbour_expansion_level=2)
        ifile.write("%-*d%-*.16e%-*.16e%.16e\n" % (10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, 30, model.h1_error, model.h1_error_rs))
        ifile.flush()
    ifile.close()

def tag():
    return "p4est"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True, nrefine=4, stress_recovery_type=1, neighbour_expansion_level=2)

timer = Timer()
print(timer)
end = time_module.time()
print("Analysis completed, time = " + str(end-start) + " s")
