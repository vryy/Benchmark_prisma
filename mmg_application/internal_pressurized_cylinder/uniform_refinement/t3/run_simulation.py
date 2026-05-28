##################################################################
## Convergence study of the internal pressurized cylinder problem
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
include_path = os.getcwd() + '/../../'
sys.path.append(include_path)
import meshxx_mmg_include
from meshxx_mmg_include import *

import simulator
import mmg_simulator
import analytical_solution

##################################################################
###  SIMULATION  #################################################
##################################################################

start = time_module.time()

def main(name="mesh10x10", logging=True, output=True,mesh_size=10.,stress_recovery_type=0,neighbour_expansion_level=1):

    mdpa_path = os.getcwd() + '/../../design_data/' + name + '.gid/'

    model = meshxx_mmg_include.Model(name,mdpa_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    params = {}
    params['write_output_per_each_step'] = output
    params['compute_condition_number'] = False
    params['analytical_solution'] = analytical_solution
    params["report_convergence_name"] = "report_el.grf"
    params["report_disp_b_name"] = "convergence_el_disp_y_0_b.grf"
    params["stress_recovery_type"] = stress_recovery_type
    params["neighbour_expansion_level"] = neighbour_expansion_level
    #####mmg parameters#####
    params["order"] = 1
    params["element_name"] = "KinematicLinear"
    params["condition_name"] = "LinePressure"
    params["initial_refinement_process"] = mmg_simulator.MmgRefinementProcessAll(mesh_size)
    ##########################
    sim = simulator.Simulator(params)
    sim.Initialize(model)
    sim.Update(model)
    sim.Run(model, logging=logging)

    return model

def test():
    model = main(logging=False, output=False, mesh_size=10., stress_recovery_type=1, neighbour_expansion_level=1)
    l2_error_ref = 5.0427222651758773e-03
    h1_error_ref = 7.2412437875972596e-02
    h1_error_rs_ref = 1.7662068635675300e-02
    print("%.16e" % model.l2_error)
    print("%.16e" % model.h1_error)
    print("%.16e" % model.h1_error_rs)
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    assert(abs(model.h1_error_rs - h1_error_rs_ref) < 1e-10)
    print("Test passed")

def study():
    ifile = open("convergence.log", "w")
    ifile.write("%-*s%-*s%-*s%-*s%s\n" % (20, "h", 10, "ndofs", 30, "l2_error", 30, "h1_error", "h1_error_rs"))
    mesh_size=20.0
    for i in range(0, 6):
        model = main(logging=False, output=False, mesh_size=mesh_size, stress_recovery_type=0, neighbour_expansion_level=2)
        ifile.write("%-*.10e%-*d%-*.16e%-*.16e%.16e\n" % (20, mesh_size, 10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, 30, model.h1_error, model.h1_error_rs))
        mesh_size *= 0.5
        ifile.flush()
    ifile.close()

def tag():
    return "mmg"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(name='mesh10x10', output=True, logging=True, mesh_size=2., stress_recovery_type=0, neighbour_expansion_level=2)

timer = Timer()
print(timer)
end = time_module.time()
print("Analysis completed, time = " + str(end-start) + " s")
