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

##################################################################
###  SIMULATION  #################################################
##################################################################

start = time_module.time()

def main(name="mesh_t3", logging=True, output=True, mesh_size=0.1, time=1.0):

    mdpa_path = os.getcwd() + '/../../design_data/' + name + '.gid/'

    model = meshxx_mmg_include.Model(name,mdpa_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    params = {}
    params['write_output_per_each_step'] = output
    params['compute_condition_number'] = False
    #####mmg parameters#####
    params["order"] = 1
    params["element_name"] = "LinearPoisson"
    params["condition_name"] = "LagrangeTemperatureConstraint"
    params["initial_refinement_process"] = mmg_simulator.MmgRefinementProcessAll(mesh_size)
    ##########################
    sim = simulator.Simulator(params)
    sim.Initialize(model)
    sim.Update(model)
    sim.Run(model, logging=logging, time=time)

    return model

def test():
    model = main(logging=False, output=False, mesh_size=0.1)
    l2_error_ref = 4.1654476207086052e-03
    print("%.16e" % model.l2_error)
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    print("Test passed")

def study():
    ifile = open("convergence.log", "w")
    ifile.write("%-*s%-*s%-*s%s\n" % (20, "h", 10, "ndofs", 30, "l2_error", "h1_error"))
    mesh_size=0.2
    time = 1.0
    for i in range(0, 6):
        model = main(logging=False, output=True, mesh_size=mesh_size, time=time)
        ifile.write("%-*.10e%-*d%-*.16e%.16e\n" % (20, mesh_size, 10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, model.h1_error))
        mesh_size *= 0.5
        time += 1.0
        ifile.flush()
    ifile.close()

def tag():
    return "mmg,thermal"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(name='mesh_t3', output=True, logging=True, mesh_size=0.1)

timer = Timer()
print(timer)
end = time_module.time()
print("Analysis completed, time = " + str(end-start) + " s")
