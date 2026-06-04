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
import meshxx_p4est_include
from meshxx_p4est_include import *

import simulator
import p4est_simulator

##################################################################
###  SIMULATION  #################################################
##################################################################

start = time_module.time()

def main(name="mesh_q4", logging=True, output=True, nrefine=0, time=1.0):

    mdpa_path = os.getcwd() + '/../../design_data/' + name + '.gid/'

    model = meshxx_p4est_include.Model(name,mdpa_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    params = {}
    params['write_output_per_each_step'] = output
    params['compute_condition_number'] = False
    #####p4est parameters#####
    params["order"] = 1
    params["element_name"] = "LinearPoisson"
    params["initial_refinement_process"] = p4est_simulator.P4RefinementProcessAll(nrefine)
    ##########################
    sim = simulator.Simulator(params)
    sim.Initialize(model)
    sim.Update(model)
    sim.Run(model, logging=logging, time=time)

    return model

def test():
    model = main(logging=False, output=False, nrefine=1)
    l2_error_ref = 8.1849408947130067e-03
    h1_error_ref = 8.8215271059204708e-02
    print("%.16e" % model.l2_error)
    print("%.16e" % model.h1_error)
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    print("Test passed")

def study():
    ifile = open("convergence.log", "w")
    ifile.write("%-*s%-*s%s\n" % (10, "ndofs", 30, "l2_error", "h1_error"))
    for i in range(0, 5):
        model = main(logging=False, output=True, nrefine=i, time=i+1)
        ifile.write("%-*d%-*.16e%.16e\n" % (10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, model.h1_error))
        ifile.flush()
    ifile.close()

def tag():
    return "p4est,thermal"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(name='mesh_q4', output=True, logging=True, nrefine=0)

timer = Timer()
print(timer)
end = time_module.time()
print("Analysis completed, time = " + str(end-start) + " s")
