##################################################################
## Convergence study of the internal pressurized cylinder problem
## Since the solution of this problem is smooth, the global convergence
## of AMR does not differ much from global refinement. This example
## is to demonstrate how to use mmg_application in the AMR context.
## Ref: https://scicomp.stackexchange.com/questions/45451/convergence-of-the-internal-pressurized-cylinder-problem-using-adaptive-mesh-ref?noredirect=1#comment93506_45451
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
include_path = os.getcwd() + '/../'
mdpa_path = os.getcwd() + '/../design_data/mesh10x10.gid/'
sys.path.append(include_path)
sys.path.append(mdpa_path)
import meshxx_p4est_include
from meshxx_p4est_include import *

import simulator
import p4est_simulator
import analytical_solution

##################################################################
###  SIMULATION  #################################################
##################################################################

start = time_module.time()

def create_params(output=True, stress_recovery_type=0, neighbour_expansion_level=1, refine_vector=[]):
    params = {}
    params['write_output_per_each_step'] = output
    params['compute_condition_number'] = False
    params['analytical_solution'] = analytical_solution
    params["report_convergence_name"] = "report_el.grf"
    params["report_disp_b_name"] = "convergence_el_disp_y_0_b.grf"
    params["stress_recovery_type"] = stress_recovery_type
    params["neighbour_expansion_level"] = neighbour_expansion_level
    #####p4est parameters#####
    params["order"] = 2
    params["element_name"] = "KinematicLinear"
    params["condition_name"] = "LinePressure"
    params["initial_refinement_process"] = p4est_simulator.P4RefinementProcessBasedOnRefineVector(refine_vector)
    ###############
    return params

def main(logging=True, output=True, time=1.0, refine_vector=[], stress_recovery_type=0, neighbour_expansion_level=1):

    model = meshxx_p4est_include.Model('mesh10x10',mdpa_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    params = create_params(output=output, stress_recovery_type=stress_recovery_type, neighbour_expansion_level=neighbour_expansion_level, refine_vector=refine_vector)
    sim = simulator.Simulator(params)
    sim.Initialize(model)
    sim.Update(model)
    sim.Run(model, time=time, logging=logging)

    return model

# find the threshold in which 40% (ratio) of the element is larger than that
def find_value(vec, ratio=0.4):
    sorted_vec = sorted(vec)
    i = int((1.-ratio)*len(vec))
    return sorted_vec[i]

# find the threshold in which the sum of largest elements exceeding 40% (ratio) of the sum
def find_value2(vec, ratio=0.4):
    sorted_vec = sorted(vec)
    summe = sum(sorted_vec)
    s = 0.0
    for i in range(0, len(sorted_vec)):
        s += sorted_vec[-i-1]
        if s > ratio*summe:
            return sorted_vec[-i-2]

def study(output=True, logging=True, ratio=0.1, nrefine=15):
    model = meshxx_p4est_include.Model('mesh10x10',mdpa_path,os.getcwd()+"/",logging=False)
    model.InitializeModel()

    params = create_params(output=output, stress_recovery_type=1, neighbour_expansion_level=2, refine_vector=[])

    sim = simulator.Simulator(params)
    sim.Initialize(model)

    ############

    if logging:
        ifile = open("convergence.log", "w")
        ifile.write("%-*s%-*s%-*s%s\n" % (10, "ndofs", 30, "l2_error", 30, "h1_error", "h1_error_rs"))
    refine_vector = []
    stress_util = RecoverStressUtility()
    for i in range(0, nrefine):
        print(f"#####################################################")
        print(f"############SOLVING the refinement step {i+1}############")
        print(f"#####################################################", flush=True)

        refine_process = p4est_simulator.P4RefinementProcessBasedOnRefineVector(refine_vector)
        refine_process.Execute(sim.p4est_model)
        sim.Update(model)
        sim.Run(model, time=i+1, logging=False)

        # compute the error and mark the element to refine
        print(f"############COMPUTING THE ZZ ESTIMATOR############")
        ee_vec = []
        e_vec = []
        for elem in model.model_part.Elements:
            ee = stress_util.ComputeZZErrorEstimation(elem, model.model_part.ProcessInfo)
            e_vec.append(elem.Id)
            ee_vec.append(ee)
        threshold = find_value2(ee_vec, ratio=ratio)
        print(f"threshold: {threshold}")
        refine_vector = []
        for i in range(0, len(e_vec)):
            if ee_vec[i] > threshold:
                refine_vector.append(e_vec[i])
        # print(f"refine_vector: {refine_vector}")
        # print(f"len(refine_vector): {len(refine_vector)}")
        # print(f"e_vec: {e_vec}")
        # print(f"len(e_vec): {len(e_vec)}")
        #
        if logging:
            ifile.write("%-*d%-*.16e%-*.16e%.16e\n" % (10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, 30, model.h1_error, model.h1_error_rs))
            ifile.flush()
        #
    if logging:
        ifile.close()

    return sim, model

def study1(output=True, logging=True, ratio=0.1, nrefine=15):
    model = meshxx_p4est_include.Model('mesh10x10',mdpa_path,os.getcwd()+"/",logging=False)
    model.InitializeModel()

    params = create_params(output=output, stress_recovery_type=0, refine_vector=[])

    sim = simulator.Simulator(params)
    sim.Initialize(model)

    ############

    if logging:
        ifile = open("convergence.log", "w")
        ifile.write("%-*s%-*s%-*s%s\n" % (10, "ndofs", 30, "l2_error", 30, "h1_error", "h1_error_rs"))
    refine_vector = []
    stress_util = RecoverStressUtility()
    for i in range(0, nrefine):
        print(f"#####################################################")
        print(f"############SOLVING the refinement step {i+1}############")
        print(f"#####################################################", flush=True)

        refine_process = p4est_simulator.P4RefinementProcessBasedOnRefineVector(refine_vector)
        refine_process.Execute(sim.p4est_model)
        sim.Update(model)
        sim.Run(model, time=i+1, logging=False)

        # compute the error and mark the element to refine
        print(f"############COMPUTING THE KELLY ERROR ESTIMATOR############")
        ee_vec = []
        e_vec = []
        icon = sim.p4est_model.ConstructInterfaces()
        stress_util.ResetLocalError(model.model_part.Elements)
        stress_util.ComputeKellyErrorEstimation(icon, DISPLACEMENT)
        for elem in model.model_part.Elements:
            ee = elem.GetValue(LOCAL_ERROR)
            e_vec.append(elem.Id)
            ee_vec.append(ee)
        threshold = find_value2(ee_vec, ratio=ratio)
        print(f"threshold: {threshold}")
        refine_vector = []
        for i in range(0, len(e_vec)):
            if ee_vec[i] > threshold:
                refine_vector.append(e_vec[i])
        # print(f"refine_vector: {refine_vector}")
        # print(f"len(refine_vector): {len(refine_vector)}")
        # print(f"e_vec: {e_vec}")
        # print(f"len(e_vec): {len(e_vec)}")
        #
        if logging:
            ifile.write("%-*d%-*.16e%-*.16e%.16e\n" % (10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, 30, model.h1_error, model.h1_error_rs))
            ifile.flush()
        #
    if logging:
        ifile.close()

    return sim, model

def test():
    sim, model = study(output=False, logging=False, ratio=0.1, nrefine=8)
    l2_error_ref = 6.9723559666458681e-05
    h1_error_ref = 4.3877765682532495e-03
    h1_error_rs_ref = 6.9957615170103508e-03
    print("%.16e" % model.l2_error)
    print("%.16e" % model.h1_error)
    print("%.16e" % model.h1_error_rs)
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    assert(abs(model.h1_error_rs - h1_error_rs_ref) < 1e-10)
    print("Test passed")

def test1():
    sim, model = study(output=True, logging=False, ratio=0.1, nrefine=2)

    util = P4estUtilities()
    util.DumpHalfEdges(sim.p4est_model)

def tag():
    return "p4est"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True, refine_vector=[], stress_recovery_type=1, neighbour_expansion_level=2)

timer = Timer()
print(timer)
end = time_module.time()
print("Analysis completed, time = " + str(end-start) + " s")
