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

def create_params(output=True, refine_vector=[], hausdorff_distance=0.2):
    params = {}
    params['write_output_per_each_step'] = output
    params['compute_condition_number'] = False
    #####mmg parameters#####
    params["order"] = 1
    params["element_name"] = "LinearPoisson"
    params["condition_name"] = "LagrangeTemperatureConstraint"
    # params["hausdorff_distance"] = hausdorff_distance
    ###############
    return params

def main(name="mesh_t3", logging=True, output=True, time=1.0, refine_vector=[]):

    mdpa_path = os.getcwd() + '/../../design_data/' + name + '.gid/'

    model = meshxx_mmg_include.Model(name,mdpa_path,os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, 0.2)

    params = create_params(output=output, refine_vector=refine_vector)
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

def study(output=True, logging=True, ratio=0.1, nrefine=10, initial_mesh_size=0.2):
    name = "mesh_t3"
    mdpa_path = os.getcwd() + '/../../design_data/' + name + '.gid/'
    model = meshxx_mmg_include.Model(name,mdpa_path,os.getcwd()+"/",logging=False)
    model.InitializeModel()

    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(NODAL_MMG_SCALAR_METRIC, initial_mesh_size)

    params = create_params(output=output, refine_vector=[], hausdorff_distance=initial_mesh_size)

    sim = simulator.Simulator(params)
    sim.Initialize(model)

    ############

    if logging:
        ifile = open("convergence.log", "w")
        ifile.write("%-*s%-*s%s\n" % (10, "ndofs", 30, "l2_error", "h1_error"))

    refine_vector = []
    stress_util = RecoverStressUtility()
    mesh_size = initial_mesh_size
    for step in range(0, nrefine):
        print(f"#####################################################")
        print(f"############SOLVING the refinement step {step+1}############")
        print(f"#####################################################", flush=True)

        refine_process = mmg_simulator.MmgRefinementProcessBasedOnRefineVector(refine_vector)
        refine_process.Execute(model.model_part)
        mesh_size *= 0.5
        sim.Update(model)
        sim.Run(model, time=step+1, logging=False)

        if logging:
            ifile.write("%-*d%-*.16e%.16e\n" % (10, model.solver.solver.builder_and_solver.GetEquationSystemSize(), 30, model.l2_error, model.h1_error))
            ifile.flush()

        # compute the error and mark the element to refine
        if step < nrefine-1:
            print(f"############COMPUTING THE KELLY ERROR ESTIMATOR############")
            ee_vec = []
            e_vec = []
            stress_util = RecoverStressUtility()
            stress_util.ComputeKellyErrorEstimation2D(model.model_part, TEMPERATURE)
            for elem in model.model_part.Elements:
                ee = elem.GetValue(LOCAL_ERROR)
                e_vec.append(elem.Id)
                ee_vec.append(ee)
            threshold = find_value2(ee_vec, ratio=ratio)
            print(f"threshold: {threshold}")
            refine_vector = []
            for i in range(0, len(e_vec)):
                if ee_vec[i] > threshold:
                    elem = model.model_part.Elements[e_vec[i]]
                    for node in elem.GetNodes():
                        refine_vector.append(node.Id)
            # also refine the boundary
            # refine_vector.extend(model.node_groups['inner'])
            # refine_vector.extend(model.node_groups['outer'])
            refine_vector = list(set(refine_vector))
            # print(f"refine_vector: {refine_vector}")
            # print(f"len(refine_vector): {len(refine_vector)}")
            # print(f"e_vec: {e_vec}")
            # print(f"len(e_vec): {len(e_vec)}")

    if logging:
        ifile.close()

    return sim, model

def test():
    sim, model = study(output=False, logging=False, ratio=0.1, nrefine=5, initial_mesh_size=0.2)
    l2_error_ref = 7.1662723071792660e-04
    h1_error_ref = 3.6299983152814685e-02
    print("%.16e" % model.l2_error)
    print("%.16e" % model.h1_error)
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)
    assert(abs(model.h1_error - h1_error_ref) < 1e-10)
    print("Test passed")

def tag():
    return "mmg"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True, refine_vector=[])

timer = Timer()
print(timer)
end = time_module.time()
print("Analysis completed, time = " + str(end-start) + " s")
