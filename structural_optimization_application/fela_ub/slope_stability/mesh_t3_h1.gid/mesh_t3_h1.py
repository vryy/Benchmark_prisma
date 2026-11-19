##################################################################
# Slope stability problem, using Drucker-Prager yield constraint
# Ref:
#   +   de Souze Neto, Computational Plasticity.
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_dir_ + "/../")
import simulation_include
try:
    from simulation_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'mesh_t3_h1'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    tol = 1e-06
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 75.0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)

    time = 1.0
    model.SolveModel(time)
    if model.solver.conic_solver.status == "success":
        print("lambda: %.16e" % (model.solver.conic_solver.gamma))
    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))
    if output:
        model.WriteOutput(time)

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(logging=False, output=False)

    lambda_ = model.solver.conic_solver.gamma
    lambda_ref = 7.6875346695739371e+00

    assert(abs(lambda_ - lambda_ref) < 1e-10)
    print("Test passed")

def tag():
    tags = "FELA"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
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
