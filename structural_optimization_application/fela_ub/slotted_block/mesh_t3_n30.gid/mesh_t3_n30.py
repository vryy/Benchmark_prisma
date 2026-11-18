##################################################################
# Slotted block problem:
# Ref:
#   +   Krabbenhoft et al, A general nonlinear optimization algorithm
#   for lower bound limit analysis, IJNME, 2003.
# Remark:
#   +   The mesh is created in GiD with option "Symmetrical Structured"
#   (Utilities/Preferences/Meshing/Structuration Type)
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_dir_ + "/../")
import simulation_include
from simulation_include import *
try:
    from mesh_t3_n30_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'mesh_t3_n30'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    L = 1.0

    tol = 1e-06
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol and (node.Y0 - 0.5*L < tol):
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)

    time = 1.0
    model.SolveModel(time)
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
    lambda_ref = 1.1703392553144596e+00

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
