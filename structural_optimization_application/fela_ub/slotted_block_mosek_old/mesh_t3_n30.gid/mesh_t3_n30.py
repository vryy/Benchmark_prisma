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
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import mesh_t3_n30_include as simulation_include
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
            node.Fix(VELOCITY_X)
        if abs(node.Y0) < tol:
            node.Fix(VELOCITY_Y)

    model_lm = MosekJ2UpperBoundLimitAnalysisModel2D("slotted_block")
    model_lm.Build(model.model_part)
    model_lm.Solve()
    load_factor = model_lm.ComputeLoadFactor(model.model_part)
    if logging:
        print("load_factor: %.16e" % (load_factor))
        print(model_lm)
    if output:
        model_lm.WriteOutput(model.model_part)
        model.WriteOutput(1.0)

    return model, model_lm

def test():
    model, model_lm = main(logging=False, output=False)
    load_factor = model_lm.ComputeLoadFactor(model.model_part)
    print("load_factor: %.16e" % (load_factor))

    load_factor_ref = 1.1705417830453886e+00

    assert(abs(load_factor - load_factor_ref) < 1e-4)
    print("Test passed")

def tag():
    tags = "FELA"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
    else:
        if not KratosStructuralOptimizationApplication.Has("Mosek"):
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
