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
import mesh_t3_h1_include as simulation_include
try:
    from mesh_t3_h1_include import *
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
            node.Fix(VELOCITY_X)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)

    cohe = 50.0
    phi = 20.0
    model_yield_constraint = MosekDruckerPragerConstraintPlaneStrain(0, cohe, phi)
    model_lm = MosekUpperBoundFelaModel2D("slope_stability", model_yield_constraint)
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

    load_factor_ref = 8.4207939801541212e+00

    assert(abs(load_factor - load_factor_ref) < 1e-10)
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
