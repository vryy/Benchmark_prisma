##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import mesh_64x64_include
try:
    from mesh_64x64_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = mesh_64x64_include.Model('mesh_64x64',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)

        if abs(node.X0-48.0) < tol:
            node.SetSolutionStepValue(FACE_LOAD_X, 0.0)
            # node.SetSolutionStepValue(FACE_LOAD_Y, 100.0/16.0)
            node.SetSolutionStepValue(FACE_LOAD_Y, 1.0)

    model_lm = MosekJ2UpperBoundLimitAnalysisModel2D("cook_membrane")
    model_lm.Build(model.model_part)
    model_lm.Solve()
    if logging:
        print(model_lm)
    if output:
        model_lm.WriteOutput(model.model_part)
        model.WriteOutput(1.0)

    return model, model_lm

def test():
    model, model_lm = main(logging=False, output=False)
    load_factor = model_lm.ComputeLoadFactor(model.model_part)
    print("load_factor: %.16e" % (load_factor))

    load_factor_ref = 1.7557311462781215e+00

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
