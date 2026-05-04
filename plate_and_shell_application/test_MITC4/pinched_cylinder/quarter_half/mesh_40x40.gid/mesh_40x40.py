##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import mesh_40x40_include as simulation_include
try:
    from mesh_40x40_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'mesh_40x40'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    ## boundary condition
    tol = 1e-06
    for node in model.model_part.Nodes:
        if (abs(node.Y0) < tol):
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Z)

        if (abs(node.Z0) < tol):
            node.Fix(DISPLACEMENT_Z)
            node.Fix(ROTATION_X)
            node.Fix(ROTATION_Y)

        if (abs(node.X0 - 300.0) < tol):
            node.Fix(DISPLACEMENT_X)
            node.Fix(ROTATION_Y)
            node.Fix(ROTATION_Z)

        if (abs(node.Y0 - 300.0) < tol):
            node.Fix(DISPLACEMENT_Y)
            node.Fix(ROTATION_X)
            node.Fix(ROTATION_Z)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    w_as = -1.826797e-5 # analytical solution for nu = 0.3
    w = model.model_part.Conditions[1].GetNodes()[0].GetSolutionStepValue(DISPLACEMENT_Z)
    print("Vertical displacement at the loaded point: %.16e, ratio = %.6e" % (w, w/w_as))

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(logging=False, output=False)

    w = model.model_part.Conditions[1].GetNodes()[0].GetSolutionStepValue(DISPLACEMENT_Z)
    w_ref = -1.8180165950553103e-05
    assert(abs(w - w_ref) / abs(w_ref) < 1e-10)
    print("Test passed")

def tag():
    tags = "unknown"
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
