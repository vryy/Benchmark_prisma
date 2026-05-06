##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import mesh_80x80_include as simulation_include
try:
    from mesh_80x80_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'mesh_80x80'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    # Abaqus reference dat file: https://docs.software.vt.edu/abaqusv2022/English/SIMACAEBMKRefMap/simabmk-c-pinchcyl.htm#simabmk-c-pinchcyl-t-exafiguresect1__simabmk-c-sxmpincyl-geom

    ## boundary condition
    tol = 1e-06
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) and (abs(node.Y0) < tol) and (abs(node.Z0) < tol):
            node.Fix(DISPLACEMENT_Y)

        if (abs(node.Y0) < tol) or (abs(node.Y0 - 600.0) < tol):
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Z)

        if (abs(node.Z0) < tol):
            node.Fix(DISPLACEMENT_Z)
            node.Fix(ROTATION_X)
            node.Fix(ROTATION_Y)

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
    w_ref = -1.8180166058906933e-05
    assert(abs(w - w_ref) / abs(w_ref) < 1e-10)
    print("Test passed")

def tag():
    tags = "MITC4"
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
