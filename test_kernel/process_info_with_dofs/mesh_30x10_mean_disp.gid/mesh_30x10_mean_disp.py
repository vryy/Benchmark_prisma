##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import mesh_30x10_mean_disp_include
from mesh_30x10_mean_disp_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = mesh_30x10_mean_disp_include.Model('mesh_30x10_mean_disp',current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    # user-defined script is used (will be appended automatically)
    # ============================================ #
    # |       USER CALCULATION SCRIPT            | #
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv #

    # use special process info to handle global multiplier dofs
    model.model_part.ProcessInfo = ProcessInfoWithDofs()

    # # add Lagrange multipliers for global constraint
    # for i in range(0, 4):
    #     model.model_part.ProcessInfo.AddDof(LAGRANGE_MULTIPLIER_CONSTRAINT)

    # loading conditions
    tol = 1e-6
    load = 0.5e-1
    for node in model.model_part.Nodes:
        if (abs(node.Y0 - 0.5) < tol) or (abs(node.Y0 + 0.5) < tol):
            node.SetSolutionStepValue(FACE_LOAD_Y, load)

    time = 0.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(logging=False, output=False)

    strain_energy = model.solver.solver.GetStrainEnergy()
    print("strain energy = %.10e" % (strain_energy))

    ######### pytesting results #########
    ref_strain_energy = 7.4769903640e-02
    assert(abs(strain_energy - ref_strain_energy) / ref_strain_energy < 1e-10)
    #####################################

def tag():
    tags = ""
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
