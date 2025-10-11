##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025-2026 by Hoang-Giang Bui (UoB)            #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Wed May  7 07:08:51 PM BST 2025
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import square1_include as simulation_include
from square1_include import *
model_name_ = 'square1'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

import spectral_utilities

def main(logging=True, output=True, output_time_series=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging)
    model.InitializeModel()

    # ============================================ #
    # |       USER CALCULATION SCRIPT            | #
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv #

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol):
            node.Fix(COMPLEX_PRESSURE)
            node.SetSolutionStepValue(COMPLEX_PRESSURE, complex(1.0, -0.3))
        if (abs(node.X0 - 1.0) < tol) or (abs(node.Y0) < tol) or (abs(node.Y0 - 1.0) < tol):
            node.Fix(COMPLEX_PRESSURE)
            node.SetSolutionStepValue(COMPLEX_PRESSURE, 0.0)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    if output_time_series:
        omega = model.model_part.Properties[1].GetValue(WAVE_NUMBER)
        time_series = spectral_utilities.CreateUniformTimeSeries(0.0, 1e-2, 1000)
        # model.SetOutputPath(name="square_in_td")
        # spectral_utilities.WriteOutput(model, time_series, omega, COMPLEX_PRESSURE)

        if logging:
            spectral_utilities.WriteNodalOutput(model, time_series, omega, COMPLEX_PRESSURE, 1454)

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 0.5) < tol) and (abs(node.Y0 - 0.5) < tol):
            mon_node = node

    pres = mon_node.GetSolutionStepValue(COMPLEX_PRESSURE)
    # print("%.16e" % (pres.real))
    # print("%.16e" % (pres.imag))
    ref_pres = complex(2.3254786506765721e-01, -6.9764359520295999e-02)
    # print(pres - ref_pres)
    assert(abs(pres - ref_pres) < 1e-10)
    print("Test passed")

def tag():
    return "helmholtz"

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
