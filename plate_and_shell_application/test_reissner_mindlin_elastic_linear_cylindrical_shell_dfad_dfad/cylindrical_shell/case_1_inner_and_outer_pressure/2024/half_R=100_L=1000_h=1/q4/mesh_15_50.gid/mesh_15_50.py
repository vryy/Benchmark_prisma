##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mo 15. Jul 17:15:35 CEST 2024
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import mesh_15_50_include
from mesh_15_50_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

import simulator

def main(logging=True, output=True, p=1e2):
    model = mesh_15_50_include.Model('mesh_15_50',current_dir_,current_dir_,logging)
    model.InitializeModel()

    model, l2_error = simulator.run(model, logging=logging, output=output, p=p, compute_l2_error=True)

    return model, l2_error

def test():
    p = 1e2
    model, l2_error = main(logging=False, output=False, p=p)

    l2_error_ref = 5.3058236242e-05
    assert(abs(l2_error - l2_error_ref) < 1e-10)

def tag():
    return "FSDT"

def print_tag():
    print("Tags: " + tag())

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
