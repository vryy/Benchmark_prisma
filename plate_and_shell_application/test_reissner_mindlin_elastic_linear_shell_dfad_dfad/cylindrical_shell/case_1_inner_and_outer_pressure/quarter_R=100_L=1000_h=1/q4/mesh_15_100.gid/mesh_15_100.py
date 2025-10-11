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
## This file is generated on So 23. Jun 07:50:25 CEST 2024
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import mesh_15_100_include
from mesh_15_100_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

import simulator

def main(logging=True, output=True, p=1e2, compute_error=False):
    model = mesh_15_100_include.Model('mesh_15_100',current_dir_,current_dir_,logging)
    model.InitializeModel()

    model = simulator.run(model, logging=logging, output=output, p=p, compute_l2_error = compute_error)

    return model

def test():
    p = 1e2
    model, l2_error = main(logging=False, output=False, p=p, compute_error=True)

    l2_error_ref = 4.5502322501e-06
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
