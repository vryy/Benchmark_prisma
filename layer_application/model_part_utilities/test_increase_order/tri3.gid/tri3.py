##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025 by Hoang-Giang Bui (UoB)                 #####
#####          2026 by Hoang-Giang Bui (DU)                  #####
##### all rights reserved                                    #####
##################################################################
## This file is generated on Sat May 23 20:53:36 BST 2026
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import tri3_include as simulation_include
try:
    from tri3_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'tri3'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    model.WriteOutput(0.0)

    model_part_utils = ModelPartUtilities()
    new_model_part = ModelPart("tri6")

    model_part_utils.CopyAndIncreaseOrder(model.model_part, new_model_part)
    print(new_model_part)

    model.SetModelPart(new_model_part)
    if output:
        model.WriteOutput(1.0)

    return model

def test():
    model = main(logging=False, output=False)

    assert(len(model.model_part.Nodes) == 121)
    assert(abs(model.model_part.Nodes[121].X0 - 1.0) < 1e-16)
    assert(abs(model.model_part.Nodes[121].Y0 - 0.1) < 1e-16)
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
