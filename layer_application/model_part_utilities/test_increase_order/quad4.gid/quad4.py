##################################################################
import sys
import os
import math
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import quad4_include
from quad4_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = quad4_include.Model('quad4',current_dir_,current_dir_,logging=logging)
    model.InitializeModel()
    model.WriteOutput(0.5)

    model_part_utils = ModelPartUtilities()
    new_model_part = ModelPart("quad9")

    model_part_utils.CopyAndIncreaseOrder(model.model_part, new_model_part)
    print(new_model_part)

    model.SetModelPart(new_model_part)
    if output:
        model.WriteOutput(1.0)

    return model

def test():
    model = main(logging=False, output=False)

    assert(len(model.model_part.Nodes) == 96)
    assert(abs(model.model_part.Nodes[96].X0 - 1.0) < 1e-16)
    assert(abs(model.model_part.Nodes[96].Y0 - 0.1) < 1e-16)
    print("Test passed")

def tag():
    tags = ""
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
