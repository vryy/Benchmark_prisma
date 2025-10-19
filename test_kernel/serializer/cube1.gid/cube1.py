##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Sa 14. Mar 00:15:32 CET 2020
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import cube1_include
from cube1_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = cube1_include.Model('cube1',current_dir_,current_dir_,logging)
    model.InitializeModel()

    return model

def test():
    model = main(logging=False, output=False)

    serializer1 = Serializer("restart_data", SerializerTraceType.SERIALIZER_TRACE_ERROR)
    serializer1.Save("model_part", model.model_part)

    serializer2 = Serializer("restart_data", SerializerTraceType.SERIALIZER_TRACE_ERROR)
    model_part1 = ModelPart()
    serializer2.Load("model_part", model_part1)
    # print(model_part1)

    assert(len(model.model_part.Nodes) == len(model_part1.Nodes))
    assert(len(model.model_part.Elements) == len(model_part1.Elements))

    print("Test passed")

def tag():
    return "restart"

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
