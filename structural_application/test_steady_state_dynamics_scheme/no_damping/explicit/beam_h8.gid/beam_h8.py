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
## This file is generated on Fr 10. Jul 16:46:15 CEST 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./beam_h8.gid')
import beam_h8_include
from beam_h8_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True):
    model = beam_h8_include.Model('beam_h8',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    rho = 2000.0
    E = 1.91e8
    nu = 0.3
    g = 9.81
    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0
    gravity[2] = -g

    model.model_part.Properties[1].SetValue(DENSITY,         rho )
    model.model_part.Properties[1].SetValue(BODY_FORCE,      ZeroVector(3) )
    model.model_part.Properties[1].SetValue(GRAVITY,         gravity )
    model.model_part.Properties[1].SetValue(YOUNG_MODULUS,      E )
    model.model_part.Properties[1].SetValue(POISSON_RATIO,         nu )
    model.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
    model.InitializeModel()

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    dx = model.model_part.Nodes[62].GetSolutionStepValue(DISPLACEMENT_X)
    dz = model.model_part.Nodes[62].GetSolutionStepValue(DISPLACEMENT_Z)
    # print("%.16e" % (dx))
    # print("%.16e" % (dz))

    ref_dx = -6.0398537421941937e-04
    ref_dz = 5.1396130293870938e-03

    assert(abs(dx - ref_dx) < 1e-10)
    assert(abs(dz - ref_dz) < 1e-10)

    print("Test passed")

def tag():
    return "ssd"

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
