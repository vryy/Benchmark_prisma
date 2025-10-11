##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Do 3. Aug 15:41:35 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./quad9_3x3.gid')
import quad9_3x3_include
from quad9_3x3_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(output=True, logging=True):
    model = quad9_3x3_include.Model('quad9_3x3',os.getcwd()+"/",os.getcwd()+"/",logging)
    model.InitializeModel()

    tol = 1e-6

    # set the initial density
    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(DENSITY, 8400.)

    # set the initial temperature
    ref_temp = 293.15

    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(TEMPERATURE, ref_temp)
        if abs(node.Y0) < tol:
            node.Fix(TEMPERATURE)

    # set the reference temperature
    values = [ref_temp]*9
    for elem in model.model_part.Elements:
        elem.SetValuesOnIntegrationPoints(REFERENCE_TEMPERATURE, values, model.model_part.ProcessInfo)

    # set the top pressure
    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(AIR_PRESSURE, 1.0e5) # Pa
        if abs(node.Y0 - 1.0) < tol:
            node.Fix(AIR_PRESSURE)
            node.SetSolutionStepValue(AIR_PRESSURE, 8.0e5) # Pa

    # boundary condition
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) or (abs(node.X0 - 1.0) < tol):
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)

    # gravity load
    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0 #-9.81
    gravity[2] = 0.0
    model.model_part.Properties[1].SetValue(GRAVITY, gravity)

    ###

    time = 0.0
    if output:
        model.WriteOutput(time)

    ###

    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    ###

    return model

def test():
    model = main(logging=False, output=False)

    assert(abs(model.model_part.Nodes[16].GetSolutionStepValue(TEMPERATURE) - 293.15) < 1e-10)
    assert(abs(model.model_part.Nodes[16].GetSolutionStepValue(DISPLACEMENT_Y) - 0.00424489795918) < 1e-10)

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
