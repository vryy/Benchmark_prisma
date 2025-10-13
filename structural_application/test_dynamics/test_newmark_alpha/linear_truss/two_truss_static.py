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
## This file is generated on Mi 17. Jun 17:23:39 CEST 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./two_beams.gid')
import two_beams_include
from two_beams_include import *
# calculate insitu-stress for geology_virgin.gid
model = two_beams_include.Model('two_beams',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

B1 = 1.0
B2 = 0.1
H = 0.7
E = 2.0e7
A = 1.0
rho = 1.0
L1 = math.sqrt(H**2 + B1**2)
r = E*A/(3.0*math.sqrt(3.0)) * (H/L1)**3
alpha = 100.0
print("r: " + str(r))

model.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model.model_part.CreateNewNode(2, -B2, H, 0.0)
model.model_part.CreateNewNode(3, B1, H, 0.0)
model.AddDofs(model.model_part)

prop1 = model.model_part.Properties[1]
prop1.SetValue(YOUNG_MODULUS, E)
prop1.SetValue(AREA, A)

prop2 = model.model_part.Properties[2]
prop2.SetValue(YOUNG_MODULUS, alpha*E)
prop2.SetValue(AREA, A)

model.model_part.CreateNewElement("TrussElement3D2N", 1, [1, 3], prop1)
model.model_part.CreateNewElement("TrussElement3D2N", 2, [1, 2], prop2)

model.model_part.CreateNewCondition("PointForce3D", 1, [3], prop2)

print(model.model_part)

model.model_part.Nodes[1].Fix(DISPLACEMENT_Y)
model.model_part.Nodes[2].Fix(DISPLACEMENT_X)
model.model_part.Nodes[2].Fix(DISPLACEMENT_Y)
model.model_part.Nodes[3].Fix(DISPLACEMENT_X)

model.model_part.Nodes[1].Fix(DISPLACEMENT_Z)
model.model_part.Nodes[2].Fix(DISPLACEMENT_Z)
model.model_part.Nodes[3].Fix(DISPLACEMENT_Z)

model.model_part.Elements[1].Initialize(model.model_part.ProcessInfo)
model.model_part.Elements[2].Initialize(model.model_part.ProcessInfo)

fact = -1.0
model.model_part.Nodes[3].SetSolutionStepValue(FORCE_X, 0.0)
model.model_part.Nodes[3].SetSolutionStepValue(FORCE_Y, fact*r)
model.model_part.Nodes[3].SetSolutionStepValue(FORCE_Z, 0.0)

time = 1.0
model.SolveModel(time)
# model.WriteOutput(time)

ux = model.model_part.Nodes[1].GetSolutionStepValue(DISPLACEMENT_X)
uy = model.model_part.Nodes[3].GetSolutionStepValue(DISPLACEMENT_Y)
print("ux: " + str(ux))
print("uy: " + str(uy))

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
