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
## This file is generated on Mi 4. Aug 11:38:04 CEST 2021 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./mesh_hyplas_q4.gid')
import mesh_hyplas_q4_include
from mesh_hyplas_q4_include import *
# calculate insitu-stress for geology_virgin.gid
model = mesh_hyplas_q4_include.Model('mesh_hyplas_q4',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

## boundary condition
ymin = 0.0
ymax = 2.666700E+01
xmin = 0.0
tol = 1.0e-6

list_xmin = []
list_ymin = []
list_ymax = []

for node in model.model_part.Nodes:
    if abs(node.X0 - xmin) < tol:
        node.Fix(DISPLACEMENT_X)
        list_xmin.append(node.Id)

    if abs(node.Y0 - ymin) < tol:
        node.Fix(DISPLACEMENT_Y)
        list_ymin.append(node.Id)

    if abs(node.Y0 - ymax) < tol:
        node.Fix(DISPLACEMENT_Y)
        list_ymax.append(node.Id)

print("list_xmin:", list_xmin)
print("list_ymin:", list_ymin)
print("list_ymax:", list_ymax)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
