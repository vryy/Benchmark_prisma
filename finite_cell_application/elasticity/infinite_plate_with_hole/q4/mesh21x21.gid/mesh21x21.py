##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
import math
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./mesh21x21.gid')
import mesh21x21_include
from mesh21x21_include import *
# calculate insitu-stress for geology_virgin.gid
model = mesh21x21_include.Model('mesh21x21',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
# =====================
# | USER SCRIPT FOR CALCULATION OF EKATE.GID |
# vvvvvvvvvvvvvvvvvvvvv
