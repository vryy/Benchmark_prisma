import math
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
################################################################################
## parameters extracted from Celia paper:
## A General Mass-Conservative Numerical Solution for the Unsaturated Flow Equation, Water Resources Research
alpha = 1.611e6
theta_s = 0.287
theta_r = 0.075
beta = 3.96
Ks = 0.00944
A = 1.175e6
gamma = 4.74
################################################################################

n = 0.5
Smax = theta_s / n
Smin = theta_r / n
density = 2e3 / 1e6         # kg/cm^3
density_water = 1e3 / 1e6   # kg/cm^3
g = 9.81e2                  # cm/s^2
pr = density_water*g*math.pow(alpha, 1.0/beta)
s1 = beta
s2 = 1.0

swcc = VanGenuchtenSWCC(s1, s2, pr, Smin, Smax)
pc = density_water*g*30.0
sw = swcc.GetValue(pc)
ds_dpc = swcc.GetDerivative(pc)
print(sw)
print(-n*ds_dpc*density_water*g)
