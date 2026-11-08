import math
import matplotlib.pyplot as plt

from KratosMultiphysics import *
from KratosMultiphysics.MortarApplication import *

# material parameters for masonry joint
# *  0: Young modulus of the unit
# *  1: Poisson ratio of the unit
# *  2: Young modulus of the mortar
# *  3: Poisson ratio of the mortar
# *  4: Joint thickness
# *  5: Tensile strength (ft)
# *  6: Compressive strength (fm)
# *  7: Fracture energy (mode I - Gf_I)
# *  8: Fracture energy (mode II - Gf_II)
# *  9: Cohesion
# *  10: Initial friction angle
# *  11: Residual friction angle
# *  12: Dilatancy angle
# *  13: Css
# *  14: kappa_p, hardening/softening paramter of the cap
# *  15: kappa_m, hardening/softening paramter of the cap
hm = 10.0
ft = 0.25
Gf = 0.018
fm = 10.5
Css = 9.0
GfII = 0.125
Eu = 16700.0
nuu = 0.15
Gu = Eu/(2*(1+nuu))
kn = 82.0
ks = 36.0
Em = kn*hm*Eu/(Eu + kn*hm)
Gm = ks*hm*Gu/(Gu + ks*hm)
num = Em/Gm/2 - 1.0
joint_params = Vector(16)
joint_params[0] = Eu
joint_params[1] = nuu
joint_params[2] = Em
joint_params[3] = num
joint_params[4] = hm
joint_params[5] = ft
joint_params[6] = fm
joint_params[7] = Gf
joint_params[8] = GfII
joint_params[9] = 1.4*ft
joint_params[10] = math.atan(0.75)/math.pi*180
joint_params[11] = math.atan(0.75)/math.pi*180
joint_params[12] = 0.0
joint_params[13] = Css
joint_params[14] = 0.09
joint_params[15] = 0.49

cap_law = CapHardSoftLaw(joint_params)
# print(cap_law)
# print(cap_law.GetValue(0.0))

##################

# npoints = 200
# dx = 0.01
# x = 0.0

# x_list = []
# y_list = []
# for i in range(0, npoints):
#     y = cap_law.GetValue(x)

#     x_list.append(x)
#     y_list.append(y)

#     x += dx

##################

x = 0.0
x_list = []
y_list = []

kappa_p = joint_params[14]
kappa_m = joint_params[15]

npoints = 30
dx = kappa_p / npoints

for i in range(0, npoints):
    y = cap_law.GetValue(x)

    x_list.append(x)
    y_list.append(y)

    x += dx

npoints = 60
dx = (kappa_m - kappa_p) / npoints

for i in range(0, npoints):
    y = cap_law.GetValue(x)

    x_list.append(x)
    y_list.append(y)

    x += dx

npoints = 50
dx = (3*kappa_m - kappa_m) / npoints

for i in range(0, npoints):
    y = cap_law.GetValue(x)

    x_list.append(x)
    y_list.append(y)

    x += dx

##################

plt.plot(x_list, y_list, marker='o')
plt.show()
