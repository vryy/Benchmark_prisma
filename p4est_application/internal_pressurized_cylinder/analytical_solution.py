import math

#REF:   Souze de Neto, Computational Plasticity, sec 7.5.1
#       R. Hill, The Mathematical Theory of Plasticity, p.106
#       J. Lubliner, Plasticity Theory, sec 4.3.4
#REMARKS:
#   +   below is the solution for Tresca criterion, do not mix it with Von Mises
class Solution:
    def __init__(self, E, nu, oy, a, b):
        self.E = E
        self.nu = nu
        self.oy = oy
        self.a = a
        self.b = b

    def computePO(self):
        return self.oy/math.sqrt(3) * (1.0 - math.pow(self.a/self.b, 2))

    def computePlim(self):
        return 2.0*self.oy/math.sqrt(3) * math.log(self.b/self.a)

## Solution before the tube's under plastic deformation (P < P0)
class ElasticSolution(Solution):
    def __init__(self, E, nu, oy, a, b):
        Solution.__init__(self, E, nu, oy, a, b)

    # REMARK: the displacement solution is only for plane strain
    def get_displacement(self, P, x, y, z = 0.0):
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
        cost = x/r
        sint = y/r
        A = (1.0 + self.nu) * (1.0 - 2.0*self.nu) * P / (self.E * (math.pow(self.b/self.a, 2) - 1.0))
        B = (1.0 + self.nu) * math.pow(self.b, 2) * P / (self.E * (math.pow(self.b/self.a, 2) - 1.0))
        ur = A*r + B/r
        ut = 0.0
        ux = ur*cost - ut*sint
        uy = ur*sint + ut*cost
        return [ux, uy, 0.0]

    def get_stress_3d(self, P, x, y, z = 0.0):
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
        cost = x/r
        sint = y/r
        orr = -P * (math.pow(self.b/r, 2) - 1.0) / (math.pow(self.b/self.a, 2) - 1.0)
        ott = P * (math.pow(self.b/r, 2) + 1.0) / (math.pow(self.b/self.a, 2) - 1.0)
        ort = 0.0
        ozz = 2.0 * self.nu * P / (math.pow(self.b / self.a, 2) - 1.0)
        oxx = orr*math.pow(cost, 2) + ott*math.pow(sint, 2) - 2.0*ort*sint*cost
        oyy = orr*math.pow(sint, 2) + ott*math.pow(cost, 2) + 2.0*ort*sint*cost
        oxy = (orr - ott)*sint*cost + ort*(math.pow(cost, 2) - math.pow(sint, 2))
        return [oxx, oyy, ozz, oxy, 0.0, 0.0]

## Solution after the tube's under plastic deformation (P > P0)
class PlasticSolution(Solution):
    def __init__(self, E, nu, oy, a, b, cond):
        Solution.__init__(self, E, nu, oy, a, b)
        self.c = a
        if cond == "closed end":
            self.alpha = 1.0 - 2.0*nu
        elif cond == "plane strain":
            self.alpha = 0.0
        elif cond == "open end":
            self.alpha = -2.0*nu
        else:
            self.alpha = 0.0 # default is plane strain

    def get_radius_plastic_front(self, P):
        c = self.a
        while True:
            f = P/(2.0*self.oy/math.sqrt(3)) - math.log(c/self.a) - 0.5*(1.0 - math.pow(c/self.b, 2))
            if abs(f) < 1.0e-10:
                break
            df = -1.0/c + c/math.pow(self.b, 2)
            c = c - f/df
        return c

    def set_load(self, P):
        self.c = self.get_radius_plastic_front(P)

    # REMARK: the displacement solution is only for plane strain
    def get_displacement(self, P, x, y, z = 0.0):
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
        cost = x/r
        sint = y/r
        Y = self.oy/math.sqrt(3)
        if r > self.c: # elastic region
            ur = (1.0+self.nu)*Y/self.E*math.pow(self.c/self.b, 2) * ((1.0-2.0*self.nu)*r + math.pow(self.b, 2) / r)
            ut = 0.0
        else: # plastic region
            orr = -Y * (1.0 - math.pow(self.c/self.b, 2) + 2.0*math.log(self.c/r))
            ur = (1.0-2.0*self.nu)*(1.0+self.nu)/self.E * r * orr + 2.0*Y*math.pow(self.c, 2) * (1.0-math.pow(self.nu, 2)) / self.E / r
            ut = 0.0
        ux = ur*cost - ut*sint
        uy = ur*sint + ut*cost
        return [ux, uy, 0.0]

    def get_stress_3d(self, P, x, y, z = 0.0):
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
        cost = x/r
        sint = y/r
        Y = self.oy / math.sqrt(3)
        if r > self.c: # elastic region
            orr = -Y * math.pow(self.c/self.b, 2) * (math.pow(self.b/r, 2) - 1.0)
            ott = Y * math.pow(self.c/self.b, 2) * (math.pow(self.b/r, 2) + 1.0)
            ort = 0.0
            ozz = 2.0*self.nu * Y * math.pow(self.c/self.b, 2) + Y*self.alpha / (math.pow(self.b/self.a, 2) - 1.0) * (1.0 - math.pow(self.c/self.b, 2) + 2.0*math.log(self.c/self.a))
        else: # plastic region
            orr = -Y * (1.0 - math.pow(self.c/self.b, 2) + 2.0*math.log(self.c/r))
            ott = Y * (1.0 + math.pow(self.c/self.b, 2) - 2.0*math.log(self.c/r))
            ort = 0.0
            ozz = 2.0*self.nu*Y*(math.pow(self.c/self.b, 2) - 2.0*math.log(self.c/r)) + Y*self.alpha / (math.pow(self.b/self.a, 2) - 1.0) * (1.0 - math.pow(self.c/self.b, 2) + 2.0*math.log(self.c/self.a))
        oxx = orr*math.pow(cost, 2) + ott*math.pow(sint, 2) - 2.0*ort*sint*cost
        oyy = orr*math.pow(sint, 2) + ott*math.pow(cost, 2) + 2.0*ort*sint*cost
        oxy = (orr - ott)*sint*cost + ort*(math.pow(cost, 2) - math.pow(sint, 2))
        return [oxx, oyy, ozz, oxy, 0.0, 0.0]


