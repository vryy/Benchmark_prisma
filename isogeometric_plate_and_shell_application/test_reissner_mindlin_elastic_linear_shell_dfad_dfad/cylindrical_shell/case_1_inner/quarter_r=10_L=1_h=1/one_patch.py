import sys
import os
import numpy as np

sys.path.append(os.getcwd() + "/..")
import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *

E = 2.6
nu = 0.3
L = 1.0
R = 10.0
h = 1.0
p = 1.0
plinear_solver = MKLPardisoSolver()
drill_stiff = 1e-3 # 0.0 #1e2

class AnalyticalSolution:
    def __init__(self, R, s, p):
        self.R = R
        self.s = s
        self.p = p

    def Interpolate(self, w):
        uref = self.p*(self.R**2)/(2*(1+self.s)) * (1 - 0.5*(1-self.s)/self.R)
        return uref

def convergence(n=5, logging=True):
    order = 2
    ny = 1
    nx = 5
    nsamplings = []
    for i in range(0, n):
        nsamplings.append([nx, ny])
        nx *= 2

    ndofs_list = []
    h_list = []
    l2_error_list = []

    s = nu/(1-nu)
    ana_sol = AnalyticalSolution(R, s, p)

    if logging:
        ifile = open("convergence.txt", "w")
        ifile.write("%-*s%-*s%-*s%-*s%-*s%s\n" % (10, "mesh", 10, "nx", 10, "ny", 10, "ndofs", 20, "l2_error", "h"))

    cnt = 1
    for nsampling in nsamplings:
        model1 = simulation_include.Model(E, nu, R, L, h, p, order=order, nsampling=nsampling, plinear_solver=plinear_solver, drill_stiff=drill_stiff)
        model1.mode = 1
        model = model1.Run(output=False)

        l2_error = simulation_include.ComputeL2Error(model.model_part, ana_sol)

        if logging:
            ifile.write("%-*d%-*d%-*d%-*d%-*e%e\n" % (10, cnt, 10, nsampling[0], 10, nsampling[1], 10, model1.ndofs, 20, l2_error, model1.h))
        ndofs_list.append(model1.ndofs)
        h_list.append(model1.h)
        l2_error_list.append(l2_error)

        cnt += 1

    slope, intercept = np.polyfit(np.log(h_list), np.log(l2_error_list), 1)
    print(f"Convergence rate: {slope}")

    if logging:
        ifile.close()

    return l2_error_list

def test():
    l2_error_list = convergence(n=4, logging = False)
    print(f"l2_error_list: {l2_error_list}")

    assert(abs(l2_error_list[0] - 9.044442956421368e-06) < 1e-7)
    assert(abs(l2_error_list[1] - 1.5520692612366603e-07) < 1e-9)
    assert(abs(l2_error_list[2] - 2.383826829212342e-09) < 1e-10)
    assert(abs(l2_error_list[3] - 3.704076262359755e-11) < 1e-12)

    print("Test passed")

def test1():
    model1 = simulation_include.Model(E, nu, R, L, h, p, order=2, nsampling=[30, 1], plinear_solver=plinear_solver, drill_stiff=drill_stiff)
    model1.mode = 1
    model = model1.Run(output=False)

    s = nu/(1-nu)
    ana_sol = AnalyticalSolution(R, s, p)
    uref = ana_sol.Interpolate(0.0)
    print(f"uref: {uref}")
    l2_error = simulation_include.ComputeL2Error(model.model_part, ana_sol)
    print("l2_error: %.6e" % (l2_error))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        order = 2
        # nsampling = [120, 1]
        nsampling = [240, 1]
        # nsampling = [2, 1]
        model1 = simulation_include.Model(E, nu, R, L, h, p, order=order, nsampling=nsampling, plinear_solver=plinear_solver, drill_stiff=drill_stiff)
        model1.mode = 1
        model1.Run(output=True)
