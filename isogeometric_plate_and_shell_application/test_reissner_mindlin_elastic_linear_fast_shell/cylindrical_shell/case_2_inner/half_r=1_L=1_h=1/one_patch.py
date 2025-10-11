import sys
import os
import numpy as np

sys.path.append(os.getcwd() + "/..")
import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *

sys.path.append(os.environ['HOME'] + "/kratos_bundle/Benchmark_kratos/python3/plate_and_shell_application/test_reissner_mindlin_elastic_linear_shell_dfad_dfad/cylindrical_shell")
import analytical_solution

E = 2.6
nu = 0.3
L = 1.0
R = 1.0
h = 1.0
p = 1.0
plinear_solver = MKLPardisoSolver()
drill_stiff = 1e-3 # set drill stiff to 0 to have optimal convergence

def test():
    l2_error_list = convergence(n=4, logging = False)
    print("l2_error_list:", l2_error_list)

    assert(abs(l2_error_list[0] - 1.706421e-02) < 1e-5)
    assert(abs(l2_error_list[1] - 8.344745e-04) < 1e-5)
    assert(abs(l2_error_list[2] - 7.074857e-05) < 1e-5)
    assert(abs(l2_error_list[3] - 9.685896e-06) < 1e-5)

    print("Test passed")

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

    ana_sol = analytical_solution.CylindricalShellSolution("/home/hbui/workspace/matlab/FSDT/2025/case2_sol_r=1-3.txt", wscale=R/h)

    if logging:
        ifile = open("convergence.log", "w")
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

def tag():
    return "FSDT"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        #
        order = 3
        nsampling = [20, 0]
        ##
        # order = 2
        # nsampling = [5, 1]
        #
        model1 = simulation_include.Model(E, nu, R, L, h, p, order=order, nsampling=nsampling, plinear_solver=plinear_solver, drill_stiff=drill_stiff)
        model1.mode = 1
        model = model1.Run(output=True)

        # ana_sol = analytical_solution.CylindricalShellSolution("/home/buih/workspace/matlab/FSDT/case2_sol_r=1-2.txt")
        # ana_sol = analytical_solution.CylindricalShellSolution("/home/hbui/workspace/matlab/FSDT/2025/case2_sol_r=1.txt")
        # ana_sol = analytical_solution.CylindricalShellSolution("/home/hbui/workspace/matlab/FSDT/2025/case2_sol_r=1-2.txt")
        ana_sol = analytical_solution.CylindricalShellSolution("/home/hbui/workspace/matlab/FSDT/2025/case2_sol_r=1-3.txt")
        l2_error = simulation_include.ComputeL2Error(model.model_part, ana_sol)
        print("Mesh size: %.10e" % (model1.h))
        print("Global displacement (L2) error: %.10e" % (l2_error))
