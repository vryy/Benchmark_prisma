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
R = 10.0
h = 1.0
p = 1.0
plinear_solver = MKLPardisoSolver()
drill_stiff = 1e-3 # 0.0 #1e2
ref_result_file = os.environ['BENCHMARK_PRISMA'] + "/isogeometric_plate_and_shell_application/test_reissner_mindlin_elastic_linear_fast_shell/cylindrical_shell/case_3_inner/case3_sol_r=10-3.txt"

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

    ana_sol = analytical_solution.CylindricalShellSolution(ref_result_file, wscale=R/h)

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

def test():
    l2_error_list = convergence(n=4, logging = False)

    assert(abs(l2_error_list[0] - 1.828290e-02) < 1e-5)
    assert(abs(l2_error_list[1] - 2.200379e-03) < 1e-5)
    assert(abs(l2_error_list[2] - 1.251050e-04) < 1e-5)
    assert(abs(l2_error_list[3] - 7.244230e-06) < 1e-5)

    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        order = 2
        # nsampling = [120, 1]
        nsampling = [240, 1]
        # nsampling = [2, 1]
        model1 = simulation_include.Model(E, nu, R, L, h, p, order=order, nsampling=nsampling, plinear_solver=plinear_solver, drill_stiff=drill_stiff)
        model1.Run(output=True)
