import sys
import os

sys.path.append(os.getcwd() + "/..")
import simulation_include
import KratosMultiphysics
from KratosMultiphysics.MKLSolversApplication import *

sys.path.append(os.environ['BENCHMARK_PRISMA'] + "/plate_and_shell_application/test_reissner_mindlin_elastic_linear_shell_dfad_dfad/cylindrical_shell")
import analytical_solution

E = 2.6
nu = 0.3
L = 1.0
R = 1.0
h = 1.0
p = 1.0
plinear_solver = MKLPardisoSolver()
drill_stiff = 0.0 # 1e1 #1e2
ref_result_file = os.getcwd() + os.sep + "case3_sol_r=1_p=1_n=1e4.txt"

def main(output=True, logging=True, order=2, nsampling=[20, 0]):
    model1 = simulation_include.Model(E, nu, R, L, h, p, order=order, nsampling=nsampling, plinear_solver=plinear_solver, drill_stiff=drill_stiff)
    model = model1.Run(output=output, logging=logging)

    ana_sol = analytical_solution.CylindricalShellSolution(ref_result_file)
    l2_error = simulation_include.ComputeL2Error(model.model_part, ana_sol)
    model.l2_error = l2_error
    # print("Mesh size: %.16e" % (model_sim.h))
    print("Global displacement (L2) error: %.16e" % (l2_error))

    return model

def test():
    model = main(logging=False, output=False)

    l2_error_ref = 5.2332692399245011e-04
    assert(abs(model.l2_error - l2_error_ref) < 1e-10)

    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=True, logging=True, order=2, nsampling=[20, 0])
        # main(output=True, logging=True, order=2, nsampling=[80, 0])
        # main(output=True, logging=True, order=2, nsampling=[80, 5])
        # main(output=True, logging=True, order=2, nsampling=[20, 1])
        # main(output=True, logging=True, order=2, nsampling=[2, 1])
