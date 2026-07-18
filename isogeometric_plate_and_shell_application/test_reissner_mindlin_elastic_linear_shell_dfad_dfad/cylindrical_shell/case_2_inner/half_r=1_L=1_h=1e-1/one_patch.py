import sys
import os

sys.path.append(os.getcwd() + "/..")
import simulation_include
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.IsogeometricApplication import *

E = 2.6
nu = 0.3
L = 1.0
R = 1.0
h = 1.e-1
p = 1.0
plinear_solver = SuperLUSolver()
# drill_stiff = 1.0
drill_stiff = 0.0

def main(logging=True, output=True, order=2, nsampling = [20, 1]):
    model1 = simulation_include.Model(E, nu, R, L, h, p, order=order, nsampling=nsampling, plinear_solver=plinear_solver, drill_stiff=drill_stiff)
    model = model1.Run(output=output)

    return model

def test():
    model = main(logging=False, output=False, order=2, nsampling=[20, 1])

    multipatch_util = MultiPatchUtility()
    [patch_id, xi] = multipatch_util.LocalCoordinates(model.mpatch, [0.0, 0.0, R], [10, 10])
    patch = model.mpatch[patch_id].GetReference()
    cgf = patch.GridFunction(CONTROL_POINT_COORDINATES)
    point = cgf.GetValue(xi)
    dgf = patch.GridFunction(DISPLACEMENT)
    disp = dgf.GetValue(xi)
    ref_disp = 6.4851693450630012e+00
    print("u: %.16e, diff: %.6e" % (disp[2], disp[2] - ref_disp))

    assert(abs(disp[2] - ref_disp) / ref_disp < 1e-12)

    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(logging=True, output=True, order=2, nsampling=[60, 1])
