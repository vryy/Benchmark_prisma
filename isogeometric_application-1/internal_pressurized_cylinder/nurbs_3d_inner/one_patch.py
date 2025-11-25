import sys
import os

import simulation_include
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

def main(output=True, logging=True, order = 2, nsampling = [9, 9]):
    # parameters taken from GeoPDE's ex_plane_strain_ring example
    E = 1.0
    nu = 0.0
    P = 1.0
    r1 = 1.0
    r2 = 2.0
    l = 1.0

    plinear_solver = SuperLUSolver()
    model = simulation_include.Model(E, nu, P, r1, r2, l, order, nsampling, plinear_solver)
    return model.Run(logging=logging, output=output)

def test():
    l2_error, h1_error = main(logging=False, output=False)

    l2_error_ref = 1.2790494188949550e-05
    h1_error_ref = 1.1579479148846048e-03

    print("l2_error: %.16e" % (l2_error))
    print("h1_error: %.16e" % (h1_error))
    assert(abs(l2_error - l2_error_ref) < 1e-10)
    assert(abs(h1_error - h1_error_ref) < 1e-10)

    print("Test passed")

def tag():
    return "IGA"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
