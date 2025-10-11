import sys

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FourCApplication import *

def main():
    # create 4C model
    fourc_problem = FourCProblem(["contact2D_self_saddlepoint.4C.yaml", "xxx"])
    fourc_problem.Run()
    return fourc_problem.GetDiscretizationNames()

def test():
    output = main()
    assert(output[0] == "structure")
    print("Test passed")

def tag():
    return "4C"

def print_tag():
    print("Tags: " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
