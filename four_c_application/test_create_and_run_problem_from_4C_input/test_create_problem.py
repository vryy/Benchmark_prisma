import sys

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
try:
    from KratosMultiphysics.FourCApplication import *
    four_c_application_is_available = True
except Exception as e:
    four_c_application_is_available = False

def main():
    if not four_c_application_is_available:
        sys.exit(1)

    # create 4C model
    fourc_problem = FourCProblem(["contact2D_self_saddlepoint.4C.yaml", "xxx"])
    fourc_problem.Run()
    return fourc_problem.GetDiscretizationNames()

def test():
    output = main()
    assert(output[0] == "structure")
    print("Test passed")

def tag():
    if four_c_application_is_available:
        return "4C"
    else:
        return "4C,untested"

def print_tag():
    print("Tags: " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
