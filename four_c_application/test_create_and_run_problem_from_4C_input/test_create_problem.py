import sys

try:
    from KratosMultiphysics import *
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.FourCApplication import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False

def main():
    if not all_modules_are_imported_successfully:
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
    if all_modules_are_imported_successfully:
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
