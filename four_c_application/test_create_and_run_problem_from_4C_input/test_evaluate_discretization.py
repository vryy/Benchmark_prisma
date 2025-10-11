import sys

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FourCApplication import *
from KratosMultiphysics.DistributedBuildersApplication import *

def main():
    # create 4C model
    fourc_problem = FourCProblem(["contact2D_self_saddlepoint.4C.yaml", "xxx"])
    # print(fourc_problem.GetDiscretizationNames())

    fourc_model = fourc_problem.GetModel()
    fourc_model.FillComplete()

    params = TrilinosParameterList()
    params.Set("action", "calc_struct_nlnstiffmass")
    fourc_model.SetZeroState("structure", 0, "displacement")
    fourc_model.SetZeroState("structure", 0, "residual displacement")
    fourc_model.Evaluate(params, "structure")

    # print(fourc_problem)

def test():
    main()
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
