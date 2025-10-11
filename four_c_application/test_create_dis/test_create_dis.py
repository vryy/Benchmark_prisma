from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FourCApplication import *

def main():
    # create 4C model
    fourc_model = FourCModel(mpi.world, 3)

    # add discretizations
    fourc_model.CreateDiscretization("dis1")
    fourc_model.CreateDiscretization("dis2")
    fourc_model.CreateDiscretization("dis3")

    # add nodes to dis1
    fourc_model.CreateNode("dis1", 1, 0.0, 0.0, 0.0)
    fourc_model.CreateNode("dis1", 2, 1.0, 0.0, 0.0)
    fourc_model.CreateNode("dis1", 3, 1.0, 1.0, 0.0)
    fourc_model.CreateNode("dis1", 4, 0.0, 1.0, 0.0)
    fourc_model.CreateNode("dis1", 5, 0.0, 0.0, 1.0)
    fourc_model.CreateNode("dis1", 6, 1.0, 0.0, 1.0)
    fourc_model.CreateNode("dis1", 7, 1.0, 1.0, 1.0)
    fourc_model.CreateNode("dis1", 8, 0.0, 1.0, 1.0)
    fourc_model.CreateElement("dis1", "SOLID", 1, [1, 2, 3, 4, 5, 6, 7, 8])
    fourc_model.FillComplete()

    # fourc_model.EvaluateSystem("dis1")

    print(fourc_model)

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
