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
    fourc_model.CreateElement("dis1", "SOLID3", 1, [1, 2, 3, 4, 5, 6, 7, 8])
    fourc_model.FillComplete()


    print(fourc_model)

def test():
    main()
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
