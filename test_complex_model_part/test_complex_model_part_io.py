from KratosMultiphysics import *

def main():

    # Define a Model
    md = Model()

    mp1 = md.CreateModelPart("ModelPart 1")
    mp2 = md.CreateComplexModelPart("ModelPart 2")
    mp3 = md.CreateGComplexModelPart("ModelPart 3")

    model_part_io = ModelPartIO("cube1")
    mp1.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part_io.ReadModelPart(mp1)
    # print(mp1)

    model_part_io2 = ComplexModelPartIO("cube1")
    mp2.AddNodalSolutionStepVariable(COMPLEX_DISPLACEMENT)
    model_part_io2.ReadModelPart(mp2)
    # print(mp2)

    model_part_io3 = GComplexModelPartIO("cube1c")
    mp3.AddNodalSolutionStepVariable(COMPLEX_DISPLACEMENT)
    model_part_io3.ReadModelPart(mp3)
    # print(mp3)


    print(f"####################{mp1.Name}")

    for node in mp1.Nodes:
        node.SetSolutionStepValue(DISPLACEMENT_X, 1.0)

    for node in mp1.Nodes:
        print(node)

    print(f"####################{mp2.Name}")

    for node in mp2.Nodes:
        node.SetSolutionStepValue(COMPLEX_DISPLACEMENT_X, 1.0)
        # node.SetSolutionStepValue(COMPLEX_DISPLACEMENT_Y, ComplexDouble(2.0, 3.0))
        node.SetSolutionStepValue(COMPLEX_DISPLACEMENT_Y, complex(2.0, 4.0)) # this also works, dont know why

    for node in mp2.Nodes:
        # print(node)
        print(node.GetSolutionStepValue(COMPLEX_DISPLACEMENT))

    print(f"####################{mp3.Name}")

    for node in mp3.Nodes:
        node.SetSolutionStepValue(COMPLEX_DISPLACEMENT_Z, complex(4.0, 3.0))

    for node in mp3.Nodes:
        print(node)
        print(node.GetSolutionStepValue(COMPLEX_DISPLACEMENT))

def test():
    main()

    print("Test passed") # Well, if it can go until here, then we are good to go

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
