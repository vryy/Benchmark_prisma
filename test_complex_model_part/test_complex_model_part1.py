from KratosMultiphysics import *

def main():
    # Define a Model
    md = Model()

    mp1 = md.CreateModelPart("ModelPart")
    print(mp1)

    mp2 = md.CreateComplexModelPart("ModelPart2")
    mp2.CreateSubModelPart("sub1")
    print(mp2)

    mp3 = md.CreateGComplexModelPart("ModelPart3")
    mp3.CreateSubModelPart("sub1")
    print(mp3)

    print(DISPLACEMENT)

def test():
    main()

    print("Test passed") # Well, if it can go until here, then we are good to go

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
