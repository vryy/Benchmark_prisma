import sys
import os

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

kernel = Kernel()   #defining kernel

def test1():
    kv = KnotArray1D([0, 0, 0, 0, 1, 1, 1, 1])
    kv1 = kv.CloneAndRefineInTheMiddle()
    kv2 = kv1.CloneAndRefineInTheMiddle()
    assert(len(kv2) == 11)
    assert(kv2[4] == 0.25)

def test2():
    kv = KnotArray1D([0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1])
    kv1 = kv.CloneAndRefineInTheMiddle()
    kv2 = kv1.CloneAndRefineInTheMiddle()
    assert(len(kv2) == 16)
    assert(kv2[4] == 0.125)
    assert(kv2[7] == 0.5)
    assert(kv2[8] == 0.5)

def test():
    test1()
    test2()
    print("Test passed")

def tag():
    return "IGA"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
