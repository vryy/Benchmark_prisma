import sys
import os

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

import geometry_factory

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

def save(logging=False, output=False):
    mpatch, dummy = geometry_factory.CreateTwoSlabs(setting=1, configuration = 1, dist = 0.0, L=1.0, w=1.0, h=1.0)
    # mpatch = Refine(mpatch)
    mpatch.Enumerate()
    print(mpatch)

    if logging:
        mpatch_export.Export(mpatch, "two_cubes.m")

    multipatch_util.CheckInterfaces(mpatch, True)

    serializer1 = Serializer("mpatch", SerializerTraceType.SERIALIZER_TRACE_ERROR)
    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        print(patch.Id)
        patch.Save(serializer1, f"patch_{patch.Id}")

def load():
    serializer2 = Serializer("mpatch", SerializerTraceType.SERIALIZER_TRACE_ERROR)
    patch1 = Patch3DSelector.RealPatch
    patch1.Load(serializer2, "patch_1")
    print(patch1)

def test():
    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
