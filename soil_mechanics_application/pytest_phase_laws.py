import math

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *

##################################################################
def test():
    alpha = 1.611e6
    theta_s = 0.9
    theta_r = 0.1
    beta = 3.96

    Smax = theta_s
    Smin = theta_r
    pr = math.pow(alpha, 1.0/beta)
    s1 = beta
    s2 = 0.7

    ######

    swcc = VanGenuchtenSWCC(s1, s2, pr, Smin, Smax)

    ds = swcc.GetDerivative(0.5)
    nds = swcc.GetNumericalDerivative(0.5, 1e-6, -1e-6)
    assert(abs(ds - nds) < 1e-10)

    d2s = swcc.GetSecondDerivative(0.5)
    nd2s = swcc.GetNumericalSecondDerivative(0.5, 1e-6, -1e-6)
    assert(abs(d2s - nd2s) < 1e-10)

    ######

    wl = MualemRelativePermeabilityWaterLaw(s2, Smin, Smax)

    dp = wl.GetDerivative(0.5)
    ndp = wl.GetNumericalDerivative(0.5, 1e-6, -1e-6)
    assert(abs(dp - ndp) < 1e-10)

    # d2p = wl.GetSecondDerivative(0.5)
    # nd2p = wl.GetNumericalSecondDerivative(0.5, 1e-6, -1e-6)
    # assert(abs(d2p - nd2p) < 1e-10)

    ######

    al = MualemRelativePermeabilityAirLaw(s2, Smin, Smax)

    dp = al.GetDerivative(0.5)
    ndp = al.GetNumericalDerivative(0.5, 1e-6, -1e-6)
    assert(abs(dp - ndp) < 1e-9)

    ######

    print("Test passed")

def tag():
    return "partially-saturated"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        pass
