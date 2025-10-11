import sys
import os

sys.path.append(os.getcwd() + "/../")

import simulator
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *

def main(output=True, logging=True, dt=1.0, num_steps=7200):
    params = {}
    params['order'] = 3
    params['nsampling'] = 8 # 4
    params['logging'] = logging
    params['output'] = output
    params['dt'] = dt
    params['num_steps'] = num_steps
    ##
    model = simulator.Run(params)
    return model

def test():
    model = main(logging=False, output=False, dt=1.0, num_steps=100)

    tol = 1.0e-6
    bottom_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Z0 - 0.0) < tol:
            bottom_nodes.append(node)

    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(WATER_FLOW, model.model_part, SuperLUSolver())
    wf = bottom_nodes[0].GetSolutionStepValue(WATER_FLOW)
    wf2 = wf[2]*60.0/0.01
    print("%.10e" % (wf2))
    ref_wf2 = -2.4973557188e-02
    assert(abs(wf2 - ref_wf2) < 1e-6)

    print("Test passed")

def tag():
    return "partially-saturated,isogeometric"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=False, logging=True)
