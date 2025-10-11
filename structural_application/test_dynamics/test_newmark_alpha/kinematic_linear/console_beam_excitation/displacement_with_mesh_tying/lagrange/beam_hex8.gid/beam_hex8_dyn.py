##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mi 8. Feb 18:05:36 CET 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./beam_hex8.gid')
import beam_hex8_include
from beam_hex8_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################
def WriteNode(ifile, step, time, exc_node, tip_node):
    dy1 = exc_node.GetSolutionStepValue(DISPLACEMENT_EINS_Y)
    dy2 = tip_node.GetSolutionStepValue(DISPLACEMENT_EINS_Y)
    ifile.write("%-*d%-*.10e%.10e\n" % (20, step, 20, dy1, dy2))
##################################################################

def main(output=True, logging=True):
    model = beam_hex8_include.Model('beam_hex8',os.getcwd()+"/",os.getcwd()+"/", analysis_type=2, logging=logging)
    model.InitializeModel()

    ## setup tying
    tying_util = MortarTyingUtility3D()
    tying_util.SetDefaultSearchParameters()
    tying_util.SetEchoLevel(2)
    tying_util.SetValue(MAXIMAL_DETECTION_DISTANCE, 1.0)
    tying_util.SetValue(GAP_TOLERANCE, 0.51)
    tying_util.SetValue(TYING_INTEGRATION_ORDER, 3)
    tying_util.InitializeSearchTree(model.model_part, 10)
    mortar_links = tying_util.SetupTyingLinkElementBased(model.model_part, 10, "tying_link_geometrical_linear_lagrange_std_mortar")

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)
            prescribed_nodes.append(node)
        if abs(node.X0 - 10.0) < tol and abs(node.Y0 - 1.0) < tol and abs(node.Z0 - 0.0) < tol:
            tip_node = node

    freq = 50.0
    alpha = 0.0
    ampl = 4e-2
    delta_time = 1e-3
    nsteps = 300
    step = 0
    time = 0
    initial_time = time

    if logging:
        exc_node = prescribed_nodes[0]
        ifile = open("displacement.txt", "w")
        ifile.write("%-*s%-*s%s\n" % (20, "step", 20, "dy_exc", "dy_tip"))
        if logging:
            WriteNode(ifile, step, time, exc_node, tip_node)

    for i in range(0, nsteps):
        step += 1
        time += delta_time

        disp = ampl * math.sin(2*math.pi*freq*(time - initial_time) + alpha)
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_NULL_Y, disp)
            node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
            node.SetSolutionStepValue(DISPLACEMENT_EINS_Y, disp)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            WriteNode(ifile, step, time, exc_node, tip_node)

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 10.0) < tol and abs(node.Y0 - 1.0) < tol and abs(node.Z0 - 0.0) < tol:
            tip_node = node

    # print("%.16e" % (tip_node.GetSolutionStepValue(DISPLACEMENT_EINS_Y)))
    ref_disp = 6.6186340279613107e-02
    assert(abs(tip_node.GetSolutionStepValue(DISPLACEMENT_EINS_Y) - ref_disp) < 1e-10)
    print("Test passed")

def tag():
    return "mesh-tying,dynamics,lagrange"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
