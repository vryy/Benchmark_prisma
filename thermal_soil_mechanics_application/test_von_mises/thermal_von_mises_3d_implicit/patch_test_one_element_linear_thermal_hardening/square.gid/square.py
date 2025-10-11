##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
import math
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./square.gid')
import square_include
from square_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True, output_last=False, total_disp=1e-5*200, delta_disp=1e-5):
    model = square_include.Model('square',os.getcwd()+"/",os.getcwd()+"/", logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            prescribed_nodes.append(node)
        if abs(node.X0 - 0.0) < tol and abs(node.Y0 - 0.0) < tol:
            reaction_node = node
        node.SetSolutionStepValue(TEMPERATURE, 0.0)

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("%-*s%-*s%s\n" % (20, "u", 20, "reaction", "temperature"))

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)
    if logging:
        reaction = reaction_node.GetSolutionStepValue(REACTION_X)
        temperature = reaction_node.GetSolutionStepValue(TEMPERATURE)
        ifile.write("%-*.10e%-*.10e%.10e\n" % (20, 0.0, 20, reaction, temperature))
        ifile.flush()

    disp = 0.0
    delta_time = delta_disp
    nsteps = int(total_disp/delta_disp)
    for i in range(0, nsteps):
        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, disp)
            # node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, delta_disp)

        print("### Time step %f started" % (disp))
        time = time + delta_time
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        if logging:
            reaction = reaction_node.GetSolutionStepValue(REACTION_X)
            temperature = reaction_node.GetSolutionStepValue(TEMPERATURE)
            ifile.write("%-*.10e%-*.10e%.10e\n" % (20, disp, 20, reaction, temperature))
            ifile.flush()

    if logging:
        ifile.close()

    if (not output) and output_last:
        model.WriteOutput(time)

    return model

def test():
    model = main(output=False, logging=False, output_last=False, total_disp=1e-5*200)

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol and abs(node.Y0 - 0.0) < tol:
            reaction_node = node

    ######pytesting######
    rx = reaction_node.GetSolutionStepValue(REACTION_X)
    print("rx: %.15e" % (rx))
    ref_rx = -1.906712506158933e+02
    assert(abs(rx - ref_rx) / abs(ref_rx) < 1e-10)
    t = reaction_node.GetSolutionStepValue(TEMPERATURE)
    print("t: %.15e" % (t))
    ref_t = 1.308128627476452e+02
    assert(abs(t - ref_t) / abs(ref_t) < 1e-10)
    #####################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        main(output=False, logging=True, output_last=True)
