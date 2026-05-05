##################################################################
# CDPM2, Single element test - uniaxial compression
# Ref: https://petergrassl.com/blog/cdpm2-single-monotonic/
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import tetra_include as simulation_include
from tetra_include import *
model_name_ = 'tetra'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def Record(ifile, time, model_part):
    disp = model_part.Nodes[2].GetSolutionStepValue(DISPLACEMENT_X)
    reac = model_part.Nodes[2].GetSolutionStepValue(REACTION_X)
    strain = model_part.Elements[1].CalculateOnIntegrationPoints(STRAIN, model_part.ProcessInfo)
    stress = model_part.Elements[1].CalculateOnIntegrationPoints(STRESSES, model_part.ProcessInfo)
    # ifile.write("%-*.10e%-*.10e%-*.10e%-*.10e%-*.10e%-*.10e%-*.10e%-*.10e%.10e\n" % (20, time, 20, disp, 20, reac \
    #     , 20, stress[0][0], 20, stress[0][1], 20, stress[0][2], 20, stress[0][3], 20, stress[0][4], stress[0][5] ))
    ifile.write("%-*.10e%-*.10e%-*.10e%-*.10e%.10e\n" % (20, time, 20, disp, 20, reac \
        , 20, strain[0][0], stress[0][0]))
    ifile.flush()

def main(logging=True, output=True, nsteps=5000, delta_time=1e-3, total_disp = -2e-3):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    model.model_part.Nodes[1].Fix(DISPLACEMENT_X)
    model.model_part.Nodes[1].Fix(DISPLACEMENT_Y)
    model.model_part.Nodes[1].Fix(DISPLACEMENT_Z)

    model.model_part.Nodes[2].Fix(DISPLACEMENT_X) # prescribed dof
    model.model_part.Nodes[2].Fix(DISPLACEMENT_Y)
    model.model_part.Nodes[2].Fix(DISPLACEMENT_Z)

    model.model_part.Nodes[3].Fix(DISPLACEMENT_X)
    model.model_part.Nodes[3].Fix(DISPLACEMENT_Z)

    model.model_part.Nodes[4].Fix(DISPLACEMENT_X)
    model.model_part.Nodes[4].Fix(DISPLACEMENT_Y)

    time = 0.0

    if logging:
        ifile = open("reaction.log", "w")
        # ifile.write("%-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s%s\n" % (20, "time", 20, "disp", 20, "reaction", 20, "sigma-xx", 20, "sigma-yy", 20, "sigma-zz", 20, "sigma-xy", 20, "sigma-yz", "sigma-xz"))
        ifile.write("%-*s%-*s%-*s%-*s%s\n" % (20, "time", 20, "disp", 20, "reaction", 20, "epsilon-xx", "sigma-xx"))
        Record(ifile, time, model.model_part)

    for i in range(0, nsteps):
        time += delta_time
        model.model_part.Nodes[2].SetSolutionStepValue(DISPLACEMENT_X, time*total_disp)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            Record(ifile, time, model.model_part)

    if logging:
        ifile.close()

    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
    return model

def test():
    model = main(logging=False, output=False, nsteps=1500)

    reac = model.model_part.Nodes[2].GetSolutionStepValue(REACTION_X)
    ref_reac = -7.4177514712512902e+03
    print("reac: %.16e, diff: %.6e" % (reac, reac - ref_reac))
    assert(abs(reac - ref_reac) < 1e-10)

    print("Test passed")

def tag():
    return "oofem,cdpm2"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, nsteps=5000, delta_time=1e-3)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
