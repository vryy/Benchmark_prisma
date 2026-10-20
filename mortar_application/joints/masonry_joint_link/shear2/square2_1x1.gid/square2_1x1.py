##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025 by Hoang-Giang Bui (UoB)                 #####
#####          2026 by Hoang-Giang Bui (DU)                  #####
##### all rights reserved                                    #####
##################################################################
## This file is generated on (date)
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
import square2_1x1_include as simulation_include
try:
    from square2_1x1_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'square2_1x1'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def GetReaction(nodes):
    reac = 0.0
    for node in nodes:
        reac += node.GetSolutionStepValue(REACTION_X)
    return reac

def main(logging=True, output=True):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    ## setup mortar tying links
    util = MortarTyingUtility2D()
    util.SetDefaultSearchParameters()
    util.InitializeSearchTree(model.model_part, 10)
    util.SetValue(TYING_INTEGRATION_ORDER, 3)
    util.SetEchoLevel(1)
    mortar_links = util.SetupTyingLinkElementBased(model.model_part, 10, "masonry_link_composite", model.model_part.Properties[2])
    # *  0: Young modulus of the unit
    # *  1: Poisson ratio of the unit
    # *  2: Young modulus of the mortar
    # *  3: Poisson ratio of the mortar
    # *  4: Joint thickness
    # *  5: Tensile strength (ft)
    # *  6: Compressive strength (fm)
    # *  7: Fracture energy (mode I - Gf_I)
    # *  8: Fracture energy (mode II - Gf_II)
    # *  9: Cohesion
    # *  10: Initial friction angle
    # *  11: Residual friction angle
    # *  12: Dilatancy angle
    # *  13: Css
    # *  14: kappa_p, hardening/softening paramter of the cap
    # *  15: kappa_m, hardening/softening paramter of the cap
    hm = 10.0
    ft = 0.25
    Eu = 16700.0
    nuu = 0.15
    Gu = Eu/(2*(1+nuu))
    kn = 82.0
    ks = 36.0
    Em = kn*hm*Eu/(Eu + kn*hm)
    Gm = ks*hm*Gu/(Gu + ks*hm)
    num = Em/Gm/2 - 1.0
    mat_params = Vector(16)
    mat_params[0] = Eu
    mat_params[1] = nuu
    mat_params[2] = Em
    mat_params[3] = num
    mat_params[4] = hm
    mat_params[5] = ft
    mat_params[6] = 10.5
    mat_params[7] = 0.018
    mat_params[8] = 0.125
    mat_params[9] = 1.4*ft
    mat_params[10] = math.atan(0.75)/math.pi*180
    mat_params[11] = math.atan(0.75)/math.pi*180
    mat_params[12] = math.atan(0.1)/math.pi*180
    mat_params[13] = 9.0
    mat_params[14] = 0.09
    mat_params[15] = 0.49
    model.model_part.Properties[2].SetValue(MATERIAL_PARAMETERS, mat_params)

    ## boundary conditions
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)
        if node.Id in [3, 1, 4, 5]:
            prescribed_nodes.append(node)

    if logging:
        ifile = open("reaction.log", "w")
        ifile.write("disp\treaction\tun\tus\tsigma\ttau\n")

    time = 0.0
    disp = 0.0
    model.SolveModel(time)
    if logging:
        reac = GetReaction(prescribed_nodes)
        for link in mortar_links:
            stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
            gap = link.Calculate(GAPS, model.model_part.ProcessInfo)
            break
        print(stress)
        print(gap)
        ifile.write("%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, reac, gap[0], gap[1], stress[0], stress[1]))

    ## loading - tension
    delta_disp = 5e-4
    delta_time = 1.0
    nsteps = 100
    for step in range(0, nsteps):
        time += delta_time
        disp += delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, disp)

        print(f"### Compression step {step} started ###", flush=True)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            reac = GetReaction(prescribed_nodes)
            for link in mortar_links:
                stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
                gap = link.Calculate(GAPS, model.model_part.ProcessInfo)
                break
            ifile.write("%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, reac, gap[0], gap[1], stress[0], stress[1]))
            ifile.flush()

    ## loading - shear
    for node in prescribed_nodes:
        node.Fix(DISPLACEMENT_Y)

    time += 100.0
    disp_y = 0.0
    delta_disp_y = 1e-3
    disp_x = disp
    delta_disp_x = 1e-4
    delta_time = 1.0
    nsteps = 500
    for step in range(0, nsteps):
        time += delta_time
        disp_x += delta_disp_x
        disp_y += delta_disp_y
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, disp_x)
            node.SetSolutionStepValue(DISPLACEMENT_Y, disp_y)

        print(f"### Shear step {step} started ###", flush=True)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            reac = GetReaction(prescribed_nodes)
            for link in mortar_links:
                stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
                gap = link.Calculate(GAPS, model.model_part.ProcessInfo)
                break
            ifile.write("%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, reac, gap[0], gap[1], stress[0], stress[1]))
            ifile.flush()

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    link = model.model_part.Conditions[3]
    stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
    gap = link.Calculate(GAPS, model.model_part.ProcessInfo)

    tau = stress[1]
    us = gap[1]
    print("us: %.16e, tau: %.16e" % (us, tau))
    ref_tau = 3.0006394691588299e-02
    assert(abs(tau - ref_tau) < 1e-10)
    print("Test passed")

def tag():
    tags = "joint"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
    return tags

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
