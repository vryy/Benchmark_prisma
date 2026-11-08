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

def main(logging=True, output=True, nsteps_shear=1000):
    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    ## setup mortar tying links
    util = MortarTyingUtility2D()
    util.SetDefaultSearchParameters()
    util.InitializeSearchTree(model.model_part, 10)
    util.SetValue(TYING_INTEGRATION_ORDER, 3)
    util.SetEchoLevel(1)
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
    mat_params[12] = math.atan(0.01)/math.pi*180
    mat_params[13] = 9.0
    mat_params[14] = 0.09
    mat_params[15] = 0.49
    model.model_part.Properties[2].SetValue(MATERIAL_PARAMETERS, mat_params)
    model.model_part.Properties[2].SetValue(THICKNESS, 1.0)
    model.model_part.Properties[2].SetValue(USE_CONSISTENT_TANGENT, True)
    HardeningLaw().Assign(CAP_HARD_SOFT_LAW, PiecewiseLinearCapHardSoftLaw(mat_params), model.model_part.Properties[2])

    mortar_links = util.SetupTyingLinkElementBased(model.model_part, 10, "masonry_link_composite", model.model_part.Properties[2])

    left_square_nodes = [6, 2, 7, 8]

    ## boundary conditions
    tol = 1e-6
    prescribed_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if (abs(node.Y0) < tol) and (node.Id in left_square_nodes):
            node.Fix(DISPLACEMENT_Y)
        if abs(node.X0 - 2.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    if logging:
        ifile = open("reaction.log", "w")
        ifile.write("disp\treaction\tun\tus\tsigma\ttau\tkappa_1\tkappa_2\tkappa_3\n")

    time = 0.0
    disp = 0.0
    model.SolveModel(time)
    if logging:
        reac = GetReaction(prescribed_nodes)
        for link in mortar_links:
            stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
            gap = link.Calculate(GAPS, model.model_part.ProcessInfo)
            kappa = link.Calculate(INTERNAL_VARIABLES, model.model_part.ProcessInfo)
            break
        # print(stress)
        # print(gap)
        ifile.write("%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, reac, gap[0], gap[1], stress[0], stress[1], kappa[0], kappa[1], kappa[2]))

    ## loading - tension
    delta_disp = 5e-4
    delta_time = 1.0
    nsteps = 100
    for step in range(0, nsteps):
        time += delta_time
        disp += delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, disp)

        print(f"### Tension step {step} started ###", flush=True)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            reac = GetReaction(prescribed_nodes)
            for link in mortar_links:
                stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
                gap = link.Calculate(GAPS, model.model_part.ProcessInfo)
                kappa = link.Calculate(INTERNAL_VARIABLES, model.model_part.ProcessInfo)
                break
            ifile.write("%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, reac, gap[0], gap[1], stress[0], stress[1], kappa[0], kappa[1], kappa[2]))
            ifile.flush()

    ## loading - shear
    time += 100.0
    disp = 0.0
    delta_disp = 1e-3
    delta_time = 1.0
    nsteps = nsteps_shear
    for step in range(0, nsteps):
        time += delta_time
        disp += delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

        print(f"### Shear step {step} started ###", flush=True)

        model.SolveModel(time)
        if output:
            model.WriteOutput(time)

        if logging:
            reac = GetReaction(prescribed_nodes)
            for link in mortar_links:
                stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
                gap = link.Calculate(GAPS, model.model_part.ProcessInfo)
                kappa = link.Calculate(INTERNAL_VARIABLES, model.model_part.ProcessInfo)
                break
            ifile.write("%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, reac, gap[0], gap[1], stress[0], stress[1], kappa[0], kappa[1], kappa[2]))
            ifile.flush()

    if logging:
        ifile.close()

    model.WriteOutput(time)

    return model

def test():
    model = main(logging=False, output=False)

    link = model.model_part.Conditions[3]
    stress = link.Calculate(STRESSES, model.model_part.ProcessInfo)
    gap = link.Calculate(GAPS, model.model_part.ProcessInfo)

    for cond in model.model_part.Conditions:
        if cond.Id >= 3:
            yield_state = cond.Calculate(YIELD_STATE, model.model_part.ProcessInfo)
            # print(yield_state)
            assert(yield_state == 2)

    tau = stress[1]
    us = gap[1]
    print("us: %.16e, tau: %.16e" % (us, tau))
    ref_us = 9.9982852394905508e-01
    ref_tau = 5.2898332517578617e-01
    assert(abs(us - ref_us) < 1e-10)
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
