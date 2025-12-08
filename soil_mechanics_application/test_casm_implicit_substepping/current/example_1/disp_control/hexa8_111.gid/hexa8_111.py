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
##################################################################
##################################################################
import hexa8_111_include
from hexa8_111_include import *
##################################################################
###  UTILITIES  ##################################################
##################################################################
def apply_load_isotropic(model, sides, load):
    if 'x' in sides:
        for c in model.layer_cond_sets['loadx']:
            cond = model.model_part.Conditions[c]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, 0.0)
            cond.SetValue(POSITIVE_FACE_PRESSURE, load)
    if 'y' in sides:
        for c in model.layer_cond_sets['loady']:
            cond = model.model_part.Conditions[c]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, 0.0)
            cond.SetValue(POSITIVE_FACE_PRESSURE, load)
    if 'z' in sides:
        for c in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[c]
            cond.SetValue(NEGATIVE_FACE_PRESSURE, 0.0)
            cond.SetValue(POSITIVE_FACE_PRESSURE, load)

def apply_disp_isotropic(model, sides, disp):
    if 'x' in sides:
        for c in model.layer_cond_sets['loadx']:
            cond = model.model_part.Conditions[c]
            for node in cond.GetNodes():
                node.SetSolutionStepValue(DISPLACEMENT_X, disp)
    if 'y' in sides:
        for c in model.layer_cond_sets['loady']:
            cond = model.model_part.Conditions[c]
            for node in cond.GetNodes():
                node.SetSolutionStepValue(DISPLACEMENT_Y, disp)
    if 'z' in sides:
        for c in model.layer_cond_sets['loadz']:
            cond = model.model_part.Conditions[c]
            for node in cond.GetNodes():
                node.SetSolutionStepValue(DISPLACEMENT_Z, disp)

def apply_disp_xyz(model, dispx, dispy, dispz):
    for c in model.layer_cond_sets['loadx']:
        cond = model.model_part.Conditions[c]
        for node in cond.GetNodes():
            node.SetSolutionStepValue(DISPLACEMENT_X, dispx)
    for c in model.layer_cond_sets['loady']:
        cond = model.model_part.Conditions[c]
        for node in cond.GetNodes():
            node.SetSolutionStepValue(DISPLACEMENT_Y, dispy)
    for c in model.layer_cond_sets['loadz']:
        cond = model.model_part.Conditions[c]
        for node in cond.GetNodes():
            node.SetSolutionStepValue(DISPLACEMENT_Z, dispz)

def get_pq(element, process_info):
    p = element.GetValuesOnIntegrationPoints(PRESSURE_P, process_info)
    q = element.GetValuesOnIntegrationPoints(PRESSURE_Q, process_info)
    pc = element.GetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, process_info)
    return {'p': p, 'q': q, 'pc': pc}

def get_strain(element, process_info):
    ev = element.GetValuesOnIntegrationPoints(EQUIVALENT_VOLUMETRIC_STRAIN, process_info)
    eq = element.GetValuesOnIntegrationPoints(EQUIVALENT_DEVIATORIC_STRAIN, process_info)
    return {'ev': ev, 'eq': eq}

def get_u(node):
    ux = node.GetSolutionStepValue(DISPLACEMENT_X)
    uy = node.GetSolutionStepValue(DISPLACEMENT_Y)
    uz = node.GetSolutionStepValue(DISPLACEMENT_Z)
    return [ux, uy, uz]

def get_e(element, process_info):
    e = element.GetValuesOnIntegrationPoints(VOID_RATIO, process_info)
    return {'e': e}

def main(output=True, logging=True):
    ##################################################################
    ###  INSITU  #####################################################
    ##################################################################
    # calculate insitu-stress for geology_virgin.gid
    model_insitu = hexa8_111_include.Model('hexa8_111',os.getcwd()+"/", logging=False)
    model_insitu.material = "linear_elastic"
    model_insitu.InitializeModel()

    ## boundary condition
    tol = 1.0e-6
    for node in model_insitu.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    time = 0.0
    model_insitu.Solve(time, 0, 0, 0, 0)

    ## apply the hydrostatic load
    load = 200000.0
    initial_load = load
    apply_load_isotropic(model_insitu, "xyz", load)

    ## apply the displacements
    #for node in model_insitu.model_part.Nodes:
    #    if abs(node.X0 - 1.0) < tol:
    #        node.Fix(DISPLACEMENT_X)
    #        node.SetSolutionStepValue(DISPLACEMENT_X, -0.95)

    time = 1.0
    model_insitu.Solve(time, 0, 0, 0, 0)
    # model_insitu.WriteOutput(time)

    # stress = model_insitu.model_part.Elements[1].GetValuesOnIntegrationPoints(STRESSES, model_insitu.model_part.ProcessInfo)
    # print(stress)

    ###################################################################
    ####  SYSTEM  #####################################################
    ###################################################################
    model = hexa8_111_include.Model('hexa8_111',os.getcwd()+"/", logging=logging)
    model.material = "casm"
    model.InitializeModel()

    print("loadx layer:", model.layer_cond_sets['loadx'])
    print("loady layer:", model.layer_cond_sets['loady'])
    print("loadz layer:", model.layer_cond_sets['loadz'])

    ## boundary condition
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        if abs(node.X0 - 1.0) < tol and abs(node.Y0 - 1.0) < tol and abs(node.Z0 - 1.0) < tol:
            tip_node = node

    ## apply the hydrostatic load
    apply_load_isotropic(model, "xyz", load)

    ## set the overconsolidation ratio
    ocr = 1.0
    model.model_part.Properties[1].SetValue(OVERCONSOLIDATION_RATIO, ocr)

    ## setting the prestress and preconsolidation pressure
    for element in model.model_part.Elements:
        stresses = model_insitu.model_part.Elements[element.Id].GetValuesOnIntegrationPoints(STRESSES, model_insitu.model_part.ProcessInfo)
        prestresses = []
        preconsolidation_pressures = []
        for stress in stresses:
            #print("stress:", stress)
            prestress = []
            for s in stress:
                prestress.append(-s)
            prestresses.append(prestress)
            preconsolidation_pressures.append(-ocr*(stress[0] + stress[1] + stress[2]) / 3)
        print("prestresses: " + str(prestresses))
        print("preconsolidation_pressures: " + str(preconsolidation_pressures))
        element.SetValuesOnIntegrationPoints(PRESTRESS, prestresses, 6, model.model_part.ProcessInfo)
        element.SetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, preconsolidation_pressures, model.model_part.ProcessInfo)
        element.ResetConstitutiveLaw()

    ## solve zero load # here the external load and prestress balances
    print("##############################################################")
    print("## SOLVING ZERO LOAD #########################################")
    print("##############################################################")
    time = 2.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)
    u = get_u(tip_node)
    print("displacement at zero: " + str(u))
    # sys.exit(0)

    # ## unload
    # print("##############################################################")
    # print("## UNLOAD ####################################################")
    # print("##############################################################")
    # delta_load = -1.0e4
    # time = 0.0
    # delta_time = 1.0

    # for i in range(0, 10):
    #     time = time + delta_time
    #     load = load + delta_load
    #     apply_load_isotropic(model, "xyz", load)
    #     model.Solve(time, 0, 0, 0, 0)
    #     # model.WriteOutput(time)

    # pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
    # print(pres['p'][0][0])
    # print(pres['q'][0][0])
    u = get_u(tip_node)
    strain = get_strain(model.model_part.Elements[1], model.model_part.ProcessInfo)
    print("displacement at start: " + str(u))
    ev_0 = strain['ev'][0][0]
    print("volumetric strain at start: " + str(ev_0))
    # sys.exit(0)

    print("##############################################################")
    print("## MIXED LOADING #############################################")
    print("##############################################################")

    ## fix the displacement for prescription
    for c in model.layer_cond_sets['loadx']:
        cond = model.model_part.Conditions[c]
        for node in cond.GetNodes():
            node.Fix(DISPLACEMENT_X)
    for c in model.layer_cond_sets['loady']:
        cond = model.model_part.Conditions[c]
        for node in cond.GetNodes():
            node.Fix(DISPLACEMENT_Y)
    for c in model.layer_cond_sets['loadz']:
        cond = model.model_part.Conditions[c]
        for node in cond.GetNodes():
            node.Fix(DISPLACEMENT_Z)

    if logging:
        ifile = open("loading_path.txt", "w")
        ifile.write("p\tq\tpc\tux\tuy\tuz\tev\teq\tea\ter\n")

        pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
        strain = get_strain(model.model_part.Elements[1], model.model_part.ProcessInfo)
        u = get_u(tip_node)
        e = get_e(model.model_part.Elements[1], model.model_part.ProcessInfo)

        ifile.write(str(pres['p'][0][0]))
        ifile.write("\t" + str(pres['q'][0][0]))
        ifile.write("\t" + str(pres['pc'][0][0]))
        ifile.write("\t" + str(u[0]))
        ifile.write("\t" + str(u[1]))
        ifile.write("\t" + str(u[2]))
        ifile.write("\t" + str(strain['ev'][0][0]))
        ifile.write("\t" + str(strain['eq'][0][0]))
        ifile.write("\t" + str(strain['ev'][0][0]/3 + strain['eq'][0][0]))
        ifile.write("\t" + str(strain['ev'][0][0]/3 - strain['eq'][0][0]/2))
        ifile.write("\t" + str(1+e['e'][0][0]))
        ifile.write("\n")
        ifile.flush()

    # print("initial p: " + str(pres['p'][0][0]))
    # print("initial pc: " + str(pres['pc'][0][0]))
    # print("initial e: " + str(e['e'][0][0]))
    # print("initial e_v: " + str(strain['ev'][0][0]))

    # sys.exit(0)

    ######################
    ## LOAD CALCULATION ##
    ######################

    import analytical_solution_disp_control_axial

    _M = model.model_part.Properties[1].GetValue(CSL_SLOPE)
    _lambda = model.model_part.Properties[1].GetValue(VIRGIN_COMPRESSION_INDEX)
    _kappa = model.model_part.Properties[1].GetValue(SWELL_INDEX)
    _nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
    _N = model.model_part.Properties[1].GetValue(SHAPE_PARAMETER)
    _R = model.model_part.Properties[1].GetValue(SPACING_RATIO)
    N = model.model_part.Properties[1].GetValue(VOID_RATIO) # intial void ratio
    ana_sol = analytical_solution_disp_control_axial.proportional_loading_drained_solution(_M, _N, _R, _lambda, _kappa, _nu)

    loadq = 0.0

    p_n = load
    q_n = loadq
    pc_n = ocr*initial_load
    e_q = 0.0
    e_v = 0.0
    k = 3.0

    e = N - _lambda*math.log(pc_n) - _kappa*math.log(p_n/pc_n)
    # print("consistent void ratio: " + str(e))
    v_n = 1+e # initial specific volume

    disp_z_list = []
    disp_x_list = []

    de_a_list = []
    for i in range(0, 400):
        de_a_list.append(1.0e-3)

    p_list = []
    for step in range(0, len(de_a_list)):
        de_a = de_a_list[step]

        # compute the trial elastic solution
        p = ana_sol.get_elas_p(v_n, p_n, q_n, k, de_a)
        q = q_n + k*(p - p_n)
        pc = pc_n

        # check the yield state
        is_yield = ana_sol.check_yield(p, q, pc)

        if is_yield == False:
            de_v = ana_sol.get_elas_de_e_v(v_n, p_n, q_n, k, de_a)
            de_q = ana_sol.get_elas_de_e_q(v_n, p_n, q_n, k, de_a)
        else:
            # compute for yield state
            eta = ana_sol.compute_eta(v_n, p_n, q_n, k, de_a)
            p = ana_sol.get_p(v_n, p_n, q_n, k, eta)
            q = q_n + k*(p - p_n)
            pc = ana_sol.get_pc(v_n, p_n, q_n, k, eta)
            de_v = ana_sol.get_de_v(v_n, p_n, q_n, k, eta)
            de_q = ana_sol.get_de_q(v_n, p_n, q_n, k, eta)

        p_list.append(p)

        e_v = e_v + de_v
        e_q = e_q + de_q
        v = v_n - de_v*v_n # Eq. (16) Peric et al, 2016, here we make it minus because the sample is compressed

        disp_z_list.append(e_v/3 + e_q)
        disp_x_list.append(e_v/3 - e_q/2)

        # update values
        p_n = p
        q_n = q
        v_n = v
        pc_n = pc

    #########

    delta_time = 1

    u = get_u(tip_node)
    initial_disp_x = u[0]
    initial_disp_y = u[1]
    initial_disp_z = u[2]

    delta_time = 1
    e_a = 0.0
    for step in range(0, len(de_a_list)):
        dispx = initial_disp_x - disp_x_list[step]
        dispy = initial_disp_y - disp_x_list[step]
        dispz = initial_disp_z - disp_z_list[step]
        apply_disp_xyz(model, dispx, dispy, dispz)

        time = time + delta_time
        print("###########################")
        print("### Loading step " + str(step) + ", time: " + str(time))
        print("###########################")
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        if logging:
            pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
            strain = get_strain(model.model_part.Elements[1], model.model_part.ProcessInfo)
            u = get_u(tip_node)
            e = get_e(model.model_part.Elements[1], model.model_part.ProcessInfo)

            ifile.write(str(pres['p'][0][0]))
            ifile.write("\t" + str(pres['q'][0][0]))
            ifile.write("\t" + str(pres['pc'][0][0]))
            ifile.write("\t" + str(u[0]))
            ifile.write("\t" + str(u[1]))
            ifile.write("\t" + str(u[2]))
            ifile.write("\t" + str(strain['ev'][0][0]))
            ifile.write("\t" + str(strain['eq'][0][0]))
            ifile.write("\t" + str(strain['ev'][0][0]/3 + strain['eq'][0][0]))
            ifile.write("\t" + str(strain['ev'][0][0]/3 - strain['eq'][0][0]/2))
            ifile.write("\t" + str(1+e['e'][0][0]))
            ifile.write("\n")
            ifile.flush()

            e_v_prescribed = -(dispz + 2.0*dispx)
            e_q_prescribed = -2.0/3*(dispz - dispx)
            e_a_prescribed = e_v_prescribed/3 + e_q_prescribed
            e_a = e_a + de_a_list[step]
            p = p_list[step]
            print("p prescribed: " + str(p))
            print("p computed: " + str(pres['p'][0][0]))
            print("ea diff: " + str(e_a_prescribed - e_a))
            print("ev diff: " + str(e_v_prescribed - strain['ev'][0][0]))
            print("eq diff: " + str(e_q_prescribed - strain['eq'][0][0]))

    if logging:
        ifile.close()

    return model

def test():

    model = main(output=False, logging=False)

    pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
    strain = get_strain(model.model_part.Elements[1], model.model_part.ProcessInfo)
    e = get_e(model.model_part.Elements[1], model.model_part.ProcessInfo)

    ref_p = 3.3290543907566910e+05
    ref_q = 3.9869477437241853e+05
    ref_pc = 6.6307793218443485e+05
    ref_ev = 3.7933119856883968e-02
    ref_eq = 3.8735562671582019e-01
    ref_e = 9.0859880688336847e-01

    mon_p = pres['p'][0][0]
    mon_q = pres['q'][0][0]
    mon_pc = pres['pc'][0][0]
    mon_ev = strain['ev'][0][0]
    mon_eq = strain['eq'][0][0]
    mon_e = e['e'][0][0]

    print("mon_p: %.16e" % mon_p)
    print("mon_q: %.16e" % mon_q)
    print("mon_pc: %.16e" % mon_pc)
    print("mon_ev: %.16e" % mon_ev)
    print("mon_eq: %.16e" % mon_eq)
    print("mon_e: %.16e" % mon_e)

    assert(abs(mon_p - ref_p) / ref_p < 1e-10)
    assert(abs(mon_q - ref_q) / ref_q < 1e-10)
    assert(abs(mon_pc - ref_pc) / ref_pc < 1e-10)
    assert(abs(mon_ev - ref_ev) < 1e-10)
    assert(abs(mon_eq - ref_eq) < 1e-10)
    assert(abs(mon_e - ref_e) < 1e-10)

    print("Test passed")

def tag():
    return "CASM"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)
