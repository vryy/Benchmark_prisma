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
sys.path.append('./hexa8_111.gid')
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

def apply_load_xyz(model, loadx, loady, loadz):
    for c in model.layer_cond_sets['loadx']:
        cond = model.model_part.Conditions[c]
        cond.SetValue(NEGATIVE_FACE_PRESSURE, 0.0)
        cond.SetValue(POSITIVE_FACE_PRESSURE, loadx)
    for c in model.layer_cond_sets['loady']:
        cond = model.model_part.Conditions[c]
        cond.SetValue(NEGATIVE_FACE_PRESSURE, 0.0)
        cond.SetValue(POSITIVE_FACE_PRESSURE, loady)
    for c in model.layer_cond_sets['loadz']:
        cond = model.model_part.Conditions[c]
        cond.SetValue(NEGATIVE_FACE_PRESSURE, 0.0)
        cond.SetValue(POSITIVE_FACE_PRESSURE, loadz)

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
    if output:
        os.system("rm *.post.bin")

    ##################################################################
    ###  INSITU  #####################################################
    ##################################################################
    # calculate insitu-stress for geology_virgin.gid
    model_insitu = hexa8_111_include.Model('hexa8_111',os.getcwd()+"/",logging=False)
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
    load = 200000
    apply_load_isotropic(model_insitu, "xyz", load)

    ## apply the displacements
    #for node in model_insitu.model_part.Nodes:
    #    if abs(node.X0 - 1.0) < tol:
    #        node.Fix(DISPLACEMENT_X)
    #        node.SetSolutionStepValue(DISPLACEMENT_X, -0.95)

    time = 1.0
    model_insitu.Solve(time, 0, 0, 0, 0)
    if output:
        model_insitu.WriteOutput(time)

    #stress = model_insitu.model_part.Elements[1].GetValuesOnIntegrationPoints(STRESSES, model_insitu.model_part.ProcessInfo)
    #print(stress)

    ###################################################################
    ####  SYSTEM  #####################################################
    ###################################################################
    model = hexa8_111_include.Model('hexa8_111',os.getcwd()+"/",logging=logging)
    model.material = "modified_cam_clay"
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

    ## setting the prestress and preconsolidation pressure
    for element in model.model_part.Elements:
        stresses = model_insitu.model_part.Elements[element.Id].GetValuesOnIntegrationPoints(STRESSES, model_insitu.model_part.ProcessInfo)
        prestresses = []
        preconsolidation_pressures = []
        for stress in stresses:
            prestress = []
            for s in stress:
                prestress.append(-s)
            prestresses.append(prestress)
            preconsolidation_pressures.append(-(stress[0] + stress[1] + stress[2]) / 3)
        print("prestresses: " + str(prestresses))
        print("preconsolidation_pressures: " + str(preconsolidation_pressures))
        element.SetValuesOnIntegrationPoints(PRESTRESS, prestresses, 6, model.model_part.ProcessInfo)
        element.SetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, preconsolidation_pressures, model.model_part.ProcessInfo)
        element.ResetConstitutiveLaw()

    ## solve zero load
    print("##############################################################")
    print("## SOLVING ZERO LOAD #########################################")
    print("##############################################################")
    time = 2.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    if logging:
        pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
        print("p: ", pres['p'][0][0])
        print("q: ", pres['q'][0][0])
        print("pc: ", pres['pc'][0][0])

        pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
        strain = get_strain(model.model_part.Elements[1], model.model_part.ProcessInfo)
        u = get_u(tip_node)
        e = get_e(model.model_part.Elements[1], model.model_part.ProcessInfo)

        ifile = open("loading_path.txt", "w")
        ifile.write("p\tq\tpc\tux\tuy\tuz\tev\teq\tea\ter\tv\n")

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

    print("##############################################################")
    print("## MIXED LOADING #############################################")
    print("##############################################################")
    delta_time = 1
    loadq = 0.0

    delta_load_list = []
    delta_loadq_list = []
    for i in range(0, 26):
        delta_load_list.append(5.0e3)
        delta_loadq_list.append(15.0e3)
    # for i in range(0, 266): #266
    #     delta_load_list.append(5.0e2)
    #     delta_loadq_list.append(15.0e2)
    # for i in range(0, 200):
    #     delta_load_list.append(5.0e2/10)
    #     delta_loadq_list.append(15.0e2/10)
    # for i in range(0, 350):
    #     delta_load_list.append(5.0e2/20)
    #     delta_loadq_list.append(15.0e2/20)
    # for i in range(0, 400):
    #     delta_load_list.append(5.0e2/40)
    #     delta_loadq_list.append(15.0e2/40)

    for step in range(0, len(delta_load_list)):
        load = load + delta_load_list[step]
        loadq = loadq + delta_loadq_list[step]/1.5
        # the coefficient 1.5 is used because (loadz - loady = 1.5*loadq)
        loadx = load - 0.5*loadq
        loady = load - 0.5*loadq
        loadz = load + loadq
        apply_load_xyz(model, loadx, loady, loadz)

        time = time + delta_time
        print("step " + str(step) + ", time: " + str(time))
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

    if logging:
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    ###### pytesting results
    pres = get_pq(model.model_part.Elements[1], model.model_part.ProcessInfo)
    strain = get_strain(model.model_part.Elements[1], model.model_part.ProcessInfo)
    # u = get_u(tip_node)
    e = get_e(model.model_part.Elements[1], model.model_part.ProcessInfo)
    print(pres)
    assert(abs(pres['q'][0][0] - 390000.0) < 1e-6)
    assert(abs(pres['p'][0][0] - 330000.0) < 1e-6)
    assert(abs(pres['pc'][0][0] - 650075.7559122) < 1e-6)
    ########################
    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)
