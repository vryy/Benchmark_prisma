##################################################################
import sys
import os
import math
import time as time_module
##################################################################
try:
    import mesh_hyplas_q4_include
    from mesh_hyplas_q4_include import *
    import KratosMultiphysics
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, disp, nodes):
    reac = 0.0
    for node in nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)
    ifile.write("%.10e\t%.15e\n" % (disp, reac))
    ifile.flush()

def main(output=True, logging=True, pod="load", number_of_pod_modes=11):

    model = mesh_hyplas_q4_include.Model('mesh_hyplas_q4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    pod_utils = POD_Utils()

    if pod == "save":
        pod_process = SnapshotCollectingProcess(model.solver.solver.linear_solver)
        model.solver.solver.builder_and_solver.SetPodProcess(pod_process)
        model.solver.solver.linear_solver = pod_process
    elif "load" in pod:
        Phi = pod_utils.ReadMat("Phi.mat", 'Phi')
        pod_process = RayleighRitzProjectionProcess(Phi)
        model.solver.solver.builder_and_solver.SetPodProcess(pod_process)
        decouple_build_and_solve = model.analysis_parameters['decouple_build_and_solve']
        if decouple_build_and_solve:
            model.solver.solver.linear_solver = pod_process
            # this must be set if decouple_build_and_solve == True
            # since pod_process is also a linear solver, this is usable and necessary

    ## boundary condition
    ymin = 0.0
    ymax = 2.666700E+01
    xmin = 0.0
    tol = 1.0e-6

    prescribed_nodes = []

    for node in model.model_part.Nodes:
        if abs(node.X0 - xmin) < tol:
            node.Fix(DISPLACEMENT_X)

        if abs(node.Y0 - ymin) < tol:
            node.Fix(DISPLACEMENT_Y)

        if abs(node.Y0 - ymax) < tol:
            node.Fix(DISPLACEMENT_Y)
            prescribed_nodes.append(node)

    model.prescribed_nodes = prescribed_nodes

    # print("prescribed_nodes:")
    # for node in prescribed_nodes:
    #     print(node.Id)
    # sys.exit(0)

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("disp\treaction\n")

    ## load increment - displacement control

    time = 0.0
    disp = 0.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)
    if logging:
        WriteLog(ifile, disp, prescribed_nodes)

    # delta_disp_list = []
    # for i in range(0, 1):
    #     delta_disp_list.append(0.5)

    delta_disp_list = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.125, 0.125, 0.125, 0.125, 0.25] # adopted from Hyplas 14_9_3_fbar.res
    # delta_disp_list = [0.5, 0.5]
    # delta_disp_list = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.25]
    # delta_disp_list = []
    # for i in range(0, 25):
    #     delta_disp_list.append(0.2)
    # delta_disp_list = []
    # for i in range(0, 5):
    #     delta_disp_list.append(0.02)
    # for i in range(0, 49):
    #     delta_disp_list.append(0.1)

    print("*********LOADING STARTED**********")

    for du in delta_disp_list:
        disp += du
        print("*********LOAD STEP " + str(disp) + " STARTED")
        for node in prescribed_nodes:
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, du)
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

        delta_time = du
        time = time + delta_time
        model.SolveModel(time)
        if output:
            model.WriteOutput(time)
        if logging:
            WriteLog(ifile, disp, prescribed_nodes)

        # print("Displacement")
        # for node in model.model_part.Nodes:
        #     print("%d  %.16e   %.16e" % (node.Id, node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y)))

    if logging:
        ifile.close()

    if pod == "save":
        # S = model.solver.solver.builder_and_solver.GetPodProcess().GetPrincipalValues()
        # print(S)
        # sys.exit(0)
        Phi = model.solver.solver.builder_and_solver.GetPodProcess().GetPrincipalComponents(number_of_pod_modes)
        pod_utils.WriteMat("Phi.mat", "Phi", Phi, False)

    print("Analysis completed")

    return model

def test1():
    model = main(output=False, logging=False, pod = "save")

    ######### pytesting results #########
    ref_reac = 2.778449543107698e+00
    reac = 0.0
    for node in model.prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)
    print("reac: %.15e" % (reac))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    #####################################
    print("Test 1 passed")

def test2():
    model = main(output=False, logging=False, pod = "load")

    ######### pytesting results #########
    ref_reac = 2.779450189964721e+00
    reac = 0.0
    for node in model.prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)
    print("reac: %.15e" % (reac))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    #####################################
    print("Test 2 passed")

def test():
    test1()
    test2()
    print("All local tests passed")

def tag():
    tags = "pod,plastic"
    try:
        test_enabled = all_modules_are_imported_successfully
        test_enabled = test_enabled and KratosErsatzAnwendung.Has('Matio')
        test_enabled = test_enabled and hasattr(KratosMultiphysics, "LapackNnlsSolver")
    except Exception as e:
        test_enabled = False
    if not test_enabled:
        tags += ",untested"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=False, logging=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
