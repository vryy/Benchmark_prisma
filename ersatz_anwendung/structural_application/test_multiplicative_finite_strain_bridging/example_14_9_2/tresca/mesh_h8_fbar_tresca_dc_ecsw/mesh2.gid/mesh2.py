##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__))
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
try:
    import mesh2_include as simulation_include
    from mesh2_include import *
    import KratosMultiphysics
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'mesh2'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def WriteLog(ifile, disp, nodes):
    reac = 0.0
    for node in nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)*4
    ifile.write("%.10e\t%.15e\n" % (disp, reac))
    ifile.flush()

def main(output=True, logging=True, pod="load", number_of_pod_modes=12, dry_run=False, check=False):

    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging)
    model.InitializeModel()

    pod_utils = POD_Utils()

    if pod == "save":
        pod_process = EcswSnapshotCollectingProcess(model.solver.solver.linear_solver)
        pod_process.SetForceTolerance(1.0)
        pod_process.SetNormalize(False)
        model.solver.solver.builder_and_solver.SetPodProcess(pod_process)
        model.solver.solver.linear_solver = pod_process
    elif pod == "load":
        # POD without hyper reduction
        Phi = pod_utils.ReadMat("Gb.mat", 'Phi')
        pod_process = RayleighRitzProjectionProcess(Phi)
        model.solver.solver.builder_and_solver.SetPodProcess(pod_process)
        decouple_build_and_solve = model.analysis_parameters['decouple_build_and_solve']
        if decouple_build_and_solve:
            model.solver.solver.linear_solver = pod_process
            # this must be set if decouple_build_and_solve == True
            # since pod_process is also a linear solver, this is usable and necessary
    elif pod == "load-ecsw":
        # POD with hyper reduction
        Phi = pod_utils.ReadMat("Gb.mat", 'Phi')
        eid = pod_utils.ReadIntVec("Gb.mat", 'eid')
        ewi = pod_utils.ReadIntVec("Gb.mat", 'ewi')
        pod_utils.ListVariables("ew.mat")
        w = pod_utils.ReadVec("ew.mat", 'x')
        pod_process = RayleighRitzProjectionProcess(Phi)
        model.solver.solver.time_scheme = ElementWeightingScheme(model.solver.solver.time_scheme)
        for i in range(len(eid)):
            model.solver.solver.time_scheme.SetElementWeight(eid[i], w[ewi[i]])
        if check:
            model.solver.solver.time_scheme.SetProjectionOperator(Phi)
        model.solver.solver.builder_and_solver.SetPodProcess(pod_process)
        decouple_build_and_solve = model.analysis_parameters['decouple_build_and_solve']
        if decouple_build_and_solve:
            model.solver.solver.linear_solver = pod_process
            # this must be set if decouple_build_and_solve == True
            # since pod_process is also a linear solver, this is usable and necessary

    ## boundary condition
    ymin = 0.0
    ymax = 2.667000E+01
    xmin = 0.0
    xmax = 6.413
    zmin = 0.0
    tol = 1.0e-6

    prescribed_nodes = []

    for node in model.model_part.Nodes:
        if abs(node.X0 - xmin) < tol:
            node.Fix(DISPLACEMENT_X)

        if abs(node.Z0 - zmin) < tol:
            node.Fix(DISPLACEMENT_Z)

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

    delta_disp_list = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.125, 0.125, 0.25, 0.5]   # 4
    delta_disp_list.extend([0.25]*4)        # 5
    delta_disp_list.extend([0.125]*7)
    delta_disp_list.extend([0.125/2]*2)     # 6
    delta_disp_list.extend([0.125/2]*13)
    delta_disp_list.extend([0.125/4]*6)     # 7

    print("*********LOADING STARTED**********")

    for du in delta_disp_list:
        disp += du
        print("*********LOAD STEP " + str(disp) + " STARTED")
        for node in prescribed_nodes:
            node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, du)
            # node.SetSolutionStepValue(DISPLACEMENT_Y, disp)

        delta_time = du
        time = time + delta_time
        if not dry_run:
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
        G = Matrix(1, 1)
        b = Vector(1)
        windex = model.solver.solver.builder_and_solver.GetPodProcess().ConstructSystem(G, b, number_of_pod_modes)
        # print(windex)
        eid = IntegerVector(len(windex))
        ewi = IntegerVector(len(windex))
        cnt = 0
        for k, i in windex.items():
            eid[cnt] = k
            ewi[cnt] = i
            cnt += 1
        # print(f"eid: {eid}")
        # print(f"ewi: {ewi}")
        pod_utils.WriteMat("Gb.mat", "G", G, False)
        pod_utils.WriteVec("Gb.mat", "b", b, True)
        pod_utils.WriteMat("Gb.mat", "Phi", Phi, True)
        pod_utils.WriteIntVec("Gb.mat", "eid", eid, True)
        pod_utils.WriteIntVec("Gb.mat", "ewi", ewi, True)

        w = Vector(1)
        LapackNnlsSolver.Solve(G, w, b, 'QR', 100, -1.0)
        pod_utils.WriteVec("ew.mat", "x", w, False)

    print("Analysis completed")

    return model

def test1():
    model = main(output=False, logging=False, pod = "save")

    ######### pytesting results #########
    ref_reac = 2.0759262156662025e+01
    reac = 0.0
    for node in model.prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)*4
    print("reac: %.16e" % (reac))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    #####################################
    print("Test 1 passed")
    # 113.52138185501099 s

def test2():
    model = main(output=False, logging=False, pod = "load", dry_run=False, check=False)

    ######### pytesting results #########
    ref_reac = 2.1281790096153781e+01
    reac = 0.0
    for node in model.prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)*4
    print("reac: %.16e" % (reac))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    #####################################
    print("Test 2 passed")
    #  46.44733119010925 s

def test3():
    model = main(output=False, logging=False, pod = "load-ecsw", dry_run=False, check=False)

    ######### pytesting results #########
    ref_reac = 2.1303193960156324e+01
    reac = 0.0
    for node in model.prescribed_nodes:
        reac += node.GetSolutionStepValue(REACTION_Y)*4
    print("reac: %.16e" % (reac))
    assert(abs(reac - ref_reac) / abs(ref_reac) < 1e-10)
    #####################################
    print("Test 3 passed")
    # 43.60981249809265 s

def test():
    test1()
    test2()
    test3()
    print("All local tests passed")

def tag():
    tags = "pod,ecsw,plastic"
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
        main(logging=True, output=True)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
