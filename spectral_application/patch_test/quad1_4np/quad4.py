import quad4_include
from quad4_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

def main(output=True, logging=True):
    model = quad4_include.Model('quad4',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    if logging:
        # export mesh to CBGeo
        import cbgeo_mpm_writer
        writer = cbgeo_mpm_writer.CbGeoMPMWriter(2)
        writer.WriteMesh(model.model_part, "mesh.txt")

    # ### DEBUGGING: observe the stiffness and mass matrix
    # mp_utils = ModelPartUtilities()
    # mp_utils.CalculateLocalSystem(model.model_part.Elements[1], model.model_part.ProcessInfo, 1)
    # mp_utils.CalculateMassMatrix(model.model_part.Elements[1], model.model_part.ProcessInfo, 1)
    # # print("----------------")
    # # mp_utils.CalculateLocalSystem(model.model_part.Elements[2], model.model_part.ProcessInfo, 1)
    # # mp_utils.CalculateMassMatrix(model.model_part.Elements[2], model.model_part.ProcessInfo, 1)
    # # print("----------------")
    # # mp_utils.CalculateLocalSystem(model.model_part.Elements[3], model.model_part.ProcessInfo, 1)
    # # mp_utils.CalculateMassMatrix(model.model_part.Elements[3], model.model_part.ProcessInfo, 1)
    # sys.exit(0)
    # ########

    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
        node.Fix(DISPLACEMENT_Z)

    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol):
            node.SetSolutionStepValue(FACE_LOAD_X, 1.0e3)
            node.SetSolutionStepValue(FACE_LOAD_Y, 0.0)

    time = 1.0
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    # time = 1.0
    # model.Solve(time, 0, 0, 0, 0)
    # model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT_X))
    #     print(node.GetSolutionStepValue(DISPLACEMENT_Y))

    return model

def test():
    model = main(logging=False, output=False)

    ######### pytesting results #########
    tol = 1.0e-6
    ref_disp_x = 4.3333333333333333e-3
    for node in model.model_part.Nodes:
        if abs(node.X0 - 1.0) < tol:
            disp_x = node.GetSolutionStepValue(DISPLACEMENT_X)
            assert(abs(disp_x - ref_disp_x) < 1e-12)
    print("Test passed")
    #####################################

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
