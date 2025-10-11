import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.EkateAuxiliaryApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.IsogeometricSoilMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.BRepApplication import *

kernel = Kernel()   #defining kernel

import model_iga_include
from model_iga_include import *

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
transfer_util = BezierPostUtility()

import geometry_factory

mpatch_import = MultiNURBSPatchGeoImporter2D()

def CreateMultiPatch(order = 2):

    start_point = [0.0, 0.0, 0.0]
    end_point = [0.2, 0.2, 1.0]
    patch_ptr = geometry_factory.CreateSlab(start_point, end_point)
    patch = patch_ptr.GetReference()
    patch.Id = 1

    ######create multipatch
    mpatch = MultiPatch3D()
    mpatch.AddPatch(patch_ptr)

    if order >= 2:
        multipatch_refine_util.DegreeElevate(mpatch[1], [0, 0, order-1])

    return mpatch

def Refine(mpatch, nsampling):
    print("###############REFINEMENT###############")

    ins_knots = []
    nsampling = nsampling
    for i in range(1, nsampling):
        ins_knots.append(float(i)/nsampling)

    multipatch_refine_util.InsertKnots(mpatch[1], [[], [], ins_knots])

    return mpatch

def CreateModel(mpatch, mpatch_sub1, mpatch_sub2, params):
    mpatch_util = MultiPatchUtility()
    element_name = "PartiallySaturatedSoilsKinematicLinearBezier3D"

    mpatch_mmp = MultiMultiPatchModelPart3D()
    mpatch_mmp.AddMultiPatch(mpatch)
    mpatch_mmp.AddMultiPatch(mpatch_sub1)
    mpatch_mmp.AddMultiPatch(mpatch_sub2)

    mpatch_mmp.BeginModelPart()
    model_part = mpatch_mmp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables(model_part)
    model_part.AddNodalSolutionStepVariable(WATER_FLOW)
    # ## for post-processing
    # model_part.AddNodalSolutionStepVariable(LOCAL_DAMAGE)
    # model_part.AddNodalSolutionStepVariable(LOCAL_PLASTICITY)
    # model_part.AddNodalSolutionStepVariable(DAMAGE_EQUATION_STRONG_FORM)
    # model_part.AddNodalSolutionStepVariable(PLASTICITY_EQUATION_STRONG_FORM)

    mpatch_mmp.CreateNodes()

    gravity = Vector(3)
    gravity[0] = 0.0
    gravity[1] = 0.0
    gravity[2] = -9.81

    #problem data
    prop = model_part.Properties[1]
    prop.SetValue(DENSITY,            0 )
    prop.SetValue(YOUNG_MODULUS,      1.3e+06 )
    prop.SetValue(POISSON_RATIO,          0.4 )
    prop.SetValue(THICKNESS, 1.0 )
    prop.SetValue(G_CONSTANT, 9.81 )
    prop.SetValue(GRAVITY, gravity )
    # prop.SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
    # print("Linear elastic model selected for Isotropic3D, description: Soil")
    prop.SetValue(CONSTITUTIVE_LAW, LinearElastic3D() )
    # dummy values to get pass initialization
    prop.SetValue(FIRST_SATURATION_PARAM, 0.0 )
    prop.SetValue(SECOND_SATURATION_PARAM, 0.0 )
    prop.SetValue(AIR_ENTRY_VALUE, 0.0 )
    prop.SetValue(DENSITY_AIR, 0.0 )
    prop.SetValue(BULK_AIR, 0.0 )
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 1)

    for patch_ptr in mpatch.Patches():
        patch = patch_ptr.GetReference()
        patch.LayerIndex = 1
        #        print(patch)

        sub_patch1_ptr = mpatch_sub1[patch.Id]
        sub_patch1 = sub_patch1_ptr.GetReference()
        sub_patch1.LayerIndex = 1

        sub_patch2_ptr = mpatch_sub2[patch.Id]
        sub_patch2 = sub_patch2_ptr.GetReference()
        sub_patch2.LayerIndex = 1

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mmp.AddElements([patch, sub_patch1, sub_patch2], element_name, last_elem_id+1, prop)

    mpatch_mmp.EndModelPart()
    #    print(mpatch_mmp)

    return mpatch_mmp

def Run(params):
    order = params['order']
    nsampling = params['nsampling']
    output = params['output']
    logging = params['logging']
    ##
    mpatch = CreateMultiPatch(order)
    mpatch_sub1 = CreateMultiPatch(order-1)
    mpatch_sub2 = CreateMultiPatch(order-1)
    mpatch = Refine(mpatch, nsampling)
    mpatch_sub1 = Refine(mpatch_sub1, nsampling)
    mpatch_sub2 = Refine(mpatch_sub2, nsampling)

    mpatch_export.Export(mpatch, "Liakopolous.m")
    # mpatch_export.Export(mpatch_sub, "Liakopolous_sub.m")

    #############VIRGIN MODEL#######################################
    virgin_mpatch_mp = CreateModel(mpatch, mpatch_sub1, mpatch_sub2, params)
    virgin_model_part = virgin_mpatch_mp.GetModelPart()

    sim_params = model_iga_include.QuasiStaticParameters()
    sim_params['decouple_build_and_solve'] = False
    sim_params['solving_scheme'] = 'monolithic'
    sim_params['stop_Newton_Raphson_if_not_converge'] = True
    sim_params['list_dof'] = True
    sim_params['convergence_criteria'] = "multiphase"
    sim_params['linear_solver'] = MKLPardisoSolver()
    sim_params['log_residuum'] = True
    sim_params['calculate_reaction'] = False
    model_virgin = model_iga_include.Model('Liakopolous', os.getcwd()+"/", virgin_model_part, sim_params)
    # for node in model.model_part.Nodes:
    #     node.AddDof(DAMAGE, DAMAGE_REACTION)
    #     node.AddDof(PLASTICITY, PLASTICITY_REACTION)
    model_virgin.InitializeModel()
    # sys.exit(0)

    ## post processing
    params_post = {}
    params_post['name'] = "Liakopolous"
    params_post['division mode'] = "uniform"
    params_post['uniform division number'] = 10
#    params_post['division mode'] = "non-uniform"
#    params_post['division number u'] = 10
#    params_post['division number v'] = 10
#    params_post['division number w'] = 1
    post_local_varibles = []
    params_post['variables list'] = [DISPLACEMENT]
    for var in post_local_varibles:
        params_post['variables list'].append(var)
    dim = 3
    def WriteOutput(time):
        for var in post_local_varibles:
            transfer_util.TransferVariablesToNodes(var, model_virgin.model_part, MKLPardisoSolver())
        virgin_mpatch_mp.SynchronizeBackward(0, DISPLACEMENT)
        # mpatch_mp.SynchronizeBackward(1, DAMAGE)
        # mpatch_mp.SynchronizeBackward(1, PLASTICITY)
        for var in post_local_varibles:
            virgin_mpatch_mp.SynchronizeBackward(0, var)
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    # WriteOutput(0.0)
    # sys.exit(0)

    ## solve the virgin model

    aux_util = SoilsAuxiliaryUtility()

    # material parameters
    for elem in model_virgin.model_part.Elements:
        elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
        elem.SetValue(DENSITY,                  2000.0)
        elem.SetValue(DENSITY_WATER,            1.0e3)
        # elem.SetValue(DENSITY_AIR,              1.295)
        # elem.SetValue(BULK_AIR,                 1.188280000e-05)
        elem.SetValue(POROSITY,                 0.2975)
        elem.SetValue(PERMEABILITY_WATER,       4.4e-6)
        elem.SetValue(PERMEABILITY_AIR,         3.2e-7)
        # elem.SetValue(FIRST_SATURATION_PARAM,   2.5)
        # elem.SetValue(SECOND_SATURATION_PARAM,  0.4)
        # elem.SetValue(AIR_ENTRY_VALUE,          3000.0)
        # aux_util.SetValue(SWCC_LAW, VanGenuchtenSWCC(2.5, 0.4, 3000.0), elem)
        # aux_util.SetValue(SWCC_LAW, FullySaturatedSWCC(), elem)
        aux_util.SetValue(SWCC_LAW, LiakopolousSWCC(), elem)
        # aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, SimpleRelativePermeabilityLaw(), elem)
        aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, LiakopolousRelativePermeabilityWaterLaw(), elem)
        aux_util.SetValue(RELATIVE_PERMEABILITY_AIR_LAW, LiakopolousRelativePermeabilityAirLaw(), elem)
        aux_util.SetValue(GAS_LAW, IdealGasLaw(1.295, 1.188280000e-05), elem)
        elem.SetValue(FIX_POROSITY,             False)

        elem.Initialize(model_virgin.model_part.ProcessInfo)

    # boundary condition
    tol = 1.0e-6
    bottom_nodes = []
    for node in model_virgin.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol or abs(node.Y0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
            bottom_nodes.append(node)

        ## uniform water flow
        node.Fix(WATER_PRESSURE)
        node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)

        ## no air flow
        node.Fix(AIR_PRESSURE)
        node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

    nsteps = 1
    time = 0.0
    delta_time = 100.0
    for step in range(0, nsteps):
        gravity = Vector(3)
        gravity[0] = 0.0
        gravity[1] = 0.0
        gravity[2] = -9.81*(step + 1) / nsteps

        model_virgin.model_part.Properties[1].SetValue(GRAVITY, gravity)

        time = time + delta_time
        model_virgin.SolveModel(time)
        # model_virgin.WriteOutput(time)

    isu = InSituStressUtility()
    isu.SetPreStressFromCurrentStress(model_virgin.model_part, model_virgin.model_part.ProcessInfo)
    for elem in model_virgin.model_part.Elements:
        elem.ResetConstitutiveLaw()

    time = time + delta_time
    model_virgin.SolveModel(time)
    # model_virgin.WriteOutput(time)

    max_disp = 0.0
    for node in model_virgin.model_part.Nodes:
        for direction in range(0,3):
            if( abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction])) > max_disp ):
                max_disp = abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction]))

    print("~~ STEP DONE (INSITU STRESS) --> residual displacement= "+str(max_disp)+"~~")
    # sys.exit(0)

    #############SYSTEM MODEL#######################################
    mpatch_mp = CreateModel(mpatch, mpatch_sub1, mpatch_sub2, params)
    model_part = mpatch_mp.GetModelPart()

    sim_params = model_iga_include.QuasiStaticParameters()
    sim_params['dissipation_radius'] = 0.1
    sim_params['decouple_build_and_solve'] = False
    sim_params['solving_scheme'] = 'monolithic'
    sim_params['stop_Newton_Raphson_if_not_converge'] = True
    sim_params['list_dof'] = True
    sim_params['convergence_criteria'] = "multiphase"
    sim_params['linear_solver'] = MKLPardisoSolver()
    sim_params['log_residuum'] = True
    sim_params['calculate_reaction'] = False
    model = model_iga_include.Model('Liakopolous', os.getcwd()+"/", model_part, sim_params)
    # for node in model.model_part.Nodes:
    #     node.AddDof(DAMAGE, DAMAGE_REACTION)
    #     node.AddDof(PLASTICITY, PLASTICITY_REACTION)
    model.InitializeModel()
    # sys.exit(0)

    ## post processing
    params_post = {}
    params_post['name'] = "Liakopolous"
    params_post['division mode'] = "uniform"
    params_post['uniform division number'] = 10
#    params_post['division mode'] = "non-uniform"
#    params_post['division number u'] = 10
#    params_post['division number v'] = 10
#    params_post['division number w'] = 1
    post_local_varibles = []
    params_post['variables list'] = [DISPLACEMENT]
    for var in post_local_varibles:
        params_post['variables list'].append(var)
    dim = 3
    def WriteOutput(time):
        for var in post_local_varibles:
            transfer_util.TransferVariablesToNodes(var, model.model_part, MKLPardisoSolver())
        mpatch_mp.SynchronizeBackward(0, DISPLACEMENT)
        # mpatch_mp.SynchronizeBackward(1, DAMAGE)
        # mpatch_mp.SynchronizeBackward(1, PLASTICITY)
        for var in post_local_varibles:
            mpatch_mp.SynchronizeBackward(0, var)
        model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

    ## solve the system

    # material parameters
    for elem in model.model_part.Elements:
        elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
        elem.SetValue(DENSITY,                  2000.0)
        elem.SetValue(DENSITY_WATER,            1.0e3)
        # elem.SetValue(DENSITY_AIR,              1.295)
        # elem.SetValue(BULK_AIR,                 1.188280000e-05)
        elem.SetValue(POROSITY,                 0.2975)
        elem.SetValue(PERMEABILITY_WATER,       4.4e-6)
        elem.SetValue(PERMEABILITY_AIR,         3.2e-7)
        # elem.SetValue(FIRST_SATURATION_PARAM,   2.5)
        # elem.SetValue(SECOND_SATURATION_PARAM,  0.4)
        # elem.SetValue(AIR_ENTRY_VALUE,          3000.0)
        # aux_util.SetValue(SWCC_LAW, VanGenuchtenSWCC(2.5, 0.4, 3000.0), elem)
        # aux_util.SetValue(SWCC_LAW, FullySaturatedSWCC(), elem)
        aux_util.SetValue(SWCC_LAW, LiakopolousSWCC(), elem)
        # aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, SimpleRelativePermeabilityLaw(), elem)
        aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, LiakopolousRelativePermeabilityWaterLaw(), elem)
        aux_util.SetValue(RELATIVE_PERMEABILITY_AIR_LAW, LiakopolousRelativePermeabilityAirLaw(), elem)
        aux_util.SetValue(GAS_LAW, IdealGasLaw(1.295, 1.188280000e-05), elem)
        elem.SetValue(FIX_POROSITY,             False)

        elem.Initialize(model.model_part.ProcessInfo)

    # boundary condition
    tol = 1.0e-6
    top_nodes = []
    bottom_nodes = []
    lateral_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            lateral_nodes.append(node)
        if abs(node.Y0 - 0.0) < tol or abs(node.Y0 - 0.2) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            lateral_nodes.append(node)
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
            bottom_nodes.append(node)
        if abs(node.Z0 - 1.0) < tol:
            top_nodes.append(node)

        node.Fix(WATER_PRESSURE)
        node.Fix(AIR_PRESSURE)

    # transfer insitu stress
    vtu = VariableTransferUtility(SuperLUSolver())
    vtu.TransferPrestressIdentically( model_virgin.model_part, model.model_part )

    # prestress = model.model_part.Elements[1].GetValuesOnIntegrationPoints(PRESTRESS, model.model_part.ProcessInfo)
    # print(prestress)
    # sys.exit(0)

    transfer_util = BezierPostUtility()
    def Record(ifile, time):
        transfer_util.TransferVariablesToNodes(WATER_FLOW, model.model_part, SuperLUSolver())
        wf = bottom_nodes[0].GetSolutionStepValue(WATER_FLOW)
        ifile.write("%.10e\t%.10e\n" % (time, wf[2]))
        ifile.flush()

    if logging:
        ifile = open("bottom_water_flow.log", "w")
        ifile.write("time\tflow\n")

    ## reset the material one more time to account for new information
    for element in model.model_part.Elements:
        element.ResetConstitutiveLaw()

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    time = 0.0
    delta_time = 1.0
    # delta_time = 100.0

    # first solve to get equilibrium
    nsteps = 1
    for step in range(0, nsteps):
        time = time + delta_time
        model.SolveModel(time)
        if output:
            WriteOutput(time)
        if logging:
            Record(ifile, time)
    max_disp = 0.0
    for node in model.model_part.Nodes:
        for direction in range(0,3):
            if( abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction])) > max_disp ):
                max_disp = abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction]))

    print(f"~~ STEP DONE (APPLICATION OF INSITU STRESS) --> residual displacement= {max_disp} ~~")
    # sys.exit(0)

    # release the water and air pressure on the model
    for node in model.model_part.Nodes:
        node.Free(WATER_PRESSURE)
        node.Free(AIR_PRESSURE)

    # but fix air pressure on top and bottom and lateral
    for node in top_nodes + bottom_nodes + lateral_nodes:
        node.Fix(AIR_PRESSURE)
        node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

    # but fix water pressure on bottom
    for node in bottom_nodes:
        node.Fix(WATER_PRESSURE)
        node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0

    delta_time = params['dt']
    num_steps = params['num_steps']
    for step in range(0, num_steps):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output and step%100 == 0:
            WriteOutput(time)
        if logging:
            Record(ifile, time)

    return model

    ##################################################################
