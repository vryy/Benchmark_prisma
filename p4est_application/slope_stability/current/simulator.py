#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.P4estApplication import *

def ConstructModelPart(p4est_model):
    # start to construct the Kratos model_part
    p4est_model.BeginModelPart()
    mp = p4est_model.GetModelPart()
    mp.SetBufferSize(2)
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( mp )
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX)

    prop = mp.Properties[1]
    prop.SetValue(CONSTITUTIVE_LAW, DummyConstitutiveLaw())
    prop.SetValue(BODY_FORCE, [0, 0, 0])

    element_name = "KinematicLinear2D9N"

    constraint_name = "LinearMasterSlaveConstraint"

    hanging_nodes = p4est_model.CreateNodes()
    import structural_solver_advanced
    structural_solver_advanced.AddDofs( mp )

    constraint_id = 1
    for hnode in hanging_nodes:
        # print(hnode)
        p4est_model.CreateConstraint(hnode, constraint_id, constraint_name, DISPLACEMENT_X)
        constraint_id = constraint_id + 1
        p4est_model.CreateConstraint(hnode, constraint_id, constraint_name, DISPLACEMENT_Y)
        constraint_id = constraint_id + 1

    layer_idx = 1
    p4est_model.AddElements(layer_idx, element_name, prop)
    p4est_model.EndModelPart()

    return mp

def Run(model, output=True, case_number=7, coarsening=False):
    # set the layer index. It is important to set the layer index.
    # The elements will be created based on layer
    for elem in model.model_part.Elements:
        elem.SetValue(LAYER_INDEX, 1)

    # specify the order of the integration
    p4est_order = P4estSecondOrder()

    # create a sample quad for the p4est tree
    p4est_quad_nodal = P4estQuadData(p4est_order.Value)
    p4est_quad_nodal.Initialize()
    p4est_quad_nodal.Register(DISPLACEMENT)
    p4est_quad_nodal.Register(REACTION)
    p4est_quad_nodal.Finalize()

    p4est_quad_int = P4estQuadData(p4est_order.Value)
    p4est_quad_int.Initialize()
    p4est_quad_int.Register(THREED_STRESSES, 6)
    p4est_quad_int.Register(PLASTICITY_INDICATOR)
    p4est_quad_int.Register(INTEGRATION_POINT_STRAIN_VECTOR, 3)
    p4est_quad_int.Register(ELASTIC_STRAIN_VECTOR, 3)
    p4est_quad_int.Finalize()

    p4est_quad = P4estQuad(p4est_quad_nodal, p4est_quad_int)

    # create the P4est model
    p4est_util = P4estUtilities()
    p4est_model = p4est_util.CreateModel(p4est_order, mpi.world, p4est_quad)
    p4est_model.Initialize(model.model_part)

    mp = ConstructModelPart(p4est_model)
    model.SetModelPart(mp)
    model.InitializeModel()

    #################################################################

    if case_number == 1:  # no refinement at all
        model.material = "mohr-coulomb"
        whole_model_refinement = 0
        max_refine_level = 0
    elif case_number == 2:  # global refinement once
        model.material = "mohr-coulomb"
        whole_model_refinement = 1
        max_refine_level = 0
    elif case_number == 3:  # global refinement twice
        model.material = "mohr-coulomb"
        whole_model_refinement = 2
        max_refine_level = 0
    elif case_number == 4:   # global refinement 3 times
        model.material = "mohr-coulomb"
        whole_model_refinement = 3
        max_refine_level = 0
    elif case_number == 5:   # local refinement from the original mesh, refine level = 1
        model.material = "mohr-coulomb"
        whole_model_refinement = 0
        max_refine_level = 1
    elif case_number == 6:   # local refinement from the original mesh, refine level = 2
        model.material = "mohr-coulomb"
        whole_model_refinement = 0
        max_refine_level = 2
    elif case_number == 7:   # local refinement from the original mesh, refine level = 3
        model.material = "mohr-coulomb"
        whole_model_refinement = 0
        max_refine_level = 3
    elif case_number == 101:   # no refinement at all  with drucker-prager
        model.material = "drucker-prager"
        whole_model_refinement = 0
        max_refine_level = 0
    elif case_number == 102:   # local refinement from the original mesh, refine level = 2, with drucker-prager
        whole_model_refinement = 0
        max_refine_level = 2
    else:
        raise Exception("Unknonwn case " + str(case_number))

    ####################################

    if case_number == 1:
        #MohrCoulomb test Without P4est, global refinement = 0, failed at 4.35  , Time = 26 sec
        delta_load_lists = [1.2578125  , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]
    elif case_number == 2:
        #MohrCoulomb test Without P4est, global refinement = 1, failed at 4.17  , Time = 116 sec
        delta_load_lists = [1.07813    , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]
    elif case_number == 3:
        #MohrCoulomb test Without P4est, global refinement = 2, failed at 4.09  , Time = 694 sec
        delta_load_lists = [1.0        , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]
    elif case_number == 4:
        #MohrCoulomb test Without P4est, global refinement = 3, failed at 4.05  , Time =4839 sec
        delta_load_lists = [0.9609375  , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.003125 ]
    elif case_number == 5:
        #MohrCoulomb test With P4est   , max_refine_level = 1, failed at 4.187 , Time = 49 sec
        delta_load_lists = [1.1015675  , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]
    elif case_number == 6:
        #MohrCoulomb test With P4est   , max_refine_level = 2, failed at 4.093 , Time = 118 sec
        delta_load_lists = [1.00390625 , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]
    elif case_number == 7:
        #MohrCoulomb test With P4est   , max_refine_level = 3, failed at 4.05  , Time = 409 sec
        delta_load_lists = [0.9609375  , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.003125 ]
    elif case_number == 101:
        #Druckerprager test Without P4est, global refinement = 0, failed at 5.03  , Time = 23 sec
        delta_load_lists = [1.953125     , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]
    elif case_number == 102:
        #Druckerprager test With P4est   , max_refine_level = 2, failed at 4.73 , Time = 118 sec
        delta_load_lists = [1.6484375    , 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0625, 0.03125, 0.03125, 0.015625, 0.0078125, 0.0078125]

    #################################################################

    balance_option = 0

    for _ in range (whole_model_refinement):
        elements_to_be_refined = []
        for elem_id in range(model.model_part.NumberOfElements()):
            elements_to_be_refined.append(elem_id+1)

        # mark and refine
        p4est_model.Mark(elements_to_be_refined, P4estRefineFlag.TO_REFINE)
        p4est_model.Refine()
        p4est_model.Repartition(balance_option)

        # create the model_part out from p4est
        mp = ConstructModelPart(p4est_model)
        model.SetModelPart(mp)
        model.InitializeModel()

    #################################################################

    refinement_count_per_solve = max_refine_level       # set to zero if you dont want to use p4est
    plasticity_indicator_refinement_threshold = 0.01    # setting a negative value will refine everything
    selection_ratio = 0.01                              # selection criteria for coarsening
    load = 0.0
    tol = 1.0e-6
    # model.WriteOutput( 0)

    for delta_load in delta_load_lists:
        load = load + delta_load
        # #################### P4est transfering data - marking - refining - transfering data to kratos ##########
        # variable_transfer = VariableTransferUtility(SuperLUSolver())
        # variable_transfer.TransferVariablesToNodes(model.model_part, THREED_STRESSES)
        # variable_transfer.TransferVariablesToNodes(model.model_part, PLASTICITY_INDICATOR)
        # variable_transfer.TransferVariablesToNodes(model.model_part, INTEGRATION_POINT_STRAIN_VECTOR)
        # variable_transfer.TransferVariablesToNodes(model.model_part, ELASTIC_STRAIN_VECTOR)
        # p4est_model.SynchronizeBackward()

        p4est_model.SynchronizeBackwardUsingExtrapolation()

        print("##############SOLVING for load = %f#################" % (load))

        for _ in range (refinement_count_per_solve):
            if not coarsening:
                # finding elements that should be refined
                elements_to_be_refined = []
                for elem_id in range(model.model_part.NumberOfElements()):
                    values_list = model.model_part.GetElement(elem_id+1).GetValuesOnIntegrationPoints(PLASTICITY_INDICATOR, model.model_part.ProcessInfo)
                    for value in values_list[0] :
                        if value >= plasticity_indicator_refinement_threshold :
                            if model.model_part.GetElement(elem_id+1).GetValue(P4EST_REFINE_LEVEL) < max_refine_level :
                                elements_to_be_refined.append(elem_id+1)
                                break
            else:
                # compute the maximum plastic indicator
                number_of_plastic_elements = 0
                max_plastic = 0.0
                for element in model.model_part.Elements:
                    values_list = element.GetValuesOnIntegrationPoints(PLASTICITY_INDICATOR, model.model_part.ProcessInfo)
                    plastic = False
                    for value in values_list:
                        if value[0] >= plasticity_indicator_refinement_threshold:
                            if value[0] > max_plastic:
                                max_plastic = value[0]
                            plastic = True
                    if plastic == True:
                        number_of_plastic_elements = number_of_plastic_elements + 1

                # finding elements that should be refined and coarsened
                elements_to_be_refined = []
                elements_to_be_coarsened = []
                for element in model.model_part.Elements:
                    values_list = element.GetValuesOnIntegrationPoints(PLASTICITY_INDICATOR, model.model_part.ProcessInfo)
                    for value in values_list:
                        if value[0] >= plasticity_indicator_refinement_threshold:
                            if value[0] > selection_ratio*max_plastic:
                                if element.GetValue(P4EST_REFINE_LEVEL) < max_refine_level :
                                    elements_to_be_refined.append(element.Id)
                                    break
                            else:
                                if elem.GetValue(P4EST_REFINE_LEVEL) < max_refine_level :
                                    elements_to_be_coarsened.append(element.Id)
                                    break

            # mark and refine
            print(f"Element to be refined at load {load} : {elements_to_be_refined}")
            p4est_model.Mark(elements_to_be_refined, P4estRefineFlag.TO_REFINE)
            p4est_model.Refine()
            if coarsening:
                print(f"Element to be coarsened at load {load} : {elements_to_be_coarsened}")
                p4est_model.Mark(elements_to_be_coarsened, P4estRefineFlag.TO_COARSEN)
                p4est_model.Coarsen()
            p4est_model.Repartition(balance_option)

            # create the model_part out from p4est
            mp = ConstructModelPart(p4est_model)
            model.SetModelPart(mp)
            model.InitializeModel()

            # synchonize the data from forest to the model_part
            p4est_model.SynchronizeForward()

        ### Boundary condition
        tol = 1.0e-6
        for node in model.model_part.Nodes:
            if abs(node.X0 - 0.0) < tol or abs(node.X0 - 75.0) < tol:
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            if abs(node.Y0 - 0.0) < tol:
                node.Fix(DISPLACEMENT_X)
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

        ### Apply Load
        model.model_part.Properties[1].SetValue(DENSITY, 2.038735 * load)
        time = load
        # export the results before solve
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            # export the results after solve
            model.WriteOutput( time)
