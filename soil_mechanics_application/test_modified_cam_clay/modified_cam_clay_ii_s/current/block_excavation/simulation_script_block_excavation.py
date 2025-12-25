##################################################################
##### simulation script for face stability simulation        #####
##### copyright by Hoang-Giang, Bui                          #####
#####          Ruhr University Bochum                        #####
##### all rights reserved                                    #####
##################################################################
import sys
import os
import math
import six
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
from KratosMultiphysics.EkateAuxiliaryApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# from KratosMultiphysics.MultithreadedSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *

kernel = Kernel()   #defining kernel
##################################################################

def GetValue(params, name, default_value):
    if name in params:
        return params[name]
    else:
        return default_value

def FixPressureNodes( model_part, free_node_list_water, free_node_list_air ):
    for node in model_part.Nodes:
        if (node.IsFixed(WATER_PRESSURE)==0):
            node.Fix(WATER_PRESSURE)
            free_node_list_water.append(node)
        if (node.IsFixed(AIR_PRESSURE)==0):
            node.Fix(AIR_PRESSURE)
            free_node_list_air.append(node)

def FreePressureNodes( free_node_list_water, free_node_list_air ):
    for item in free_node_list_water:
        item.Free(WATER_PRESSURE)
    for item in free_node_list_air:
        item.Free(AIR_PRESSURE)

def ApplyInsituWaterPressure( model_part, y_zero, gravity_y ):
    water_density = 1000.0;
    for node in model_part.Nodes:
        water_pressure= water_density*gravity_y*(y_zero-(node.Y-node.GetSolutionStepValue(DISPLACEMENT_Y, 0)))
        if( water_pressure < 1.0 ):
            water_pressure = 1.0
        node.SetSolutionStepValue(WATER_PRESSURE, water_pressure)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, water_pressure)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, water_pressure)
    for node in model_part.Nodes:
        node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

def SetReferenceWaterPressure( model_part ):
    SetReferenceWaterPressureForElements(model_part.Elements, model_part.ProcessInfo)

def SetReferenceWaterPressureForElements( elements, process_info ):
    for element in elements:
        water_pressures = element.GetValuesOnIntegrationPoints( WATER_PRESSURE, process_info )
        pressure_list = []
        for item in water_pressures:
            pressure_list.append( item[0] )
        element.SetValuesOnIntegrationPoints( REFERENCE_WATER_PRESSURE, pressure_list, process_info )

class BlockExcavationSimulator:
    def __init__( self, params ):
        self.params = params
        sys.path.append(params['path']+'/'+params['name']+'.gid')

        if not "account_for_water" in self.params:
            self.params["account_for_water"] = False

    def Run( self ):
        model_virgin = BlockExcavationSimulator.PrepareVirgin(self.params)
        model1 = BlockExcavationSimulator.PrepareSystem(model_virgin, self.params)
        BlockExcavationSimulator.RunExcavationAnalysis(model_virgin, model1, self.params)

    @staticmethod
    def ApplyBC( model, params, tol ):
        """ apply the displacement boundary condition
        """
        x_min = params['x_min']
        x_max = params['x_max']
        y_min = params['y_min']
        y_max = params['y_max']
        for node in model.model_part.Nodes:
            if abs(node.X0 - x_min) < tol:
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            if abs(node.X0 - x_max) < tol:
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            if abs(node.Y0 - y_min) < tol:
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

    @staticmethod
    def ApplyWaterBC(model, params, tol):
        """ apply the water pressure boundary condition
        """
        x_min = params['x_min']
        x_max = params['x_max']
        y_min = params['y_min']
        y_max = params['y_max']
        # fix the water pressure at front and back side
        for node in model.model_part.Nodes:
            if abs(node.X0 - xmin) < tol:
                node.Fix(WATER_PRESSURE)
                if node.Has(LAGRANGE_WATER_PRESSURE):
                    node.Fix(LAGRANGE_WATER_PRESSURE)
                    node.SetSolutionStepValue(LAGRANGE_WATER_PRESSURE, 0.0)
            if abs(node.X0 - xmax) < tol:
                node.Fix(WATER_PRESSURE)
                if node.Has(LAGRANGE_WATER_PRESSURE):
                    node.Fix(LAGRANGE_WATER_PRESSURE)
                    node.SetSolutionStepValue(LAGRANGE_WATER_PRESSURE, 0.0)

        # fix the water pressure on the top surface
        for node in model.model_part.Nodes:
            if abs(node.Y0 - ymax) < tol:
                node.Fix(WATER_PRESSURE)

    @staticmethod
    def Activate( model, layer_name ):
        if layer_name in model.layer_sets:
            for elem_id in model.layer_sets[layer_name]:
                elem = model.model_part.Elements[elem_id]
                elem.Set(ACTIVE, True)
                elem.SetValue(IS_INACTIVE, False)

    @staticmethod
    def Deactivate( model, layer_name ):
        if layer_name in model.layer_sets:
            for elem_id in model.layer_sets[layer_name]:
                elem = model.model_part.Elements[elem_id]
                elem.Set(ACTIVE, False)
                elem.SetValue(IS_INACTIVE, True)

    @staticmethod
    def ExtractSurfaceSettlement( model, y_top, tol ):
        disp_list = []
        for node in model.model_part.Nodes:
            if abs(node.Y0 - y_top) < tol:
                disp_list.append([node.X0, node.GetSolutionStepValue(DISPLACEMENT_Y)])
        disp_list = sorted(disp_list, key = lambda x: x[0])
        return disp_list

    @staticmethod
    def WriteDisp( filename, disp_list ):
        ifile = open(filename, "w")
        ifile.write("x\tdy\n")
        for v in disp_list:
            ifile.write(str(v[0]) + '\t' + str(v[1]) + '\n')
        ifile.close()

    @staticmethod
    def PrepareVirgin( params ):
        geology_virgin_include = __import__(params['name']+"_include")
        if not os.path.isdir(params['path']+"/"+params['name']+".gid/virgin_results/"):
            os.mkdir(params['path']+"/"+params['name']+".gid/virgin_results/")
        model_virgin = geology_virgin_include.Model(params['name'],params['path']+"/"+params['name']+".gid/",params['path']+"/"+params['name']+".gid/virgin_results/",logging=False)
        model_virgin.InitializeModel()

        isu = InSituStressUtility()

        inital_gravity = Vector(3)
        inital_gravity = model_virgin.model_part.Properties[1].GetValue(GRAVITY)
        print("inital_gravity: " + str(inital_gravity))

        from soil_properties_utility import SoilPropertiesUtility
        spu = SoilPropertiesUtility(params['matfile_virgin'])
        # spu.SetDimension(2)
        # spu.SetEchoLevel(1)
        spu.search_type = "by_searching"
        for element in model_virgin.layer_sets['ground']:
           spu.SetMaterialProperties( model_virgin.model_part, model_virgin.model_part.Elements[element] )
        for layer_name in model_virgin.layer_sets:
            if 'excavation' in layer_name:
                for element in model_virgin.layer_sets[layer_name]:
                   spu.SetMaterialProperties( model_virgin.model_part, model_virgin.model_part.Elements[element] )

        BlockExcavationSimulator.ApplyBC(model_virgin, params, 1.0e-3)

        load_virgin_from_restart = False # turn this on or off
        stop_virgin = False

        if load_virgin_from_restart == True:
            delta_time = 100.0
            time = delta_time * (params['number_of_loading_steps_virgin'] + 1)
            model_virgin.LoadRestartFile(time)
            isu.SetPreStressFromCurrentStress(model_virgin.model_part, model_virgin.model_part.ProcessInfo)
            time = time + delta_time
            model_virgin.WriteOutput(time)
        else:
            ##### INSITU-STRESS MODEL #####
            time = 0.0
            delta_time = 100.0

            number_of_load_steps = params['number_of_loading_steps_virgin']

            #Prescribe Hydrostratic Pressure
            free_node_list_water_virgin = []
            free_node_list_air_virgin = []
            z_coord_of_groundwater_table = params['ground_water_table']
            FixPressureNodes(model_virgin.model_part, free_node_list_water_virgin, free_node_list_air_virgin)
            ApplyInsituWaterPressure(model_virgin.model_part, z_coord_of_groundwater_table, 9.81)

            print("##### FIRST STAGE: SOLVE WITH FULL GRAVITY #####", flush=True)
            for step in range(0, number_of_load_steps):
                g_vect = float(step + 1) / number_of_load_steps * inital_gravity
                print("###############virgin phase: loading gravity vector = " + str(g_vect) + "#############")
                for prop in model_virgin.model_part.Properties:
                    prop.SetValue(GRAVITY, g_vect)
                time = time + delta_time
                if not params['dry_run']:
                    model_virgin.SolveModel(time)
                model_virgin.WriteOutput(time)

            print("##### SECOND STAGE: SOLVE FOR INSITU_STRESS #####", flush=True)
            isu = InSituStressUtility()
            isu.SetPreStressFromCurrentStress(model_virgin.model_part, model_virgin.model_part.ProcessInfo)
            for elem in model_virgin.model_part.Elements:
                elem.ResetConstitutiveLaw()
            time = time + delta_time
            if not params['dry_run']:
                model_virgin.SolveModel(time)
            model_virgin.WriteOutput(time)
            print("~~~~~~~~~~~~~~ STEP DONE: APPLICATION OF INSITU STRESS ~~~~~~~~~~~~~~", flush=True)

            ##Model is initialized with InSitu Stress, the next step is only to check the residual displacements
            print("##### THIRD STAGE: CHECKING #####")
            max_disp = 0.0
            for node in model_virgin.model_part.Nodes:
                for direction in range(0,3):
                    if( abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction])) > max_disp ):
                        max_disp = abs(float(node.GetSolutionStepValue(DISPLACEMENT)[direction]))

            print("~~ STEP DONE (INSITU STRESS) --> residual displacement= "+str(max_disp)+"~~")
            print("")
            print("")
            print("")
            print("#######################################################")
            print("#######################################################")
            print("############### insitu_stress done ####################")
            print("#######################################################")
            print("#######################################################")
            print("")
            print("")
            print("")
#            model_virgin.WriteRestartFile(time)
            isu.SetInSituStressFromCurrentStress(model_virgin.model_part, model_virgin.model_part.ProcessInfo)
            if stop_virgin:
                sys.exit(0)

        model_virgin.spu = spu

        return model_virgin

    @staticmethod
    def TransferPrestress(model_virgin, model1, params):
        print("##### TRANSFERRING INSITU STRESS #####")
        # for node in model_virgin.model_part.Nodes:
        #     print("node " + str(node.Id) + " PRESTRESS: " + str(node.GetSolutionStepValue(PRESTRESS)))
        vtu = VariableTransferUtility(MKLPardisoSolver())
        if params['transfer_method'] == "identical":
            vtu.TransferPrestressIdentically( model_virgin.model_part, model1.model_part )
        elif params['transfer_method'] == "normal":
#            vtu.TransferPrestress( model_virgin.model_part, model1.model_part )
            vtu.TransferSpecificVariableWithComponents( model_virgin.model_part, model1.model_part, PRESTRESS, 6 )

#        isu.SetInSituStressFromCurrentStress(model1.model_part, model1.model_part.ProcessInfo)
        model1.WriteOutput(0.0)

        soil_type = params['soil_type']
        if soil_type == "critical state":
            print("##### SETTING PRECONSOLIDATION PRESSURE #####")
            ocr_top = params['OCR_top']
            ocr_bottom = params['OCR_bottom']
            model1.model_part.Properties[1].SetValue(OVERCONSOLIDATION_RATIO, ocr_top)
            # model1.model_part.Properties[1].SetValue(STRESS_UNIT_SCALE, 1.0) # for unit conversion: Pa -> kPa
            model1.model_part.Properties[2].SetValue(OVERCONSOLIDATION_RATIO, ocr_bottom)
            # model1.model_part.Properties[2].SetValue(STRESS_UNIT_SCALE, 1.0) # for unit conversion: Pa -> kPa

            ## reset the material one more time to account for new information
            for element in model1.soil_elems:
                element.ResetConstitutiveLaw()

    @staticmethod
    def PrepareSystem( params ):
        print("##### SET UP MODEL #####")

        viscous_damping = GetValue(params, 'viscous_damping', 0.0)
        abs_tol = GetValue(params, 'abs_tol', 1e-13)
        rel_tol = GetValue(params, 'rel_tol', 1e-13)
        local_error_tolerance = GetValue(params, 'local_error_tolerance', 1e-6)

        system_include = __import__(params['name']+"_include")
        model1 = system_include.Model(params['name'],params['path']+"/"+params['name']+".gid/",params['path']+"/"+params['name']+".gid/",logging=params['logging'],abs_tol=abs_tol,rel_tol=rel_tol)
        model1.InitializeModel(sub_steps_range=params['sub_steps_range'], viscous_damping=viscous_damping, local_error_tolerance=local_error_tolerance)

        ##### INITIALIZE SOIL PROPERTIES UTILITY #####
        soil_type = params['soil_type']
        if soil_type == "critical state":
            import critical_state_properties_utility
            from critical_state_properties_utility import CriticalStatePropertiesUtility as SystemMaterialPropertiesUtility
        elif soil_type == "elastoplastic":
            import soil_properties_utility
            from soil_properties_utility import SoilPropertiesUtility as SystemMaterialPropertiesUtility
        else:
            print("Invalid soil_type: " + str(soil_type))
            sys.exit(0)
        spu = SystemMaterialPropertiesUtility(params['matfile'])
        spu.search_type = "by_searching"

        ##### BOUNDARY CONDITION #####
        BlockExcavationSimulator.ApplyBC(model1, params, 1.0e-3)

        ##### SET UP SOIL PROPERTIES #####
        model1.soil_elems = []
        for element in model1.layer_sets['ground']:
            spu.SetMaterialProperties( model1.model_part, model1.model_part.Elements[element] )
            model1.soil_elems.append( model1.model_part.Elements[element] )
        for layer_name in model1.layer_sets:
            if 'excavation' in layer_name:
                for element in model1.layer_sets[layer_name]:
                    spu.SetMaterialProperties( model1.model_part, model1.model_part.Elements[element] )
                    model1.soil_elems.append( model1.model_part.Elements[element] )

        model1.spu = spu

        ##### SPECIAL FLAGS #####

        if 'first_yielding_compute_mode' in params:
            first_yielding_compute_mode = params['first_yielding_compute_mode']
            model1.model_part.Properties[1].SetValue(FIRST_YIELDING_COMPUTE_MODE, first_yielding_compute_mode)
            model1.model_part.Properties[2].SetValue(FIRST_YIELDING_COMPUTE_MODE, first_yielding_compute_mode)

        ###################################

        return model1

    @staticmethod
    def RunExcavationAnalysis( model_virgin, model1, params ):

        ##################################################################
        ## Transfer stress ###############################################
        ##################################################################

        BlockExcavationSimulator.TransferPrestress(model_virgin, model1, params)

        ###### SOLVE THE FIRST TIME ###### to solve for initial zero state

        #Prescribe Hydrostratic Pressure
        free_node_list_water1 = []
        free_node_list_air1 = []
        z_coord_of_groundwater_table = params['ground_water_table']
        FixPressureNodes(model1.model_part, free_node_list_water1, free_node_list_air1)
        ApplyInsituWaterPressure(model1.model_part, z_coord_of_groundwater_table, 9.81)
        SetReferenceWaterPressure(model1.model_part)

        time = 100.0

        if not params['dry_run']:
            model1.SolveModel(time)
        if params['output']:
            model1.WriteOutput(time)

        ##################################################################
        ## DEACTIVATION ##################################################
        ##################################################################

        ## deactivate the soil
        delta_time = params['time_excavation']
        for step in range(0, params['number_of_steps']):
            for i in range(0, params['number_of_excavation_layers_per_step']):
                excavation_layer = "excavation_" + str(step*params['number_of_excavation_layers_per_step'] + i + 1)
                BlockExcavationSimulator.Deactivate(model1, excavation_layer)

            time = time + delta_time
            if not params['dry_run']:
                model1.SolveModel(time)
            if params['output']:
                model1.WriteOutput(time)

        print("EXCAVATION ANALYSIS COMPLETED.")
