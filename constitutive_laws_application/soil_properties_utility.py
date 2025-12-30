import sys
import os
import ast

from Kratos import *
from KratosStructuralApplication import *
from KratosSoilMechanicsApplication import *
from KratosConstitutiveLawsApplication import *

def ComputeElementCenter(element):
    elem_center = [0.0, 0.0, 0.0]
    for node in element.GetNodes():
        elem_center[0] = elem_center[0] + node.X0
        elem_center[1] = elem_center[1] + node.Y0
        elem_center[2] = elem_center[2] + node.Z0
    nnodes = len(element.GetNodes())
    elem_center[0] = elem_center[0] / nnodes
    elem_center[1] = elem_center[1] / nnodes
    elem_center[2] = elem_center[2] / nnodes
    return elem_center

class Material:
    def __init__( self, Id, polygon, params ):
        self.Id = Id
        self.polygon = polygon
        self.params = params

    # note: polygons are defined in x-z plane
    def IsInside( self, point ):
        n = len(self.polygon)
        inside = False
        # print(self.polygon)
        p1x,p1y = self.polygon[0]
        for i in range(n+1):
            p2x,p2y = self.polygon[i % n]
            if point[2] > min(p1y,p2y):
                if point[2] <= max(p1y,p2y):
                    if point[0] <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (point[2]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or point[0] <= xinters:
                            inside = not inside
            p1x,p1y = p2x,p2y
        return inside

    # note: polygons are defined in x-z plane
    def IsInside2D( self, point ):
        n = len(self.polygon)
        inside = False
        # print(self.polygon)
        p1x,p1y = self.polygon[0]
        for i in range(n+1):
            p2x,p2y = self.polygon[i % n]
            if point[1] > min(p1y,p2y):
                if point[1] <= max(p1y,p2y):
                    if point[0] <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (point[1]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or point[0] <= xinters:
                            inside = not inside
            p1x,p1y = p2x,p2y
        return inside

class SoilPropertiesUtility:
    """ Utility to assign the material properties for integration point
        If search_type = "by_searching", the search on pylygon will be performed. The polygon is defined in XZ plane.
        If search_type = "by_name", the specific material will be set
        If search_type = "by_prop_id", the Properties Id of the element will govern the material
    """
    def __init__( self, matfile ):
        print("initialized SoilPropertiesUtility")
        self.ReadMaterials( matfile )
        self.echo_level = 0
        self.dimension = 3

        self.umat3e = Umat3e()
        self.udsme = UDSMeImplicit()

        self.aux_util = SoilsAuxiliaryUtility()

    def SetEchoLevel( self, level ):
        self.echo_level = level

    def SetDimension( self, dim ):
        self.dimension = dim

    def SetSearchingBySoilLayer(self, list_layers, layer_dict):
        """ Set the search mode to searching on the soil layer
            The layer surface must be defined by BRep and not intersecting
            The layer_dict shall also include the above and below layer and map the layer identification
            to corresponding material name
            E.g. layer_dict = {-1: "A1", 0: "A1", 2: "A2", 3: "A2"}
        """
        self.soil_layers = list_layers
        self.layer_dict = layer_dict
        self.search_type = "by_searching_on_soil_layer"

        # map from layer name to a unique identifier and vice versa
        all_layers = list(set(self.layer_dict.values()))
        self.layer_prop_id = {}
        self.layer_prop_matname = {}
        it = 1
        for name in all_layers:
            self.layer_prop_id[name] = it
            self.layer_prop_matname[it] = name
            it = it + 1
        print("SoilPropertiesUtility, layer_prop_id: ", self.layer_prop_id)

    def MaterialType( self ):
        return "elastoplastic"

    def ReadMaterials( self, matfile ):
        self.materials = {}
        print("reading materials from file: "+str(matfile))
        matdata = open( matfile, 'r' ).readlines()
        matblock = False
        polygonblock = False
        parametersblock = False
        matname = ""
        coords = []
        parameters = {}
        for line in matdata:
            if "begin material" in line:
                matblock = True
                matname = line.split()[2]
            elif "end material" in line:
                self.materials[matname] = Material(matname, coords, parameters)
                coords = []
                parameters = {}
                matblock = False
            elif "begin polygon" in line and matblock:
                polygonblock = True
            elif "end polygon" in line and matblock:
                polygonblock = False
            elif polygonblock and matblock:
                point_data = line.split()
                coords.append( (float(point_data[0]),float(point_data[2])) )
            elif "begin parameters" in line and matblock and not polygonblock:
                parametersblock = True
            elif "end parameters" in line:
                parametersblock = False
            elif parametersblock and matblock and not polygonblock:
                param = line.split()
                parameters[param[0]] = param[1]
        print("list of materials:")
        for m in self.materials.keys():
            print(" " + str(m))

    ### Set the material for element. Before calling this subroutine, search_type shall be set. In the case search_type="by_searching", mat_type has to also set
    def SetMaterialProperties( self, model_part, element ):
        if self.search_type == "by_searching" or self.search_type == "by_searching excluding fictitious":
            integration_points = element.GetIntegrationPointsInReferenceFrame()
        else:
            integration_points = element.GetIntegrationPointsLocalCoordinates()
        if self.echo_level > 0:
            print("setting material for element " + str(element.Id) + ", number of integration points = " + str(len(integration_points)))

        if len(integration_points) == 0:
            return

        cl_pointers = []
        youngs_moduli = []
        poisson_ratios = []
        densities = []
        cohesions = []
        hardening_moduli = []
        K0_values = []
        internal_friction_angles = []
        dilatancy_angles = []
        transition_angles = []
        apex_control_params = []
        primary_hydration_times = []
        primary_hydration_time_gradients = []
        stiffness_ratios = []

        parent_element_indices = []
        integration_point_indices = []

        for counter in range(0, len(integration_points) ):
            if self.search_type == "by_searching":
                mat = self.GetMaterialID(integration_points[counter])
            elif self.search_type == "by_name":
                mat = self.mat_type
                if not(mat in self.materials):
                    raise Exception("Material " + mat + " does not exist in the material file")
            elif self.search_type == "by_searching_on_soil_layer":
                mat = self.GetMaterialLayer(element, integration_points[counter])
            elif self.search_type == "by_prop_id":
                mat = self.mat_type[element.Properties.Id]
                if not(mat in self.materials):
                    raise Exception("Material " + mat + " does not exist in the material file")
            elif self.search_type == "by_searching excluding fictitious":
                # for this mode the brep and fictitious_mat_type has to be set
                mat = self.GetMaterialID(integration_points[counter])
                if not self.brep.IsInside(integration_points[counter]):
                    mat = self.fictitious_mat_type
            else:
                raise Exception("Unknown search_type " + self.search_type)

            parent_element_indices.append(element.Id)
            integration_point_indices.append(counter)

            if( mat != "dummy" ): # the material name must not be dummy
                material = self.materials[mat]
                if self.echo_level > 1:
                    print("setting material " + mat + ", model_type: " + material.params['model_type'] + " for element " + str(element.Id) + ", point " + str(counter))
                if( material.params['model_type'] == "umat" ):
                    cl_pointers.append( self.umat3e.Clone() )
                    cl = cl_pointers[-1]
                    cl.SetValue(ABAQUS_LIBRARY_NAME, material.params['library_name'], model_part.ProcessInfo)
                    cl.SetValue(UMAT_NAME, material.params['name'], model_part.ProcessInfo)
                    tmp = ast.literal_eval(material.params['umat_params'])
                    mat_params = Vector(len(tmp))
                    for i in range(0, len(tmp)):
                        mat_params[i] = float(tmp[i])
                    # print("mat_params:" + str(mat_params))
                    cl.SetValue(MATERIAL_PARAMETERS, mat_params, model_part.ProcessInfo)
                    cl.SetValue(UMAT_NDI, int(material.params['ndi']), model_part.ProcessInfo)
                    cl.SetValue(UMAT_NSHR, int(material.params['nshr']), model_part.ProcessInfo)
                    cl.SetValue(UMAT_NSTATV, int(material.params['nstatv']), model_part.ProcessInfo)
                    cl.SetValue(UMAT_CMNAME, material.params['cmname'], model_part.ProcessInfo)
                    youngs_moduli.append( 0.0 )
                    poisson_ratios.append( 0.0 )
                    K0_values.append( float(material.params['K0']) )
                    cohesions.append( 0.0 )
                    hardening_moduli.append( 0.0 )
                    internal_friction_angles.append( 0.0 )
                    dilatancy_angles.append( 0.0 )
                    densities.append( float(material.params['density']) )
                    transition_angles.append( 0.0 )
                    apex_control_params.append( 0.0 )
                    primary_hydration_times.append( 0.0 )
                    primary_hydration_time_gradients.append( 0.0 )
                    stiffness_ratios.append( 0.0 )
                elif( material.params['model_type'] == "udsm" ):
                    cl_pointers.append( self.udsme.Clone() )
                    cl = cl_pointers[-1]
                    cl.SetValue(PLAXIS_LIBRARY_NAME, material.params['library_name'], model_part.ProcessInfo)
                    cl.SetValue(USERMOD_NAME, material.params['name'], model_part.ProcessInfo)
                    tmp = ast.literal_eval(material.params['udsm_params'])
                    mat_params = Vector(len(tmp))
                    for i in range(0, len(tmp)):
                        mat_params[i] = float(tmp[i])
                    # print("mat_params:" + str(mat_params))
                    cl.SetValue(MATERIAL_PARAMETERS, mat_params, model_part.ProcessInfo)
                    cl.SetValue(SOIL_MODEL_NUMBER, int(material.params['model_number']), model_part.ProcessInfo)
                    cl.SetValue(IS_UNDRAINED, bool(material.params['is_undr']), model_part.ProcessInfo)
                    cl.SetValue(BULK_W, float(material.params['bulk_w']), model_part.ProcessInfo)
                    youngs_moduli.append( 0.0 )
                    poisson_ratios.append( 0.0 )
                    K0_values.append( float(material.params['K0']) )
                    cohesions.append( 0.0 )
                    hardening_moduli.append( 0.0 )
                    internal_friction_angles.append( 0.0 )
                    dilatancy_angles.append( 0.0 )
                    densities.append( float(material.params['density']) )
                    transition_angles.append( 0.0 )
                    apex_control_params.append( 0.0 )
                    primary_hydration_times.append( 0.0 )
                    primary_hydration_time_gradients.append( 0.0 )
                    stiffness_ratios.append( 0.0 )
                else:
                    print("ERROR: Material not defined")
                    sys.exit(0)
        if( len(cl_pointers) != len(integration_points) ):
            print("len(cl_pointers):", len(cl_pointers))
            print("len(integration_points):", len(integration_points))
            print("ERROR: not all points are defined")
            sys.exit(0)
        if self.echo_level > 2:
            print(" YOUNG_MODULUS: " + str(youngs_moduli))
            print(" POISSON_RATIO: " + str(poisson_ratios))
        element.SetValuesOnIntegrationPoints( CONSTITUTIVE_LAW, cl_pointers, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( K0, K0_values, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( YOUNG_MODULUS, youngs_moduli, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( POISSON_RATIO, poisson_ratios, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( COHESION, cohesions, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( ISOTROPIC_HARDENING_MODULUS, hardening_moduli, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( INTERNAL_FRICTION_ANGLE, internal_friction_angles, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( DILATANCY_ANGLE, dilatancy_angles, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( TRANSITION_ANGLE, transition_angles, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( APEX_CONTROL_PARAMETER, apex_control_params, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( PRIMARY_HYDRATION_TIME, primary_hydration_times, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( PRIMARY_HYDRATION_TIME_GRADIENT, primary_hydration_time_gradients, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( STIFFNESS_RATIO, stiffness_ratios, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( PARENT_ELEMENT_ID, parent_element_indices, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( INTEGRATION_POINT_INDEX, integration_point_indices, model_part.ProcessInfo )

        if self.search_type == "by_searching" or self.search_type == "by_searching excluding fictitious":
            # search for the material
            elem_center = ComputeElementCenter(element)
            mat = self.GetMaterialID(elem_center)
            self.mat_type = mat
        elif self.search_type == "by_name":
            mat = self.mat_type
            if not(mat in self.materials):
                raise Exception("Material " + mat + " does not exist in the material file")
        elif self.search_type == "by_searching_on_soil_layer":
            elem_center = ComputeElementCenter(element)
            mat = self.GetMaterialLayer(element, elem_center)
            self.mat_type = mat
        elif self.search_type == "by_prop_id":
            mat = self.mat_type[element.Properties.Id]
            if not(mat in self.materials):
                raise Exception("Material " + mat + " does not exist in the material file")
        else:
            raise Exception("Unknown search_type " + self.search_type)

        if 'density' in self.materials[mat].params:
            density = float(self.materials[mat].params['density'])
            element.SetValue( DENSITY, density )

        if self.echo_level > 2:
            print(" DENSITY: " + str(density))

        if 'density_water' in self.materials[mat].params:
            density_water = float(self.materials[mat].params['density_water'])
            element.SetValue( DENSITY_WATER, density_water )

        if 'porosity' in self.materials[mat].params:
            porosity = float(self.materials[mat].params['porosity'])
            element.SetValue( POROSITY, float(porosity) )
        if 'fix_porosity' in self.materials[mat].params:
            if self.materials[mat].params['fix_porosity'] == "true":
                element.Properties.SetValue(FIX_POROSITY, True)
                element.SetValue(FIX_POROSITY, True)
            else:
                element.Properties.SetValue(FIX_POROSITY, False)
                element.SetValue(FIX_POROSITY, False)
        else:
            element.SetValue(FIX_POROSITY, True) # default value is True

        if 'thickness' in self.materials[mat].params:
            element.Properties.SetValue(THICKNESS, float(self.materials[mat].params['thickness']))

        if 'permeability' in self.materials[mat].params:
            permeability = float(self.materials[mat].params['permeability'])
            element.SetValue( PERMEABILITY_WATER, permeability )

        if 'permeability_1_day' in self.materials[mat].params:
            permeability_1_day = float(self.materials[mat].params['permeability_1_day'])
            element.SetValue( PERMEABILITY_1_DAY, permeability_1_day )

        if 'permeability_28_days' in self.materials[mat].params:
            permeability_28_days = float(self.materials[mat].params['permeability_28_days'])
            element.SetValue( PERMEABILITY_28_DAYS, permeability_28_days )

        if 'permeability_transition' in self.materials[mat].params:
            permeability_transition = float(self.materials[mat].params['permeability_transition'])
            element.SetValue( PERMEABILITY_TRANSITION, permeability_transition )

        if 'permeability_air' in self.materials[mat].params:
            permeability_air = float(self.materials[mat].params['permeability_air'])
            element.SetValue( PERMEABILITY_AIR, permeability_air )

        if 'density_air' in self.materials[mat].params:
            density_air = float(self.materials[mat].params['density_air'])
            element.SetValue( DENSITY_AIR, density_air )

        if 'bulk_air' in self.materials[mat].params:
            bulk_air = float(self.materials[mat].params['bulk_air'])
            element.SetValue( BULK_AIR, bulk_air )

        if 'swcc' in self.materials[mat].params:
            swcc = self.materials[mat].params['swcc']
        else:
            swcc = "van_genuchten"

        if swcc == "van_genuchten":
            if 'air_entry_value' in self.materials[mat].params:
                air_entry_value = float(self.materials[mat].params['air_entry_value'])
                element.SetValue( AIR_ENTRY_VALUE, air_entry_value )

            if 'first_saturation_param' in self.materials[mat].params:
                first_saturation_param = float(self.materials[mat].params['first_saturation_param'])
                element.SetValue( FIRST_SATURATION_PARAM, first_saturation_param )

            if 'second_saturation_param' in self.materials[mat].params:
                second_saturation_param = float(self.materials[mat].params['second_saturation_param'])
                element.SetValue( SECOND_SATURATION_PARAM, second_saturation_param )

            if 'min_saturation' in self.materials[mat].params:
                min_saturation = float(self.materials[mat].params['min_saturation'])
                element.SetValue( MIN_SATURATION, min_saturation )

            if 'max_saturation' in self.materials[mat].params:
                max_saturation = float(self.materials[mat].params['max_saturation'])
                element.SetValue( MAX_SATURATION, max_saturation )
        elif swcc == "liakopoulous":
            self.aux_util.SetValue(SWCC_LAW, LiakopolousSWCC(), element)
            self.aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, LiakopolousRelativePermeabilityWaterLaw(), element)
            self.aux_util.SetValue(RELATIVE_PERMEABILITY_AIR_LAW, LiakopolousRelativePermeabilityAirLaw(), element)
        else:
            raise Exception(f"Invalid SWCC {swcc}")

        ####

        element.ResetConstitutiveLaw()
        if self.echo_level > 0:
            print("element " + str(element.Id) + " is set with " + str(mat) + " for " + str(len(cl_pointers)) + " integration points")

    ### Get the list of material name at each Gauss point
    def GetMaterialName( self, model_part, element ):
        if self.search_type == "by_searching":
            integration_points = element.GetIntegrationPointsInReferenceFrame()
        else:
            integration_points = element.GetIntegrationPointsLocalCoordinates()
        if self.echo_level > 0:
            print("setting material for element " + str(element.Id) + ", number of integration points = " + str(len(integration_points)))

        if len(integration_points) == 0:
            return []

        mat_name = []
        for counter in range(0, len(integration_points) ):
            if self.search_type == "by_searching":
                mat = self.GetMaterialID(integration_points[counter])
            elif self.search_type == "by_name":
                mat = self.mat_type
                if not(mat in self.materials):
                    raise Exception("Material " + mat + " does not exist in the material file")
            elif self.search_type == "by_searching_on_soil_layer":
                mat = self.GetMaterialLayer(element, integration_points[counter])
            elif self.search_type == "by_prop_id":
                mat = self.mat_type[element.Properties.Id]
                if not(mat in self.materials):
                    raise Exception("Material " + mat + " does not exist in the material file")
            else:
                raise Exception("Unknown search_type " + self.search_type)
            mat_name.append(mat)

        return mat_name


    ### Set the initial stress for K0 procedure. z_list denotes the vertical location of each soil layer
    ### It is assumed that the material layer shall be horizontal
    ### Example:
    ###     z_list = [0.0, -2.0, -25.0, -40.0]
    ###     d_list = [1.95, 2.2, 2.4]
    ###     k0_list = [0.5, 0.6, 0.65]
    def SetInitialStress2D( self, model_part, element, g, d_list, z_list, k0_list ):
        points = element.GetValuesOnIntegrationPoints( INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo )
        stresses = []
        for i in range(0, len(points)):
            # print("setting for point " + str(points[i][1]))
            if (points[i][1] > z_list[0]) or (points[i][1] < z_list[-1]):
                raise Exception("The point lies outside the material range")
            for j in range(1, len(z_list)):
                if points[i][1] >= z_list[j]:
                    r = j-1
                    break
            k0 = k0_list[r]
            density = d_list[r]
            sigmav_top = 0.0
            for j in range(0, r):
                sigmav_top = sigmav_top + g*d_list[j]*(z_list[j]-z_list[j+1])
            # print("k0:", k0)
            # print("density:", density)
            # print("sigmav_top:", sigmav_top)
            initial_stress = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            initial_stress[1] = density*g*(z_list[r]-points[i][1]) + sigmav_top # yy
            initial_stress[0] = k0*initial_stress[1] # xx
            initial_stress[2] = k0*initial_stress[1] # zz
            initial_stress[3] = 0.0 # xy
            initial_stress[4] = 0.0 # yz
            initial_stress[5] = 0.0 # xz
            stresses.append(initial_stress)
        # print("stresses:", stresses)
        element.SetValuesOnIntegrationPoints( INITIAL_STRESS, stresses, 6, model_part.ProcessInfo )
        # sys.exit(0)

    #
    # set the initial stress for K0 procedure
    #
    def SetInitialStress( self, model_part, element, g, z0 ):
        points = element.GetValuesOnIntegrationPoints( INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo )
        stresses = []
        for i in range(0, len(points)):
            # print("setting for point " + str(points[i][1]))
            if (points[i][1] > z_list[0]) or (points[i][1] < z_list[-1]):
                raise Exception("The point lies outside the material range")
            for j in range(1, len(z_list)):
                if points[i][1] >= z_list[j]:
                    r = j-1
                    break
            k0 = k0_list[r]
            density = d_list[r]
            sigmav_top = 0.0
            for j in range(0, r):
                sigmav_top = sigmav_top + g*d_list[j]*(z_list[j]-z_list[j+1])
            initial_stress = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            initial_stress[2] = density*g*(z_list[r]-points[i][1]) + sigmav_top # yy
            initial_stress[0] = k0*initial_stress[2] # xx
            initial_stress[1] = k0*initial_stress[2] # zz
            initial_stress[3] = 0.0 # xy
            initial_stress[4] = 0.0 # yz
            initial_stress[5] = 0.0 # xz
            stresses.append(initial_stress)
        element.SetValuesOnIntegrationPoints( INITIAL_STRESS, stresses, 6, model_part.ProcessInfo )

    #
    # set the preconsolidation pressure based on OCR and vertical stress
    #
    def SetPreconsolidationPressure( self, model_part, element, OCR ):
        prestress = element.GetValuesOnIntegrationPoints( PRESTRESS, model_part.ProcessInfo )
        ocr_pressures = []
        for pres in prestress:
            ocr_pressures.append(OCR*(pres[2]))
        element.SetValuesOnIntegrationPoints( PRECONSOLIDATION_PRESSURE, ocr_pressures, model_part.ProcessInfo )

    #
    # get the material only in the region defined by polygon
    #
    def GetMaterialID( self, point ):
        if self.dimension == 3:
            for material in self.materials:
                if len(self.materials[material].polygon) != 0:
                    if self.materials[material].IsInside( point ):
                        return material
        elif self.dimension == 2:
            for material in self.materials:
                if len(self.materials[material].polygon) != 0:
                    if self.materials[material].IsInside2D( point ):
                        return material
        return("dummy")

    #
    # get the material only in the region defined by soil layer
    #
    def GetMaterialLayer( self, elem, local_point ):
        stat = []
        current_layer = -1 # start first with layer -1
        for layer in self.soil_layers:
            stat = layer.IsInside(elem, local_point, 0) # search in reference configuration

            if stat == False:
                break

            current_layer = current_layer + 1

        return self.layer_dict[current_layer]
