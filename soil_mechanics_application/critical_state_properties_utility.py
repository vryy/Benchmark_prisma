import sys
import os
import math

#importing Kratos main library
from Kratos import *

from KratosStructuralApplication import *
from KratosSoilMechanicsApplication import *

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

class CriticalStatePropertiesUtility:
    def __init__( self, matfile ):
        print("initialized CriticalStatePropertiesUtility")
        self.ReadMaterials( matfile )
        self.echo_level = 0
        self.dimension = 3

        # self.isotropic3dpointer = LinearElastic3D()
        self.isotropic3dpointer = Isotropic3D()
        self.planestrainpointer = LinearElasticPlaneStrain()
        self.planestresspointer = LinearElasticPlaneStress()
        self.camclaypointer = CamClay3D()
        self.modifiedcamclayib_pointer = ModifiedCamClay3D_IB()
        self.modifiedcamclayik_pointer = ModifiedCamClay3D_IK()
        self.modifiedcamclayii_pointer = ModifiedCamClay3dImplicit_II()
        self.modifiedcamclayiis_pointer = ModifiedCamClay3dImplicit_II_s()
        self.modifiedcamclayiias_pointer = ModifiedCamClay3dImplicit_II_as()
        self.modifiedcamclayiics_pointer = ModifiedCamClayIIConstantStiffness()
        self.modifiedcamclayiii_pointer = ModifiedCamClay3dImplicit_III()
        self.casm_pointer1 = ClayAndSandExplicit()
        self.casm_pointer2 = ClayAndSandImplicitSubstepping()
        self.casm_pointer3 = ClayAndSandConstantStiffness()
        self.casm_pointer4 = ClayAndSandImplicitAdaptiveSubstepping()

    def SetEchoLevel( self, level ):
        self.echo_level = level

    def SetDimension( self, dim ):
        self.dimension = dim

    def MaterialType( self ):
        return "critical state"

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

    ### Set the material for element. Before calling this subroutine, search_type and mat_type shall be be set
    def SetMaterialProperties( self, model_part, element ):
        integration_points = element.GetIntegrationPointsInReferenceFrame()
        if self.echo_level > 0:
            print("setting material for element " + str(element.Id) + ", number of integration points = " + str(len(integration_points)))

        if len(integration_points) == 0:
            return

        cl_pointers = []
        youngs_moduli = []
        poisson_ratios = []
        densities = []
        K0_values = []
        density = 0.0
        csl_slopes = []
        virgin_compression_indices = []
        swell_indices = []
        void_ratios = []
        spacing_ratios = []
        shape_parameters = []
        shear_modulus_evolutions = []
        associativities = []
        parent_element_indices = []
        integration_point_indices = []
        for counter in range(0, len(integration_points) ):
            if self.search_type == "by_searching":
                mat = self.GetMaterialID(integration_points[counter])
            elif self.search_type == "by_name":
                mat = self.mat_type
                if not(mat in self.materials):
                    raise Exception("Material " + mat + " does not exist in the material file")
            elif self.search_type == "by_searching excluding fictitious":
                # for this mode the brep and fictitious_mat_type has to be set
                mat = self.GetMaterialID(integration_points[counter])
                if not self.brep.IsInside(integration_points[counter]):
                    mat = self.fictitious_mat_type
            else:
                raise Exception("Unknown search_type " + self.search_type)

            if( mat != "dummy" ): # the material name must not be dummy
                material = self.materials[mat]
                if self.echo_level > 1:
                    print("setting material " + mat + ", model_type: " + material.params['model_type'] + " for element " + str(element.Id) + ", point " + str(counter))
                if( (material.params['model_type'] == "isotropic3d") or (material.params['model_type'] == "plane_strain") or (material.params['model_type'] == "plane_stress") ):
                    if( material.params['model_type'] == "isotropic3d" ):
                        cl_pointers.append( self.isotropic3dpointer.Clone() )
                    elif( material.params['model_type'] == "plane_strain" ):
                        cl_pointers.append( self.planestrainpointer.Clone() )
                    elif( material.params['model_type'] == "plane_stress" ):
                        cl_pointers.append( self.planestresspointer.Clone() )
                    youngs_moduli.append( float(self.materials[mat].params['youngs_modulus']) )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    densities.append( float(self.materials[mat].params['density']) )
                    csl_slopes.append( 0.0 ) # not required
                    virgin_compression_indices.append( 0.0 ) # not required
                    swell_indices.append( 0.0 ) # not required
                    void_ratios.append( 0.0 ) # not required
                    spacing_ratios.append( 0.0 ) # not required
                    shape_parameters.append( 0.0 ) # not required
                    shear_modulus_evolutions.append( 0.0 ) # not required
                    associativities.append( 0.0 ) # not required
                elif( (material.params['model_type'] == "cam_clay") or (material.params['model_type'] == "modified_cam_clay_ik") or (material.params['model_type'] == "modified_cam_clay_ib") ):
                    if material.params['model_type'] == "cam_clay":
                        cl_pointers.append( self.camclaypointer.Clone() )
                    elif material.params['model_type'] == "modified_cam_clay_ik":
                        cl_pointers.append( self.modifiedcamclayik_pointer.Clone() )
                    elif material.params['model_type'] == "modified_cam_clay_ib":
                        cl_pointers.append( self.modifiedcamclayib_pointer.Clone() )
                    youngs_moduli.append( 0.0 )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    densities.append( float(self.materials[mat].params['density']) )
                    csl_slopes.append( float(self.materials[mat].params['csl_slope']) )
                    virgin_compression_indices.append( float(self.materials[mat].params['virgin_compression_index']) )
                    swell_indices.append( float(self.materials[mat].params['swell_index']) )
                    void_ratios.append( float(self.materials[mat].params['void_ratio']) )
                    spacing_ratios.append( 0.0 ) # not required
                    shape_parameters.append( 0.0 ) # not required
                    shear_modulus_evolutions.append( 0.0 ) # not required
                    associativities.append( 0.0 ) # not required
                elif( "modified_cam_clay_ii" in material.params['model_type'] ):
                    if( material.params['model_type'] == "modified_cam_clay_ii" ):
                        cl_pointers.append( self.modifiedcamclayii_pointer.Clone() )
                    elif( material.params['model_type'] == "modified_cam_clay_ii_s" ):
                        cl_pointers.append( self.modifiedcamclayiis_pointer.Clone() )
                    elif( material.params['model_type'] == "modified_cam_clay_ii_as" ):
                        cl_pointers.append( self.modifiedcamclayiias_pointer.Clone() )
                    elif( material.params['model_type'] == "modified_cam_clay_ii_constant_stiffness" ):
                        cl_pointers.append( self.modifiedcamclayiics_pointer.Clone() )
                    youngs_moduli.append( 0.0 )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    densities.append( float(self.materials[mat].params['density']) )
                    csl_slopes.append( float(self.materials[mat].params['csl_slope']) )
                    virgin_compression_indices.append( float(self.materials[mat].params['virgin_compression_index']) )
                    swell_indices.append( float(self.materials[mat].params['swell_index']) )
                    void_ratios.append( float(self.materials[mat].params['void_ratio']) )
                    spacing_ratios.append( 0.0 ) # not required
                    shape_parameters.append( 0.0 ) # not required
                    shear_modulus_evolutions.append( 0.0 ) # not required
                    associativities.append( 0.0 ) # not required
                elif( material.params['model_type'] == "modified_cam_clay_iii" ):
                    cl_pointers.append( self.modifiedcamclayiii_pointer.Clone() )
                    youngs_moduli.append( 0.0 )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    densities.append( float(self.materials[mat].params['density']) )
                    csl_slopes.append( float(self.materials[mat].params['csl_slope']) )
                    virgin_compression_indices.append( float(self.materials[mat].params['virgin_compression_index']) )
                    swell_indices.append( float(self.materials[mat].params['swell_index']) )
                    void_ratios.append( 0.0 ) # not required
                    spacing_ratios.append( 0.0 ) # not required
                    shape_parameters.append( 0.0 ) # not required
                    shear_modulus_evolutions.append( float(self.materials[mat].params['shear_modulus_evolution']) )
                    associativities.append( float(self.materials[mat].params['associativity']) )
                elif( "casm" in material.params['model_type'] ):
                    if( material.params['model_type'] == "casm_explicit" ):
                        cl_pointers.append( self.casm_pointer1.Clone() )
                    elif( material.params['model_type'] == "casm_implicit_substepping" ):
                        cl_pointers.append( self.casm_pointer2.Clone() )
                    elif( material.params['model_type'] == "casm_implicit_adaptive_substepping" ):
                        cl_pointers.append( self.casm_pointer4.Clone() )
                    elif( material.params['model_type'] == "casm_constant_stiffness" ):
                        cl_pointers.append( self.casm_pointer3.Clone() )
                    youngs_moduli.append( 0.0 )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    densities.append( float(self.materials[mat].params['density']) )
                    csl_slopes.append( float(self.materials[mat].params['csl_slope']) )
                    virgin_compression_indices.append( float(self.materials[mat].params['virgin_compression_index']) )
                    swell_indices.append( float(self.materials[mat].params['swell_index']) )
                    void_ratios.append( float(self.materials[mat].params['void_ratio']) )
                    spacing_ratios.append( float(self.materials[mat].params['spacing_ratio']) )
                    shape_parameters.append( float(self.materials[mat].params['shape_parameter']) )
                    shear_modulus_evolutions.append( 0.0 ) # not required
                    associativities.append( 0.0 ) # not required
                else:
                    raise Exception("ERROR: Material " + material.params['model_type'] + " not defined")
            density = float(self.materials[mat].params['density'])
            porosity = float(self.materials[mat].params['porosity'])
            permeability = float(self.materials[mat].params['permeability'])
            parent_element_indices.append(element.Id)
            integration_point_indices.append(counter)
        if( len(cl_pointers) != len(integration_points) ):
            print("ERROR: not all points are defined")
            sys.exit(0)
        for cl in cl_pointers:
            cl.Set(ACTIVE, True)
        element.SetValuesOnIntegrationPoints( CONSTITUTIVE_LAW, cl_pointers, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( K0, K0_values, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( YOUNG_MODULUS, youngs_moduli, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( POISSON_RATIO, poisson_ratios, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( CSL_SLOPE, csl_slopes, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( VIRGIN_COMPRESSION_INDEX, virgin_compression_indices, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( SWELL_INDEX, swell_indices, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( VOID_RATIO, void_ratios, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( SPACING_RATIO, spacing_ratios, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( SHAPE_PARAMETER, shape_parameters, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( SHEAR_MODULUS_EVOLUTION, shear_modulus_evolutions, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( ASSOCIATIVITY, associativities, model_part.ProcessInfo )
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
        else:
            raise Exception("Unknown search_type " + self.search_type)

        if 'density' in self.materials[mat].params:
            density = float(self.materials[mat].params['density'])
            element.SetValue( DENSITY, density )

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
        if 'fix_void_ratio' in self.materials[mat].params:
            if self.materials[mat].params['fix_void_ratio'] == "true":
                element.Properties.SetValue(FIX_VOID_RATIO, True)
            else:
                element.Properties.SetValue(FIX_VOID_RATIO, False)

        if 'thickness' in self.materials[mat].params:
            element.Properties.SetValue(THICKNESS, float(self.materials[mat].params['thickness']))

        if 'local_error_tolerance' in self.materials[mat].params:
            element.Properties.SetValue(LOCAL_ERROR_TOLERANCE, float(self.materials[mat].params['local_error_tolerance']))

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
            permeability_28_days = float(self.materials[mat].params['permeability_transition'])
            element.SetValue( PERMEABILITY_TRANSITION, permeability_28_days )

        ####

        element.ResetConstitutiveLaw()
        if self.echo_level > 0:
            print("element " + str(element.Id) + " is set with " + str(len(cl_pointers)) + " integration points")

    def SetIsotropicPrestress( self, element, prestress, process_info ):
        integration_points = element.GetIntegrationPoints()
        prestresses = []
        for point in integration_points:
            prestresses.append([prestress, prestress, prestress, 0, 0, 0]);
        element.SetValuesOnIntegrationPoints( PRESTRESS, prestresses, 6, process_info )

    def SetPrestressFactor( self, element, factor, process_info ):
        integration_points = element.GetIntegrationPoints()
        prestress_factors = []
        for point in integration_points:
            prestress_factors.append(factor);
        element.SetValuesOnIntegrationPoints( PRESTRESS_FACTOR, prestress_factors, process_info )

    def SetPreconsolidationPressure2D( self, element, ocr, process_info ):
        prestress = element.GetValuesOnIntegrationPoints( PRESTRESS, process_info )
        ocr_pressures = []
        inactive_indices = []
        cnt = 0
        for pres in prestress:
            if pres[1] < 0.0:
                inactive_indices.append(cnt)
                cnt = cnt + 1
                ocr_pressures.append(0.0)
            else:
                ocr_pressures.append(ocr*(pres[1]))
        element.SetValuesOnIntegrationPoints( PRECONSOLIDATION_PRESSURE, ocr_pressures, process_info )
        if len(inactive_indices) > 0:
            print("!!!WARNING!!! at element " + str(element.Id) + ": " + str(len(inactive_indices)) + " is deactivated")
            cl_pointers = element.GetValuesOnIntegrationPoints( CONSTITUTIVE_LAW, process_info )
            for i in inactive_indices:
                cl_pointers[i].Set(ACTIVE, False)

    def SetPreconsolidationPressure( self, element, ocr, process_info ):
        prestress = element.GetValuesOnIntegrationPoints( PRESTRESS, process_info )
        ocr_pressures = []
        inactive_indices = []
        cnt = 0
        for pres in prestress:
            if pres[2] < 0.0:
                inactive_indices.append(cnt)
                cnt = cnt + 1
                ocr_pressures.append(0.0)
                print("!!!WARNING!!! at element " + str(element.Id) + " negative preconsolidation pressure is detected")
            else:
                ocr_pressures.append(ocr*(pres[2]))
        element.SetValuesOnIntegrationPoints( PRECONSOLIDATION_PRESSURE, ocr_pressures, process_info )
        if len(inactive_indices) > 0:
            print("!!!WARNING!!! at element " + str(element.Id) + ": " + str(len(inactive_indices)) + " is deactivated")
            cl_pointers = element.GetValuesOnIntegrationPoints( CONSTITUTIVE_LAW, process_info )
            for i in inactive_indices:
                cl_pointers[i].Set(ACTIVE, False)

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

#####BEGIN REPORTING
def Report1(ifile, model, prestress):
    strains = model.model_part.Elements[1].GetValuesOnIntegrationPoints(STRAIN, model.model_part.ProcessInfo)
    stresses = model.model_part.Elements[1].GetValuesOnIntegrationPoints(STRESSES, model.model_part.ProcessInfo)
    pcs = model.model_part.Elements[1].GetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, model.model_part.ProcessInfo)
    e_xx = strains[0][0]
    e_yy = strains[0][1]
    e_zz = strains[0][2]
    e_xy = strains[0][3]
    e_yz = strains[0][4]
    e_xz = strains[0][5]
    e_p = e_xx + e_yy + e_zz
    e_q = math.sqrt(math.pow(e_xx-e_p/3, 2) + math.pow(e_yy-e_p/3, 2) + math.pow(e_zz-e_p/3, 2) + 2.0*(math.pow(e_xy, 2) + math.pow(e_yz, 2) + math.pow(e_xz, 2)))

#    print("stresses[0]: ", stresses[0])
    o_xx = stresses[0][0]
    o_yy = stresses[0][1]
    o_zz = stresses[0][2]
    o_xy = stresses[0][3]
    o_yz = stresses[0][4]
    o_xz = stresses[0][5]
    p_o = (o_xx+ o_yy + o_zz)/3
    p = p_o - prestress
    q = math.sqrt(1.5) * math.sqrt(math.pow(o_xx-p_o, 2) + math.pow(o_yy-p_o, 2) + math.pow(o_zz-p_o, 2) + 2.0*(math.pow(o_xy, 2) + math.pow(o_yz, 2) + math.pow(o_xz, 2)))

    pc = pcs[0][0]

#    psigma = model.model_part.Elements[1].GetValuesOnIntegrationPoints(PRINCIPAL_STRESS, model.model_part.ProcessInfo)
#    o_1 = psigma[0][0]
#    o_2 = psigma[0][1]
#    o_3 = psigma[0][2]
    if o_xx >= o_yy and o_xx >= o_zz:
        o_1 = o_xx
        if o_yy >= o_zz:
            o_2 = o_yy
            o_3 = o_zz
        else:
            o_2 = o_zz
            o_3 = o_yy
    elif o_yy >= o_xx and o_yy >= o_zz:
        o_1 = o_yy
        if o_xx >= o_zz:
            o_2 = o_xx
            o_3 = o_zz
        else:
            o_2 = o_zz
            o_3 = o_xx
    elif o_zz >= o_xx and o_zz >= o_yy:
        o_1 = o_zz
        if o_xx >= o_yy:
            o_2 = o_xx
            o_3 = o_yy
        else:
            o_2 = o_yy
            o_3 = o_xx

#    ifile.write(str(e_yy) + "\t" + str(e_p) + "\t" + str(e_q) + "\t" + str(-p) + "\t" + str(q) + "\n")
    ifile.write(str(e_zz) + "\t" + str(e_p) + "\t" + str(e_q) + "\t" + str(-p) + "\t" + str(q) + "\t" + str(o_1-o_3) + "\t" + str(pc) + "\n")

def Report2(ifile, element, pidx, prestress, process_info):
    strains = element.GetValuesOnIntegrationPoints(STRAIN, process_info)
    stresses = element.GetValuesOnIntegrationPoints(STRESSES, process_info)
    pcs = element.GetValuesOnIntegrationPoints(PRECONSOLIDATION_PRESSURE, process_info)
    e_xx = strains[pidx][0]
    e_yy = strains[pidx][1]
    e_zz = strains[pidx][2]
    e_xy = strains[pidx][3]
    e_yz = strains[pidx][4]
    e_xz = strains[pidx][5]
    e_p = e_xx + e_yy + e_zz
    e_q = math.sqrt(math.pow(e_xx-e_p/3, 2) + math.pow(e_yy-e_p/3, 2) + math.pow(e_zz-e_p/3, 2) + 2.0*(math.pow(e_xy, 2) + math.pow(e_yz, 2) + math.pow(e_xz, 2)))

    o_xx = stresses[pidx][0]
    o_yy = stresses[pidx][1]
    o_zz = stresses[pidx][2]
    o_xy = stresses[pidx][3]
    o_yz = stresses[pidx][4]
    o_xz = stresses[pidx][5]
    p_o = (o_xx+ o_yy + o_zz)/3
    p = p_o - prestress
    q = math.sqrt(1.5) * math.sqrt(math.pow(o_xx-p_o, 2) + math.pow(o_yy-p_o, 2) + math.pow(o_zz-p_o, 2) + 2.0*(math.pow(o_xy, 2) + math.pow(o_yz, 2) + math.pow(o_xz, 2)))

    pc = pcs[pidx][0]

    ifile.write(str(e_xx) + "\t" + str(e_yy) + "\t" + str(e_zz) + "\t" + str(e_p) + "\t" + str(e_q) + "\t" + str(-p) + "\t" + str(q) + "\t" + str(pc) + "\n")
#####END REPORTING

