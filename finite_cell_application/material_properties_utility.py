import sys
import os
##################################################################
##################################################################

#importing Kratos main library
from Kratos import *
#kernel = Kernel()   #defining kernel

from KratosStructuralApplication import *
from KratosSoilMechanicsApplication import *

class Material:
    def __init__( self, Id, params ):
        self.Id = Id
        self.params = params

class MaterialPropertiesUtility:
    def __init__( self, matfile ):
        print("initialized MaterialPropertiesUtility")
        self.ReadMaterials( matfile )
#        self.isotropic3dpointer = Isotropic3D()
#        self.planestrainpointer = PlaneStrain()
#        self.planestresspointer = PlaneStress()
        self.isotropic3dpointer = LinearElastic3D()
        self.planestrainpointer = LinearElasticPlaneStrain()
        self.planestresspointer = LinearElasticPlaneStress()
        self.druckerpragerpointer = DruckerPragerNonAssociative()
        self.mohrcoulombpointer = MohrCoulomb3dImplicit()
        self.vonmisespointer = VonMises3dImplicitDC()
        self.trescapointer = Tresca3dImplicitDC()

    def ReadMaterials( self, matfile ):
        self.materials = {}
        print("reading materials from file: "+str(matfile))
        matdata = open( matfile, 'r' ).readlines()
        matblock = False
        parametersblock = False
        matname = ""
        parameters = {}
        for line in matdata:
            if "begin material" in line:
                matblock = True
                matname = line.split()[2]
            elif "end material" in line:
                self.materials[matname] = Material(matname, parameters)
                parameters = {}
                matblock = False
            elif "begin parameters" in line and matblock:
                parametersblock = True
            elif "end parameters" in line:
                parametersblock = False
            elif parametersblock and matblock:
                param = line.split()
                parameters[param[0]] = param[1]

    def SetMaterialProperties( self, model_part, element ):
        integration_points = element.GetIntegrationPoints()
        if len(integration_points) == 0:
            return
        cl_pointers = []
        youngs_moduli = []
        poisson_ratios = []
        densities = []
        thicknesses = []
        cohesions = []
        hardening_moduli = []
        K0_values = []
        tensile_strengths = []
        internal_friction_angles = []
        dilatancy_angles = []
        density = 0.0
        for counter in range(0, len(integration_points) ):
            if self.search_type == "by_searching":
                mat = self.GetMaterialID(integration_points[counter])
            elif self.search_type == "by_name":
                mat = self.mat_type
                if not(mat in self.materials):
                    raise Exception("Material " + mat + " does not exist in the material file")
            else:
                raise Exception("Unknown search_type " + self.search_type)

            if( mat != "dummy" ):
                material = self.materials[mat]
                if( material.params['model_type'] == "isotropic3d" or material.params['model_type'] == "plane_strain" or material.params['model_type'] == "plane_stress" ):
                    if( material.params['model_type'] == "isotropic3d" ):
                        cl_pointers.append( self.isotropic3dpointer.Clone() )
                    elif( material.params['model_type'] == "plane_strain" ):
                        cl_pointers.append( self.planestrainpointer.Clone() )
                    elif( material.params['model_type'] == "plane_stress" ):
                        cl_pointers.append( self.planestresspointer.Clone() )
                    youngs_moduli.append( float(self.materials[mat].params['youngs_modulus']) )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    thicknesses.append( float(self.materials[mat].params['thickness']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    tensile_strengths.append( 0.0 )
                    cohesions.append( 0.0 )
                    hardening_moduli.append( 0.0 )
                    internal_friction_angles.append( 0.0 )
                    dilatancy_angles.append( 0.0 )
                    densities.append( float(self.materials[mat].params['density']) )
                elif( material.params['model_type'] == "von_mises" ):
                    cl_pointers.append( self.vonmisespointer.Clone() )
                    youngs_moduli.append( float(self.materials[mat].params['youngs_modulus']) )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    thicknesses.append( float(self.materials[mat].params['thickness']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    tensile_strengths.append( float(self.materials[mat].params['tensile_strength']) )
                    hardening_moduli.append( float(self.materials[mat].params['hardening_modulus']) )
                    cohesions.append( 0.0 )
                    internal_friction_angles.append( 0.0 )
                    dilatancy_angles.append( 0.0 )
                    densities.append( float(self.materials[mat].params['density']) )
                elif( material.params['model_type'] == "tresca" ):
                    cl_pointers.append( self.trescapointer.Clone() )
                    youngs_moduli.append( float(self.materials[mat].params['youngs_modulus']) )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    thicknesses.append( float(self.materials[mat].params['thickness']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    tensile_strengths.append( float(self.materials[mat].params['tensile_strength']) )
                    hardening_moduli.append( float(self.materials[mat].params['hardening_modulus']) )
                    cohesions.append( 0.0 )
                    internal_friction_angles.append( 0.0 )
                    dilatancy_angles.append( 0.0 )
                    densities.append( float(self.materials[mat].params['density']) )
                elif( material.params['model_type'] == "drucker_prager" ):
                    cl_pointers.append( self.druckerpragerpointer.Clone() )
                    youngs_moduli.append( float(self.materials[mat].params['youngs_modulus']) )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    thicknesses.append( float(self.materials[mat].params['thickness']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    tensile_strengths.append( 0.0 )
                    cohesions.append( float(self.materials[mat].params['cohesion']) )
                    hardening_moduli.append( float(self.materials[mat].params['hardening_modulus']) )
                    internal_friction_angles.append( float(self.materials[mat].params['internal_friction_angle']) )
                    dilatancy_angles.append( float(self.materials[mat].params['internal_friction_angle']) ) # by default, only associative plasticity is supported
                    densities.append( float(self.materials[mat].params['density']) )
                elif( material.params['model_type'] == "mohr_coulomb" ):
                    cl_pointers.append( self.mohrcoulombpointer.Clone() )
                    youngs_moduli.append( float(self.materials[mat].params['youngs_modulus']) )
                    poisson_ratios.append( float(self.materials[mat].params['poisson_ratio']) )
                    thicknesses.append( float(self.materials[mat].params['thickness']) )
                    K0_values.append( float(self.materials[mat].params['K0']) )
                    tensile_strengths.append( 0.0 )
                    cohesions.append( float(self.materials[mat].params['cohesion']) )
                    hardening_moduli.append( float(self.materials[mat].params['hardening_modulus']) )
                    internal_friction_angles.append( float(self.materials[mat].params['internal_friction_angle']) )
                    dilatancy_angles.append( float(self.materials[mat].params['dilatancy_angle']) )
                    densities.append( float(self.materials[mat].params['density']) )
                else:
                    print("ERROR: Material not defined")
                    sys.exit(0)
            density = float(self.materials[mat].params['density'])
            porosity = float(self.materials[mat].params['porosity'])
            permeability = float(self.materials[mat].params['permeability'])
        if( len(cl_pointers) != len(integration_points) ):
            print("ERROR: not all points are defined")
            sys.exit(0)
        element.SetValuesOnIntegrationPoints( CONSTITUTIVE_LAW, cl_pointers, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( K0, K0_values, model_part.ProcessInfo )
        element.SetValue( DENSITY, float(density) )
        element.SetValue( POROSITY, float(porosity) )
        element.SetValue( PERMEABILITY_WATER, float(permeability) )
        element.SetValuesOnIntegrationPoints( YOUNG_MODULUS, youngs_moduli, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( POISSON_RATIO, poisson_ratios, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( THICKNESS, thicknesses, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( TENSILE_STRENGTH, tensile_strengths, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( COHESION, cohesions, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( ISOTROPIC_HARDENING_MODULUS, hardening_moduli, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( INTERNAL_FRICTION_ANGLE, internal_friction_angles, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( DILATANCY_ANGLE, dilatancy_angles, model_part.ProcessInfo )
        element.ResetConstitutiveLaw()

    def SetExtrapolatedMaterial( self, model_part, element, cl, representative_weight = 0.0, rule = "explicit" ):
        integration_points = element.GetIntegrationPoints()
        if rule == "implicit":
            extrapolated_constitutive_law_pointer = ExtrapolatedConstitutiveLawImplicit()
        elif rule == "explicit":
            extrapolated_constitutive_law_pointer = ExtrapolatedConstitutiveLawExplicit()
        elif rule == "implex":
            extrapolated_constitutive_law_pointer = ExtrapolatedConstitutiveLawImplex()
        else:
            print("Undefined extrapolated_constitutive_law rule:", rule)
            sys.exit(1)
        cl_pointers = []
        weights = []
        for counter in range(0, len(integration_points) ):
            new_constitutive_law = extrapolated_constitutive_law_pointer.Clone()
            new_constitutive_law.SetValue( CONSTITUTIVE_LAW, cl, model_part.ProcessInfo )
            cl_pointers.append( new_constitutive_law )
            weights.append(representative_weight / len(integration_points))
        element.SetValuesOnIntegrationPoints( CONSTITUTIVE_LAW, cl_pointers, model_part.ProcessInfo )
        element.SetValuesOnIntegrationPoints( REPRESENTATIVE_WEIGHT, weights, model_part.ProcessInfo )
