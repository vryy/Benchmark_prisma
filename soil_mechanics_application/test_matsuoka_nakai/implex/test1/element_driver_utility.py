##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import math
#import numpy
import csv
import sys
import os
import time as time_module
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *

"""
This utility provides function to perform load control and displacement control uniaxial, biaxial and triaxial tests
"""

def free_displacement_bc(model):
	""" free displacement boundary conditions in all the nodes"""
	for i in range(1, len(model.model_part.Nodes)+1):
	    model.model_part.Nodes[i].SetSolutionStepValue(DISPLACEMENT_X, 0.0 )
	    model.model_part.Nodes[i].Free(DISPLACEMENT_X )
	    model.model_part.Nodes[i].SetSolutionStepValue(DISPLACEMENT_Y, 0.0 )
	    model.model_part.Nodes[i].Free(DISPLACEMENT_Y )
	    model.model_part.Nodes[i].SetSolutionStepValue(DISPLACEMENT_Z, 0.0 )
	    model.model_part.Nodes[i].Free(DISPLACEMENT_Z )

def assing_layers(model,tol,L):
	""" create layers: """
	model.layer_nodes_sets['fix_X']=[]
	model.layer_nodes_sets['fix_Y']=[]
	model.layer_nodes_sets['fix_Z']=[]
	model.layer_nodes_sets['surface_X']=[]
	model.layer_nodes_sets['surface_Y']=[]
	model.layer_nodes_sets['surface_Z']=[]
	for i in range(1, len(model.model_part.Nodes)+1):
	        x0= model.model_part.Nodes[i].X0
	        y0= model.model_part.Nodes[i].Y0
	        z0= model.model_part.Nodes[i].Z0
	        if( x0 < tol and x0 > -tol):
	              model.layer_nodes_sets['fix_X'].append(i)
	        if( y0 < tol and y0 > -tol):
	              model.layer_nodes_sets['fix_Y'].append(i)
	        if( z0 < tol and z0 > -tol):
	              model.layer_nodes_sets['fix_Z'].append(i)
	        if( x0 < L+tol and x0 > L-tol):
	              model.layer_nodes_sets['surface_X'].append(i)
	        if( y0 < L+tol and y0 > L-tol):
	              model.layer_nodes_sets['surface_Y'].append(i)
	        if( z0 < L+tol and z0 > L-tol):
	              model.layer_nodes_sets['surface_Z'].append(i)


def assing_fix_displacement_bc(model):
	""" assigning displament bc: """
	for i in range(0, len(model.layer_nodes_sets['fix_X'])):
	        model.model_part.Nodes[model.layer_nodes_sets['fix_X'][i]].SetSolutionStepValue(DISPLACEMENT_X, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_X'][i]].SetSolutionStepValue(DISPLACEMENT_NULL_X, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_X'][i]].SetSolutionStepValue(DISPLACEMENT_EINS_X, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_X'][i]].Fix(DISPLACEMENT_X )


	for i in range(0, len(model.layer_nodes_sets['fix_Y'])):
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Y'][i]].SetSolutionStepValue(DISPLACEMENT_Y, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Y'][i]].SetSolutionStepValue(DISPLACEMENT_NULL_Y, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Y'][i]].SetSolutionStepValue(DISPLACEMENT_EINS_Y, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Y'][i]].Fix(DISPLACEMENT_Y )

	for i in range(0, len(model.layer_nodes_sets['fix_Z'])):
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Z'][i]].SetSolutionStepValue(DISPLACEMENT_Z, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Z'][i]].SetSolutionStepValue(DISPLACEMENT_NULL_Z, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Z'][i]].SetSolutionStepValue(DISPLACEMENT_EINS_Z, 0.0 )
	        model.model_part.Nodes[model.layer_nodes_sets['fix_Z'][i]].Fix(DISPLACEMENT_Z )

def assign_initial_prestress(model,initial_prestress):
	""" an initial prestress is assigned to the gauss points """
	vector_size = 6
	# Get the PRESTRESS computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the PRESTRESS on each integration points:
		prestresses = element.GetValuesOnIntegrationPoints( PRESTRESS, model.model_part.ProcessInfo )
		prestress_list = []
		for item in prestresses:
			initial_prestress_list = []
			initial_prestress_list.append(initial_prestress["Sxx"])
			initial_prestress_list.append(initial_prestress["Syy"])
			initial_prestress_list.append(initial_prestress["Szz"])
			initial_prestress_list.append(initial_prestress["Sxy"])
			initial_prestress_list.append(initial_prestress["Syz"])
			initial_prestress_list.append(initial_prestress["Sxz"])
			prestress_list.append(initial_prestress_list)
		# Set the initial PRESTRESS on each integration points:
		element.SetValuesOnIntegrationPoints( PRESTRESS, prestress_list, vector_size,model.model_part.ProcessInfo )

def assign_isotropic_pressure(model,p):
	""" isotropic pressure is applied to each face """
	for cond in model.model_part.Conditions:
		if cond.Has(POSITIVE_FACE_PRESSURE):
			cond.SetValue(POSITIVE_FACE_PRESSURE,p)


def assign_sigma_x_stress(model,sigma_x):
	""" pressure is applied to the top face perpendicular to the x face"""
	for cond in model.model_part.Conditions:
		if cond.Has(POSITIVE_FACE_PRESSURE):
			if cond.Id == 2: # x face id
				cond.SetValue(POSITIVE_FACE_PRESSURE,sigma_x)

def assign_sigma_y_stress(model,sigma_y):
	""" pressure is applied to the top face perpendicular to the y face"""
	for cond in model.model_part.Conditions:
		if cond.Has(POSITIVE_FACE_PRESSURE):
			if cond.Id == 1: # y face id
				cond.SetValue(POSITIVE_FACE_PRESSURE,sigma_y)

def assign_sigma_z_stress(model,sigma_z):
	""" pressure is applied to the top face perpendicular to the z face"""
	for cond in model.model_part.Conditions:
		if cond.Has(POSITIVE_FACE_PRESSURE):
			if cond.Id == 3: # z face id
				cond.SetValue(POSITIVE_FACE_PRESSURE,sigma_z)


def assing_fix_displacement_on_surface_x(model,u_x):
        """ assigning displament on surface x: """
        for i in range(0, len(model.layer_nodes_sets['surface_X'])):
                model.model_part.Nodes[model.layer_nodes_sets['surface_X'][i]].SetSolutionStepValue(DISPLACEMENT_X, u_x )
                model.model_part.Nodes[model.layer_nodes_sets['surface_X'][i]].SetSolutionStepValue(DISPLACEMENT_NULL_X, u_x )
                model.model_part.Nodes[model.layer_nodes_sets['surface_X'][i]].SetSolutionStepValue(DISPLACEMENT_EINS_X, u_x )
                model.model_part.Nodes[model.layer_nodes_sets['surface_X'][i]].Fix(DISPLACEMENT_X )


def assing_fix_displacement_on_surface_y(model,u_y):
	""" assigning displament on surface y: """
	for i in range(0, len(model.layer_nodes_sets['surface_Y'])):
	        model.model_part.Nodes[model.layer_nodes_sets['surface_Y'][i]].SetSolutionStepValue(DISPLACEMENT_Y, u_y )
	        model.model_part.Nodes[model.layer_nodes_sets['surface_Y'][i]].SetSolutionStepValue(DISPLACEMENT_NULL_Y, u_y )
	        model.model_part.Nodes[model.layer_nodes_sets['surface_Y'][i]].SetSolutionStepValue(DISPLACEMENT_EINS_Y, u_y )
	        model.model_part.Nodes[model.layer_nodes_sets['surface_Y'][i]].Fix(DISPLACEMENT_Y )

def assing_fix_displacement_on_surface_z(model,u_z):
        """ assigning displament on surface z: """
        for i in range(0, len(model.layer_nodes_sets['surface_Z'])):
                model.model_part.Nodes[model.layer_nodes_sets['surface_Z'][i]].SetSolutionStepValue(DISPLACEMENT_Z, u_z )
                model.model_part.Nodes[model.layer_nodes_sets['surface_Z'][i]].SetSolutionStepValue(DISPLACEMENT_NULL_Z, u_z )
                model.model_part.Nodes[model.layer_nodes_sets['surface_Z'][i]].SetSolutionStepValue(DISPLACEMENT_EINS_Z, u_z )
                model.model_part.Nodes[model.layer_nodes_sets['surface_Z'][i]].Fix(DISPLACEMENT_Z )


def compute_pressure(sigma_x,sigma_y,sigma_z):
	return (sigma_x+sigma_y+sigma_z)/3.0


def make_implex_implicit(model):
	for element in model.model_part.Elements:
		#-> Get the integration points:
		integration_points = element.GetIntegrationPointsLocalCoordinates()
		item_list = []
		for item in range(len(integration_points)):
			item_list.append(True)
		if item == []:
			print("VARIABLE FORCE_IMPLICIT NOT REGISTERED IN THE ELEMENT")
		else:
			#-> Set the flag on each integration points:
			element.SetValuesOnIntegrationPoints( FORCE_IMPLICIT, item_list, model.model_part.ProcessInfo )


def set_temperature(model,temperature):
	for element in model.model_part.Elements:
		#-> Get the integration points:
		integration_points = element.GetIntegrationPointsLocalCoordinates()
		item_list = []
		for item in range(len(integration_points)):
			item_list.append(temperature)
		if item == []:
			print("VARIABLE TEMPERATURE NOT REGISTERED IN THE ELEMENT")
		else:
			#-> Set the flag on each integration points:
			element.SetValuesOnIntegrationPoints( TEMPERATURE, item_list, model.model_part.ProcessInfo )

# post-processing variables

def get_log_equivalent_volumetric_stress(model):
	""" it assumed that the output is returned in Pa"""
	""""the logarithm of the pressure is computed after scale the value to kPa"""
	unit_scale_kpa = 1000
	output_variable = 0.0
	for element in model.model_part.Elements:
		output_variable_list = element.GetValuesOnIntegrationPoints( PRESSURE_P, model.model_part.ProcessInfo )
		output_variable = output_variable_list[0][0]
	return numpy.log(float(output_variable)/unit_scale_kpa)

def get_equivalent_volumetric_stress(model):
	output_variable = 0.0
	for element in model.model_part.Elements:
		output_variable_list = element.GetValuesOnIntegrationPoints( PRESSURE_P, model.model_part.ProcessInfo )
		output_variable = output_variable_list[0][0]
	return float(output_variable)

def get_equivalent_deviatoric_stress(model):
	output_variable = 0.0
	for element in model.model_part.Elements:
		output_variable_list = element.GetValuesOnIntegrationPoints( PRESSURE_Q, model.model_part.ProcessInfo )
		output_variable = output_variable_list[0][0]
	return float(output_variable)

def get_equivalent_volumetric_strain(model):
	output_variable = 0.0
	for element in model.model_part.Elements:
		output_variable_list = element.GetValuesOnIntegrationPoints( EQUIVALENT_VOLUMETRIC_STRAIN, model.model_part.ProcessInfo )
		output_variable = output_variable_list[0][0]
	return -float(output_variable)

def get_equivalent_deviatoric_strain(model):
	output_variable = 0.0
	for element in model.model_part.Elements:
		output_variable_list = element.GetValuesOnIntegrationPoints( EQUIVALENT_DEVIATORIC_STRAIN, model.model_part.ProcessInfo )
		output_variable = output_variable_list[0][0]
	return float(output_variable)


def get_strain_xx(model):
	""" compute the strain_xx at the gauss points """
	vector_size = 6
	strain_xx = 0.0
	# Get the STRAIN computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the STRAIN on each integration points:
		strains = element.GetValuesOnIntegrationPoints( STRAIN, model.model_part.ProcessInfo )
		for item in strains:
			strain_xx += item[0] 
	return float(strain_xx/len(strains))

def get_strain_yy(model):
	""" compute the strain_yy at the gauss points """
	vector_size = 6
	strain_yy = 0.0
	# Get the STRAIN computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the STRAIN on each integration points:
		strains = element.GetValuesOnIntegrationPoints( STRAIN, model.model_part.ProcessInfo )
		for item in strains:
			strain_yy += item[1] 
	return float(strain_yy/len(strains))

def get_strain_zz(model):
	""" compute the strain_zz at the gauss points """
	vector_size = 6
	strain_zz = 0.0
	# Get the STRAIN computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the STRAIN on each integration points:
		strains = element.GetValuesOnIntegrationPoints( STRAIN, model.model_part.ProcessInfo )
		for item in strains:
			strain_zz += item[2] 
	return float(strain_zz/len(strains))

def get_stress_xx(model):
	""" compute the stress_xx at the gauss points """
	vector_size = 6
	stress_xx = 0.0
	# Get the STRESSES computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the STRESSES on each integration points:
		stresses = element.GetValuesOnIntegrationPoints( STRESSES, model.model_part.ProcessInfo )
		for item in stresses:
			stress_xx += item[0] 
	return float(stress_xx/len(stresses))

def get_stress_yy(model):
	""" compute the stress_yy at the gauss points """
	vector_size = 6
	stress_yy = 0.0
	# Get the STRESSES computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the STRESSES on each integration points:
		stresses = element.GetValuesOnIntegrationPoints( STRESSES, model.model_part.ProcessInfo )
		for item in stresses:
			stress_yy += item[1] 
	return float(stress_yy/len(stresses))

def get_stress_zz(model):
	""" compute the stress_zz at the gauss points """
	vector_size = 6
	stress_zz = 0.0
	# Get the STRESSES computed on each integration points of the elements:
	for element in model.model_part.Elements:
		# Get the STRESSES on each integration points:
		stresses = element.GetValuesOnIntegrationPoints( STRESSES, model.model_part.ProcessInfo )
		for item in stresses:
			stress_zz += item[2] 
	return float(stress_zz/len(stresses))


def create_csv_file(file_name,data_x,data_y):
	with open(file_name, 'w') as f:
		# create the csv file
		for i in range(0,len(data_x)):
			f.write(str(data_x[i])+',' +str(data_y[i]))
			f.write('\n')