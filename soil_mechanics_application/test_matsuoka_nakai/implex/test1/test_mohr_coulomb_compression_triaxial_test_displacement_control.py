import sys
import os
import math
import time as time_module
#importing Kratos main library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *
#from KratosMultiphysics.BRepApplication import *
#from KratosMultiphysics.EkateAuxiliaryApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
#from KratosMultiphysics.FreezingSoilApplication import *
kernel = Kernel()   #defining kernel
##################################################################
# TODO: ADD REFERENCE DIANA TEST PLOT RESULTS AUTOMATICALLY SAMPLE DATA FROM DIANA AN PLAXIS INSIDE THE SCRIPT

##################################################################
import element_driver_utility
from element_driver_utility import *
import element_driver_utility as edu

# import plotting_utility as pltu

path = os.getcwd()
name = "element_driver_model"

sys.path.append(path+"/"+name+'.gid')
print(sys.path)
print("...READING SYSTEM...")
system_include = __import__(name+"_include")

##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################
def main(output=True, logging=True, dilatancy_angle=15.0):

    model1 = system_include.Model(name,path+"/"+name+".gid/",path+"/"+name+".gid/",logging=logging,dilatancy_angle=dilatancy_angle)
    model1.InitializeModel()

    #------------------------------------------------------------------------------------------------
    #-> Parameters:
    #------------------------------------------------------------------------------------------------
    test_parameters = {}
    test_parameters['tolerance']    = 1E-6
    test_parameters['L']            = 1.0
    #------------------------------------------------------------------------------------------------
    #-> Displacement boundary condition are set to free:
    edu.free_displacement_bc(model1)
    #------------------------------------------------------------------------------------------------
    print("~~~~~~~~ Assign nodes sets/layers ~~~~~~~~")
    edu.assing_layers(model1,test_parameters['tolerance'],test_parameters['L'] )
    print("~~~~~~~~ Assign displacement boundary condition ~~~~~~~~")
    edu.assing_fix_displacement_bc(model1)

    print ("++++++++++ << Displacement boundary condition assigned: >> ++++++++++")

    ###################################################################################################
    print ("++++++++++ << Start triaxial test: >> ++++++++++")
    ###################################################################################################
    data_e_axial    = []
    data_q_dev      = []
    data_e_vol      = []
    #---------------------------------------------
    # loading: isotropic
    #---------------------------------------------
    # input:
    step_initial    = 0
    step_final      = 1
    time            = 0
    delta_time      = 1
    p_initial       = 1e6
    p_final         = 1e6
    delta_p         = (p_final -p_initial)/step_final
    p               = p_initial
    initial_prestress = {}
    initial_prestress["Sxx"] = 1e6
    initial_prestress["Sxy"] = 0.0
    initial_prestress["Sxz"] = 0.0
    initial_prestress["Syy"] = 1e6
    initial_prestress["Syz"] = 0.0
    initial_prestress["Szz"] = 1e6

    #edu.make_implex_implicit(model1)
    #---------------------------------------------------------------
    # we solve the corresponding displacements for isotropic pressure
    #---------------------------------------------------------------
    for step in range(step_initial,step_final):
        edu.assign_initial_prestress(model1,initial_prestress)
        edu.assign_isotropic_pressure(model1,p)
        # solve for the current time step:
        time = time + delta_time
        print("----------------------------------------------------------------------------------------------------------")
        print("Solving for the given pressure: " + str(p)+ " at time = " +str(time)+ " s")
        print("----------------------------------------------------------------------------------------------------------")

        model1.Solve(time,0,0,0,0)
        if output:
            model1.WriteOutput(time)
        # update:
        p += delta_p
        # prepare plot output
        e_axial = -edu.get_strain_zz(model1)
        q_dev = -(edu.get_stress_zz(model1) - edu.get_stress_xx(model1))
        e_vol = -edu.get_equivalent_volumetric_strain(model1)
        data_e_axial.append(e_axial)
        data_q_dev.append(q_dev)
        data_e_vol.append(e_vol)
    #---------------------------------------------
    # loading: triaxial
    #---------------------------------------------
    # input:
    step_initial    = 1
    step_final      = 100
    time            = 0.0
    delta_time      = 1.0
    u_z_initial     = 0.0
    u_z_final       = -4.0e-2
    delta_u_z       = (u_z_final - u_z_initial)/(step_final-step_initial)
    u_z             = u_z_initial + delta_u_z
    #---------------------------------------------------------------
    # we solve the corresponding displacements for the triaxial test
    #---------------------------------------------------------------
    for step in range(step_initial,step_final):

        edu.assing_fix_displacement_on_surface_z(model1,u_z)
        # solve for the current time step:
        time = time + delta_time
        print("----------------------------------------------------------------------------------------------------------")
        print("Solving for the given u_z: " + str(u_z)+ " at time = " +str(time)+ " s")
        print("----------------------------------------------------------------------------------------------------------")

        model1.Solve(time,0,0,0,0)
        if output:
            model1.WriteOutput(time)
        # update:
        u_z  += delta_u_z
        # prepare plot output
        e_axial = -edu.get_strain_zz(model1)
        q_dev = edu.get_equivalent_deviatoric_stress(model1) #edu.get_stress_zz(model1) - edu.get_stress_xx(model1)
        e_vol = -edu.get_equivalent_volumetric_strain(model1)
        data_e_axial.append(e_axial)
        data_q_dev.append(q_dev)
        data_e_vol.append(e_vol)

    if logging:
        #pltu.plot_strain_axial_vs_deviatoric_stress(data_e_axial,data_q_dev)
        #pltu.plot_volumetric_strain_vs_deviatoric_stain(data_e_axial,data_e_vol)
        #---------------------------------------------------------------
        # create csv file
        #---------------------------------------------------------------
        analysis_type = "triaxial_test_axial_strain_deviatoric_stress"
        file_name = analysis_type + str('.csv')
        edu.create_csv_file(file_name,data_e_axial,data_q_dev)
        analysis_type = "triaxial_test_axial_strain_volumetric_strain"
        file_name = analysis_type + str('.csv')
        edu.create_csv_file(file_name,data_e_axial,data_e_vol)

    return data_e_axial, data_q_dev, data_e_vol

def test():
    data_e_axial, data_q_dev, data_e_vol = main(logging=False, output=False, dilatancy_angle=0.0)
    print("%.16e" % (data_e_vol[-1]))
    assert(abs(data_e_vol[-1] - (-6.6666799999074502e-03)) < 1e-10)

    data_e_axial, data_q_dev, data_e_vol = main(logging=False, output=False, dilatancy_angle=15.0)
    print("%.16e" % (data_e_vol[-1]))
    assert(abs(data_e_vol[-1] - (7.3012474485577641e-03)) < 1e-10)

    data_e_axial, data_q_dev, data_e_vol = main(logging=False, output=False, dilatancy_angle=30.0)
    print("%.16e" % (data_e_vol[-1]))
    assert(abs(data_e_vol[-1] - (3.3333320000447010e-02)) < 1e-10)

    print("Test passed")

def tag():
    return "matsuoka-nakai"

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
