##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020,     #####
#####     2021, 2022 by Hoang-Giang Bui for SFB837           #####
#####     2023 by Hoang-Giang Bui                            #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Sa 30. Sep 15:29:31 CEST 2023
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./two_cubes_32_hex8.gid')
import two_cubes_32_hex8_include
from two_cubes_32_hex8_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

import tsi_partitioning_strategy

def WriteAnalyticalSolution(n, max_contact_force, ac, k):

    ifile = open("analytical_solution.txt", "w")
    ifile.write("contact_force\ttemp-slave\ttemp-master\n")

    for i in range(0, n):
        pn = i*max_contact_force/(n-1)
        eta = ac*pn/k
        temp_slave = (40.0*eta+20.0*(1+eta)) / (1.0+2.0*eta)
        temp_master = (20.0*eta+40.0*(1+eta)) / (1.0+2.0*eta)
        ifile.write("%.10e\t%.10e\t%.10e\n" % (pn, temp_slave, temp_master))

    ifile.close()

def WriteLog(ifile, slave_nodes, master_nodes, disp, ac, k):
    contact_force_slave = 0.0
    for node in slave_nodes:
        contact_force_slave += node.GetSolutionStepValue(CONTACT_FORCE_Z)

    temp_slave = 0.0
    for node in slave_nodes:
        temp_slave += node.GetSolutionStepValue(TEMPERATURE)

    contact_force_master = 0.0
    for node in master_nodes:
        contact_force_master += node.GetSolutionStepValue(CONTACT_FORCE_Z)

    temp_master = 0.0
    for node in master_nodes:
        temp_master += node.GetSolutionStepValue(TEMPERATURE)

    print("contact force slave: %.10e, master: %.10e" % (contact_force_slave, contact_force_master))
    print("contact force slave (avg): %.10e, master (avg): %.10e" % (contact_force_slave / len(slave_nodes), contact_force_master / len(master_nodes)))

    pn = contact_force_slave / len(slave_nodes)

    ## analytical solution is taken from
    ## Wriggers et al, Contact constraints within coupled thermomechanical analysis - A finite element model, 1994
    eta = ac*pn/k
    temp_slave_ana = (40.0*eta+20.0*(1+eta)) / (1.0+2.0*eta)
    temp_master_ana = (20.0*eta+40.0*(1+eta)) / (1.0+2.0*eta)

    ifile.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (disp, pn, temp_slave / len(slave_nodes), temp_master / len(master_nodes), temp_slave_ana, temp_master_ana))

    return pn

### partitioning alrithm

### end of partitioning algorithm

def main(output=True, logging=True):

    model = two_cubes_32_hex8_include.Model('two_cubes_32_hex8',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    tol = 1.0e-6
    slave_nodes = []
    master_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
            node.Fix(TEMPERATURE)
            node.SetSolutionStepValue(TEMPERATURE, 40.0)
        if abs(node.Z0 - 1.0) < tol:
            master_nodes.append(node)
        if abs(node.Z0 - 1.1) < tol:
            slave_nodes.append(node)

    prescribed_nodes = []
    for node_id in model.layer_nodes_sets['load']:
        node = model.model_part.Nodes[node_id]
        node.Fix(DISPLACEMENT_Z)
        node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, 20.0)
        prescribed_nodes.append(node)

    ## contact parameters
    model.solver.solver.contact_tying_indices = {}
    # model.solver.solver.contact_tying_indices[10] = "contact_link_kinematic_linear_penalty"
    # model.solver.solver.contact_tying_indices[10] = "contact_link_kinematic_linear_penalty_no_linearized"
    model.solver.solver.contact_tying_indices[10] = Tm_Contact_Link_Kinematic_Linear_Penalty_No_Linearized_3D()
    model.solver.solver.Parameters['penalty'] = {10: 1.0e6} #{10: 1.0e7}
    model.solver.solver.Parameters['friction_coefficient'] = {10: 0.0}
    # model.solver.solver.Parameters['solution_strategy'] =
    model.solver.solver.Parameters['dimension'] = 3
    model.solver.solver.Parameters['gap_tolerance'] = 1.0e-10
    model.solver.solver.Parameters['tying_integration_order'] = 4
    model.solver.solver.Parameters['predict_local_point_method'] = 0
    model.solver.solver.Parameters['maximal_detection_distance'] = 1e-3
    model.solver.solver.Parameters['polygon_offset'] = 0.01
    model.solver.solver.Parameters['stop_active_set_if_not_converged'] = False
    model.solver.solver.Parameters['stop_Newton_Raphson_if_not_converged'] = False
    model.solver.solver.contact_tying_values = {}
    model.solver.solver.contact_tying_values[10] = {THERMAL_CONTACT_CONDUCTANCE: 10.0}
    #model.solver.solver.Parameters['test_linearization'] = False
    #model.solver.solver.Parameters['test_linearization_disp'] = 1.0e-7
    #model.solver.solver.Parameters['test_linearization_tol'] = 1.0e-6
    #model.solver.solver.Parameters['visualize_contact_pairs'] = False
    model.solver.solver.Parameters['compute_contact_force'] = True
    #model.solver.solver.Parameters['active_set_rel_tol'] = 1.0e-8
    #model.solver.solver.Parameters['active_set_abs_tol'] = 1.0e-12
    #model.solver.solver.Parameters['max_active_set_iter'] = 10
    model.solver.solver.InitializeContact()

    #################################################

    if logging:
        ifile = open("monitoring.log", "w")
        ifile.write("disp\tcontact_force\ttemp-slave\ttemp-master\ttemp-slave-ana\ttemp-master-ana\n")

    #################################################

    ac = 10.0
    k = model.model_part.Properties[1].GetValue(THERMAL_CONDUCTIVITY)
    solve_tolerance = 1e-6

    time = 0.0
    tsi_partitioning_strategy.Timeloop(model, time, tol=solve_tolerance)
    if output:
        model.WriteOutput(time)
    # sys.exit(0)

    disp = 0.0
    delta_disp = 0.1
    delta_time = 0.1

    time = time + delta_time
    disp = disp + delta_disp
    for node in prescribed_nodes:
        node.SetSolutionStepValue(DISPLACEMENT_Z, -disp)
    tsi_partitioning_strategy.Timeloop(model, time, tol=solve_tolerance)
    if output:
        model.WriteOutput(time)
    if logging:
        pn = WriteLog(ifile, slave_nodes, master_nodes, disp, ac, k)

    delta_disp = 1e-3
    delta_time = delta_disp
    for i in range(0, 15):
        print("############STARTING STEP", i, "###############")
        time = time + delta_time
        disp = disp + delta_disp
        for node in prescribed_nodes:
            node.SetSolutionStepValue(DISPLACEMENT_Z, -disp)
            # print(node.GetSolutionStepValue(DISPLACEMENT))
        tsi_partitioning_strategy.Timeloop(model, time, tol=solve_tolerance)
        if output:
            model.WriteOutput(time)
        if logging:
            pn = WriteLog(ifile, slave_nodes, master_nodes, disp, ac, k)

    if logging:
        WriteAnalyticalSolution(100, pn, ac, k)
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False)

    ###### pytesting results
    tol = 1.0e-6
    slave_nodes = []
    master_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Z0 - 1.0) < tol:
            master_nodes.append(node)
        if abs(node.Z0 - 1.1) < tol:
            slave_nodes.append(node)
    slave_temp = 0.0
    for node in slave_nodes:
        slave_temp += node.GetSolutionStepValue(TEMPERATURE)
    slave_temp /= len(slave_nodes)
    master_temp = 0.0
    for node in master_nodes:
        master_temp += node.GetSolutionStepValue(TEMPERATURE)
    master_temp /= len(master_nodes)
    print("slave_temp: ", slave_temp)
    print("master_temp: ", master_temp)
    ref_slave_temp = 29.19952063565615
    ref_master_temp = 30.79721851939952
    assert(abs(slave_temp - ref_slave_temp) / abs(ref_slave_temp) < 1e-10)
    assert(abs(master_temp - ref_master_temp) / abs(ref_master_temp) < 1e-10)
    ########################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
