import sys
import os
import math
##################################################################
##################################################################
simulator_path = os.environ['BENCHMARK_PRISMA'] + "/material_point_application/Bouncing_balls"
sys.path.append(simulator_path)
import system_include
from system_include import *

mesh_division = [24, 24]
model = system_include.Model('mesh_' + str(mesh_division[0]) + "x" + str(mesh_division[1]), os.getcwd()+"/", os.getcwd()+"/", mesh_division)
model.material = "plane_strain"
model.InitializeModel()
# model.WriteOutput(0.0)

##################################################################
###  SIMULATION  #################################################
##################################################################

sys.path.append(os.environ['BENCHMARK_PRISMA'] + "/finite_cell_application")
import material_properties_utility

import simulator_multi

matfile_path = simulator_path + "/matfile.dat"

def GetParams():
    params = {}
    params["quadrature_method"] = "quadtree"
    params["number_of_samplings"] = 4
    params["qt_depth"] = 0
    params["cut_cell_quadrature_method"] = 0x03
    params["small_weight"] = 0.0 #1.0e-9
    params["action_on_small_cut_cell"] = "eliminate"
    params["number_of_minimum_physical_points"] = 1
    params["material_properties_utility"] = material_properties_utility.MaterialPropertiesUtility(matfile_path)
    params["material_properties_utility"].search_type = "by_name"
    params["material_properties_utility"].mat_type = "steel"
    # params["mat_type"] = "soil"
    params["mesh_division_x"] = mesh_division[0]
    params["mesh_division_y"] = mesh_division[1]
    params["write_output_per_each_step"] = False
    params["export_physical_integration_point"] = False
    params["physical_integration_point_prop_id"] = 11
    #
    params["number_of_balls"] = 2
    params["friction_coefficient"] = 0.0
    params["hit_angle"] = 15.0
    params["particle_definition"] = "solid_particle_large_strain"
    params["write_cbgeo"] = False
    params["delta_time"] = 1e-3
    params["total_time"] = 5.0e1
    params["sample_output"] = 100
    params["write_output"] = True
    params["write_energy"] = True
    params["velocity_magnitude"] = 0.1*math.sqrt(2.0)
    params["mpm_scheme"] = "USF"
    params["search_tolerance"] = 1e-10
    return params

def main(output=True, logging=True, total_time=50.0):

    params = GetParams()
    params["write_output"] = output
    params["write_energy"] = logging
    params["total_time"] = total_time

    sim = simulator_multi.Simulator(params)
    # sim.Initialize(model)
    return sim.Run(model)

def test():

    particle_manager_list, energies = main(output=False, logging=False, total_time=7.0)

    monitoring_particle_id = 11
    cx = particle_manager_list[0].Particles[monitoring_particle_id].X
    cy = particle_manager_list[0].Particles[monitoring_particle_id].Y
    energy = energies[0] + energies[1]

    ref_cx = 5.0146713487104866e-01
    ref_cy = 1.8263372235848863e-01
    ref_energy = 3.6642091963003378e+00

    print("cx: %.16e, cy: %.16e" % (cx, cy))
    print("energy: %.16e" % (energy))

    assert(abs(cx - ref_cx) < 1e-10)
    assert(abs(cy - ref_cy) < 1e-10)
    assert(abs(energy - ref_energy) < 1e-10)

    print("Test passed")

def tag():
    return "mpm"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(output=True)
