import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import mesh21x21_include
from mesh21x21_include import *
##################################################################
###  SIMULATION  #################################################
##################################################################

import analytical_solution

import material_properties_utility

import simulator

#sys.path.append("./")
#import quadrature_data

def compute_fit():
    model = mesh21x21_include.Model('mesh21x21',os.getcwd()+"/",os.getcwd()+"/",logging=False)
    model.InitializeModel()

    params = {}
    params["qt_depth"] = 4
    params["number_of_samplings"] = 4
    params["cut_cell_quadrature_method"] = 0x03
    params["fitting_space_dimension"] = 2
    params["fitting_function_degree"] = 2
    params["fit_mode"] = "multithread"
    params["integrator_quadrature_method"] = 0x02
    params["fit_mode"] = "multithread" #"serial"
    params["fit_solver_type"] = "direct" #"nnls" #
    params["fit_echo_level"] = 3
    params["fit_small_weight"] = 0.0
    params["export_physical_integration_point"] = False
    params["physical_integration_point_prop_id"] = 88
    params["export_quadtree_cell"] = False
    params["sample_quad_element_name"] = "DummySurfaceElement2D4N"
    params["write_quadrature_to_file"] = True
    params["quad_filename"] = "quadrature_data_3x3_p2_d" + str(params["qt_depth"]) + ".py"
    # params["quad_filename"] = "quadrature_data_3x3_p3_d" + str(params["qt_depth"]) + ".py"
    #params["quad_filename"] = "quadrature_data_4x4_d" + str(params["qt_depth"]) + ".py"
    params["quad_filetype"] = "python"
    params["quad_accuracy"] = 20

    start = time_module.time()
    sim = simulator.Simulator(params)
    sim.MomentFit(model.model_part.Elements)
    end = time_module.time()
    print("Fitting time: " + str(end-start) + " s")

def main(logging=True, output=True):
    model = mesh21x21_include.Model('mesh21x21',os.getcwd()+"/",os.getcwd()+"/",logging=logging,post_mode="binary")
    model.InitializeModel()

    import quadrature_data_3x3_p2_d4 as quadrature_data

    params = {}
    params["analytical_solution"] = analytical_solution
    params["quadrature_method"] = "moment-fit quadtree"
    params["number_of_samplings"] = 4
    params["qt_depth"] = 4
    params["cut_cell_quadrature_method"] = 0x03
    params["quadrature_data"] = quadrature_data
    params["material_properties_utility"] = material_properties_utility.MaterialPropertiesUtility("matfile.dat")
    params["material_properties_utility"].search_type = "by_name"
    params["material_properties_utility"].mat_type = "elastic"
    params["export_physical_integration_point"] = True
    params["physical_integration_point_prop_id"] = 88
    params["export_quadtree_cell"] = False
    params["sample_quad_element_name"] = "DummySurfaceElement2D4N"
    params["write_quadrature_to_file"] = False
    params["quad_filename"] = "quadrature_data.py"
    params["quad_filetype"] = "python"
    params["quad_accuracy"] = 20
    params["enable_ghost_penalty"] = False
    params["space_dim"] = 2
    params["estimated_number_of_neighbours"] = 10
    params["sample_ghost_penalty_condition"] = GhostPenaltyDisplacementGradientCondition()
    params["ghost_penalty_properties"] = model.model_part.Properties[4]
    model.model_part.Properties[4].SetValue(GHOST_PENALTY_STABILIZATION_FACTOR, 1.0e-1)
    model.model_part.Properties[4].SetValue(GHOST_PENALTY_STABILIZATION_ORDER, 2)

    sim = simulator.Simulator(params)
    sim.Run(model)

    if output:

        #################### TRANSFER THE RESULT TO NODES ####################
        sim.variable_projection_utility.TransferVariablesToNodes(VON_MISES_STRESS, model.model_part.ProcessInfo)

        model.WriteOutput(1.0)

        sim.GetVMSOnBoundaryByInterpolation(model.model_part, 500, "vms.txt")
        sim.GetAnalyticalVMSOnBoundary(200, "vms_ana.txt")

        #################### TRANSFER THE RESULT TO OTHER MESH FOR VISUALIZATION ####################
        #model_output_path = os.getcwd() + "/../../../../infinite_plate_with_hole/post_meshes/mesh_005.gid"
        #sys.path.append(model_output_path)
        #import mesh_005_include
        #from mesh_005_include import *
        #model_output = mesh_005_include.Model('mesh_005', model_output_path+"/", output_path)
        #model_output.InitializeModel()

        sys.path.append(os.getcwd() + "/mesh_001.gid")
        import mesh_001_include
        model_output = mesh_001_include.Model('mesh_001', os.getcwd() + "/mesh_001.gid/", os.getcwd(), logging=False)
        model_output.InitializeModel()

        sim.variable_transfer_utility.TransferVariablesFromNodeToNode(model.model_part, model_output.model_part, VON_MISES_STRESS)
        sim.variable_transfer_utility.TransferVariablesFromNodeToNode(model.model_part, model_output.model_part, DISPLACEMENT)
        model_output.WriteOutput(1.0)

        sim.GetVMSOnBoundary(model_output.model_part, 1.0e-2, "vms_post.txt")

    if output or logging:
        sim.ComputeError(model)

    return sim, model

def test():
    sim, model = main(logging=False, output=False)

    strain_energy_exact, l2_error, h1_error = sim.ComputeError(model)

    ref_strain_energy = 8.182504326516675e-03
    ref_l2_error = 4.053315664006885e-03
    ref_h1_error = 3.423123473875946e-02

    assert(abs(strain_energy_exact - ref_strain_energy) < 1e-10)
    assert(abs(l2_error - ref_l2_error) < 1e-10)
    assert(abs(h1_error - ref_h1_error) < 1e-10)

    print("Test passed")

def tag():
    return "finite-cell"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
