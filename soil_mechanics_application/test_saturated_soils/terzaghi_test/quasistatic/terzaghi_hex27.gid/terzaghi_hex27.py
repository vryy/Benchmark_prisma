##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./terzaghi_hex27.gid')
import terzaghi_hex27_include
from terzaghi_hex27_include import *

def main(output=True, logging=True):
    model = terzaghi_hex27_include.Model('terzaghi_hex27',os.getcwd()+"/",logging)
    model.InitializeModel()

    for e in model.layer_sets['Layer0']:
        elem = model.model_part.Elements[e]
        elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
        elem.SetValue(DENSITY_WATER,            1000.0)
        elem.SetValue(DENSITY_AIR,              0.0)
        elem.SetValue(K0,                       1.0)
        elem.SetValue(POROSITY,                 0.2)
        elem.SetValue(PERMEABILITY_WATER,       1.0e-6)
        elem.SetValue(PERMEABILITY_AIR,         0.0)
        elem.SetValue(FIRST_SATURATION_PARAM,   2.0)
        elem.SetValue(SECOND_SATURATION_PARAM,  0.0)
        elem.SetValue(AIR_ENTRY_VALUE,          1.0)

    tol = 1.0e-6
    top_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        if abs(node.Y0 - 0.0) < tol or abs(node.Y0 - 1.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        if abs(node.Z0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)
        if abs(node.Z0 - 7.0) < tol:
            top_nodes.append(node)

    pressure = -1.0e3
    for node in top_nodes:
        node.Fix(WATER_PRESSURE)
        node.SetSolutionStepValue(WATER_PRESSURE, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, 0.0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, 0.0)
        node.SetSolutionStepValue(FACE_LOAD_Z, pressure)

    time = 0.001
    model.Solve(time, 0, 0, 0, 0)
    if output:
        model.WriteOutput(time)

    delta_time = 0.001
    for step in range(0,10):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 0.01
    for step in range(0,10):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 0.1
    for step in range(0,10):
        time = time + delta_time
        model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 0.5
    for step in range(0,5):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 1.0
    for step in range(0,5):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 5.0
    for step in range(0,5):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 10.0
    for step in range(0,5):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 50.0
    for step in range(0,5):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 100.0
    for step in range(0,5):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    delta_time = 500.0
    for step in range(0,25):
        for inner_step in range(0,10):
            time = time + delta_time
            model.Solve( time, 0, 0, 0, 0 )
        if output:
            model.WriteOutput( time )
    print("Analysis completed")

    return model

def test():
    model = main(logging=False, output=False)

    ### pytesting results
    tol = 1.0e-6
    top_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.Z0 - 7.0) < tol:
            top_nodes.append(node)
    disp_z = top_nodes[0].GetSolutionStepValue(DISPLACEMENT_Z)
    ref_disp_z = -6.4968886835e-04
    print("%.10e" % disp_z)
    assert(abs((disp_z - ref_disp_z) / ref_disp_z) < 1e-10)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True)
