##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014-2022 by Hoang-Giang Bui (SFB837)         #####
#####          2023-2024 by Hoang-Giang Bui (Hereon)         #####
#####          2025 by Hoang-Giang Bui (UoB)                 #####
#####          2026 by Hoang-Giang Bui (DU)                  #####
##### all rights reserved                                    #####
##################################################################
# Reference:
# + geometry: Sepulveda-Cano et al, Numerical analysis of soil desaturation by an air injection method
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import simulation_include
try:
    from simulation_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'ground1_t6'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

aux_util = SoilsAuxiliaryUtility()
vtu = VariableTransferUtility(SuperLUSolver())

Smax = 0.99
Smin = 0.01
n = 8.696
m = 1 - 1/n
alpha = 2.237e-4
pr = 1.0 / alpha

theta = 0.421

rho_w = 1e3
rho_a = 1.28
mu_w = 1e-3
mu_a = 1.81e-5
kint = 2.04e-11
g = 9.81

kw = kint / mu_w * rho_w * g
ka = kint / mu_a * rho_a * g

g = 9.81

def SetMaterialProperties(elem):
    elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
    elem.SetValue(DENSITY,                  2000.0)
    elem.SetValue(DENSITY_WATER,            rho_w)
    elem.SetValue(POROSITY,                 theta)
    elem.SetValue(PERMEABILITY_WATER,       kw)
    elem.SetValue(PERMEABILITY_AIR,         ka)
    aux_util.SetValue(SWCC_LAW, InversedVanGenuchtenSWCC1(n, m, pr, Smin, Smax).Create(1e-3), elem)
    aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, MualemRelativePermeabilityWaterLaw(m, Smin, Smax), elem)
    aux_util.SetValue(RELATIVE_PERMEABILITY_AIR_LAW, MualemRelativePermeabilityAirLaw(m, Smin, Smax), elem)
    aux_util.SetValue(GAS_LAW, IdealGasLaw(rho_a, rho_a*1e-5), elem)
    elem.SetValue(FIX_POROSITY,             True)

def main(output=True, logging=True, total_time=100.0, ramp_time = 1e2, ramp_steps = 60, pinj=90e3, \
    solution_strategy="implicit_Newton_Raphson", analysis_type=1, dissipation_radius=0.1):

    ## solve the system

    model = simulation_include.Model(model_name_,os.getcwd()+"/",os.getcwd()+"/", \
        logging=logging, \
        solution_strategy=solution_strategy, \
        analysis_type=analysis_type, \
        dissipation_radius=dissipation_radius, \
        stop_Newton_Raphson_if_not_converge=True)
    model.InitializeModel()

    # material parameters
    for e in model.layer_sets['Ground']:
        elem = model.model_part.Elements[e]
        SetMaterialProperties(elem)
        elem.Initialize(model.model_part.ProcessInfo)

    # boundary condition
    tol = 1.0e-6
    density_water = rho_w
    y_max = 15.0
    top_nodes = []
    bottom_nodes = []
    lateral_nodes = []
    for node in model.model_part.Nodes:
        if abs(node.X0 - 0.0) < tol or abs(node.X0 - 30.0) < tol:
            lateral_nodes.append(node)
        if abs(node.Y0 - 0.0) < tol:
            bottom_nodes.append(node)
        if abs(node.Y0 - 15.0) < tol:
            top_nodes.append(node)

        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)

        node.Fix(SATURATION)
        node.Fix(AIR_PRESSURE)

        node.SetSolutionStepValue(SATURATION, Smax)
        node.SetSolutionStepValue(SATURATION_EINS, Smax)
        node.SetSolutionStepValue(SATURATION_NULL, Smax)

        # node.SetSolutionStepValue(AIR_PRESSURE, p)
        # node.SetSolutionStepValue(AIR_PRESSURE_EINS, p)
        # node.SetSolutionStepValue(AIR_PRESSURE_NULL, p)

    for element in model.model_part.Elements:
        water_pressures = element.GetValuesOnIntegrationPoints( WATER_PRESSURE, model.model_part.ProcessInfo )
        ipoints = element.GetValuesOnIntegrationPoints( INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model.model_part.ProcessInfo )
        pressure_list = []
        for item in ipoints:
            p = density_water*g*(y_max - item[1])
            pressure_list.append( p )
        element.SetValuesOnIntegrationPoints( REFERENCE_WATER_PRESSURE, pressure_list, model.model_part.ProcessInfo )

    ## reset the material one more time to account for new information
    for element in model.model_part.Elements:
        element.ResetConstitutiveLaw()

    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    # release the saturation and air pressure on the model
    for node in model.model_part.Nodes:
        node.Free(SATURATION)
        node.Free(AIR_PRESSURE)

    # but fix air pressure on top and lateral
    for node in top_nodes + lateral_nodes:
        node.Fix(AIR_PRESSURE)

    # # fix water pressure on top and lateral
    # TODO we need to add the constraint here
    # for node in top_nodes + lateral_nodes:
    #     node.Fix(WATER_PRESSURE)

    # for node in top_nodes:
    for node in top_nodes:
        node.Fix(SATURATION)

    # pressure ramp-up
    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    time = 0.0
    delta_time = ramp_time / ramp_steps
    delta_pres = pinj / ramp_steps
    for i in range(0, ramp_steps):

        p = (i+1)*delta_pres
        print(f"### Applying air pressure {p}", flush=True)

        # apply the injection pressure
        for node_id in model.node_groups['Injection']:
            node = model.model_part.Nodes[node_id]
            node.Fix(AIR_PRESSURE)
            node.SetSolutionStepValue(AIR_PRESSURE, p)
            node.SetSolutionStepValue(AIR_PRESSURE_NULL, p)
            node.SetSolutionStepValue(AIR_PRESSURE_EINS, p)

        time += delta_time
        model.SolveModel(time)

        # if output:
        #     vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)
        #     vtu.TransferVariablesToNodes(model.model_part, AIR_FLOW)

        if output:
            model.WriteOutput(time)

        model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0


    cnt = 1
    while True:
        time = time + delta_time

        if time > total_time:
            break

        print(f"### Air pressure propagating phase, time = {time}", flush=True)

        converged = model.SolveModel( time )

        # if output:
        #     vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)
        #     vtu.TransferVariablesToNodes(model.model_part, AIR_FLOW)

        if not converged:
            raise Exception(f"The time stepping does not converge at time step {time}")

        if output and (cnt%10 == 0):
            model.WriteOutput( time )

        cnt += 1

    return model

def test():
    model = main(logging=False, output=False, pinj=90e3, ramp_steps=40, total_time=130.0, analysis_type=4, dissipation_radius=1.0)

    mon_node = model.model_part.Nodes[299]

    pa = mon_node.GetSolutionStepValue(AIR_PRESSURE)
    ref_pa = 8.3321926417441719e+04
    print("pa: %.16e, diff: %.16e" % (pa, pa - ref_pa))
    assert(abs(pa - ref_pa) / abs(ref_pa) < 1e-10)

    print("Test passed")

def tag():
    if all_modules_are_imported_successfully:
        return "partially-saturated"
    else:
        return ""

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=True, pinj=90e3, ramp_steps=50, total_time=1000.0, analysis_type=4, dissipation_radius=1.)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
