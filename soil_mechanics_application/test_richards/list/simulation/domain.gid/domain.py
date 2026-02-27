##################################################################
## Performing the example 2 as in the reference paper
## Reference:
##  +   List et al, A study on iterative methods for solving Richards’ equation
##  +   https://github.com/PMFlow/RichardsEquation.git
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + "/"
import domain_include as simulation_include
try:
    from domain_include import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################
model_name_ = 'domain'
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

aux_util = SoilsAuxiliaryUtility()
vtu = VariableTransferUtility(SuperLUSolver())

#silt loam parameters
Smax = 0.396
Smin = 0.131
n = 2.06
m = 1 - 1/n
pr = 1.0/0.423

def SetMaterialProperties(elem):
    elem.SetValue(USE_DISTRIBUTED_PROPERTIES, 1)
    elem.SetValue(PERMEABILITY_WATER,       4.96e-2)
    aux_util.SetValue(SWCC_LAW, VanGenuchtenSWCC(n, m, pr, Smin, Smax), elem)
    aux_util.SetValue(RELATIVE_PERMEABILITY_WATER_LAW, MualemRelativePermeabilityWaterHeadLaw(n, m, pr), elem)

def ReadMatlabVector(filename):
    ifile = open(filename, 'r')

    values = []
    for line in ifile.readlines():
        if line.startswith('#'):
            continue

        words = line.split()

        vec = []
        for w in words:
            vec.append(float(w))

        values.append(vec)

    ifile.close()

    return values

def main(logging=True, output=True):

    model = simulation_include.Model(model_name_,current_dir_,current_dir_,logging=logging, \
        solution_strategy="implicit_Newton_Raphson", \
        analysis_type=4, \
        # analysis_type=2, \
        # dissipation_radius=0.5, \
        # dissipation_radius=0.9, \
        dissipation_radius=1.0, \
        stop_Newton_Raphson_if_not_converge=True)
    model.InitializeModel()

    # material parameters
    for elem in model.model_part.Elements:
        SetMaterialProperties(elem)
        elem.Initialize(model.model_part.ProcessInfo)

    # initial condition
    for node in model.model_part.Nodes:
        h0 = 1.0 - node.Y0
        node.SetSolutionStepValue(WATER_PRESSURE, -h0)
        node.SetSolutionStepValue(WATER_PRESSURE_NULL, -h0)
        node.SetSolutionStepValue(WATER_PRESSURE_EINS, -h0)

    tol = 1e-6
    for node in model.model_part.Nodes: # on Gamma_D1
        if (node.X0 < 1.0 + tol) and abs(node.Y0 - 3.0) < tol:
            h0 = -2.0
            node.SetSolutionStepValue(WATER_PRESSURE_NULL, -h0)

    # boundary condition on Gamma_D2
    for node in model.model_part.Nodes:
        if abs(node.X0 - 2.0) < tol and (node.Y0 > -tol) and (node.Y0 < 1.0 + tol):
            h1 = 1.0 - node.Y0
            node.Fix(WATER_PRESSURE)
            node.SetSolutionStepValue(WATER_PRESSURE, -h1)
            node.SetSolutionStepValue(WATER_PRESSURE_EINS, -h1)
            node.SetSolutionStepValue(WATER_PRESSURE_NULL, -h1)
            # print("node %d WATER_PRESSURE is fixed" % (node.Id))

    T = 3.0/16
    dt = 1.0/48
    # dt = 1e-2
    t1 = T/3
    output_period = 1

    if output:
        model.WriteOutput(0.0)

    # time stepping
    model.model_part.ProcessInfo[FIRST_TIME_STEP] = 1

    t = 0.0
    i = 0
    while t <= T:
        t += dt
        print(f"### Time step {t}", flush=True)

        # time-dependent boundary condition
        if t <= t1:
            h0 = -2.0 + 2.2*t/t1
        else:
            h0 = 0.2
        for node in model.model_part.Nodes: # on Gamma_D1
            if (node.X0 < 1.0 + tol) and abs(node.Y0 - 3.0) < tol:
                node.Fix(WATER_PRESSURE)
                node.SetSolutionStepValue(WATER_PRESSURE_EINS, -h0)

        model.SolveModel(t)

        if output:
            # vtu.TransferVariablesToNodes(model.model_part, WATER_FLOW)
            if i%output_period == 0:
                model.WriteOutput(t)

        model.model_part.ProcessInfo[FIRST_TIME_STEP] = 0
        i += 1

        # if i == 1:
        #     break

    if output:
        ref_results = ReadMatlabVector(os.environ['HOME'] + '/workspace2/RichardsEquation/2D/Richards_2D/p.txt')
        # print(ref_results)

        error = 0.0
        for i in range(0, 21):
            for j in range(0, 31):
                x = (i / 20.) * 2
                y = (j / 30.) * 3

                found = False
                for node in model.model_part.Nodes:
                    if abs(node.X0 - x) < tol and abs(node.Y0 - y) < tol:
                        h = node.GetSolutionStepValue(WATER_PRESSURE)
                        hr = ref_results[j][i]
                        print('x: %e, y: %e, h: %.16e, ref: %e, diff: %e' % (x, y, h, hr, h+hr))
                        node.SetSolutionStepValue(TEMPERATURE, h + hr)
                        error += abs(h+hr)
                        found = True
                        break

                if not found:
                    raise Exception("Cannot find node at %e, %e" % (x, y))

        print("error: %e" % (error))

        model.WriteOutput(t)

    return model

def test():
    model = main(logging=False, output=False)

    tol = 1e-6
    for node in model.model_part.Nodes:
        if (abs(node.X0 - 1.0) < tol) and abs(node.Y0 - 2.5) < tol:
            mon_node = node

    wh = mon_node.GetSolutionStepValue(WATER_PRESSURE)
    print("wh: %.16e" % (wh))
    ref_wh = 8.3331156939562168e-01
    assert(abs(wh - ref_wh) < 1e-12)

    print("Test passed")

def tag():
    tags = "richards"
    if not all_modules_are_imported_successfully:
        tags += ",untested"
    return tags

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
