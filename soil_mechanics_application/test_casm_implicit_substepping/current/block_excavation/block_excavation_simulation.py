import sys

import simulation_script_block_excavation
from simulation_script_block_excavation import *

def main(logging=True, output=True, number_of_excavation_layers_per_step=8, number_of_sub_steps=1, viscous_damping=0.0):
    params = {}
    params['path'] = os.getcwd() + '/'
    params['name'] = 'block_excavation'
    params['matfile_virgin'] = 'matfile.dat.plane_strain'
    params['mattype_virgin'] = 'elastoplastic'
    params['matfile'] = 'matfile.dat.casm_s'
    params['soil_type'] = 'critical state'
    params['number_of_loading_steps_virgin'] = 1
    params['OCR_top'] = 1.7
    params['OCR_bottom'] = 1.0
    params['x_min'] = 0.0
    params['x_max'] = 5.0
    params['y_min'] = 0.0
    params['y_max'] = 4.0
    params['ground_water_table'] = 16.014
    params['number_of_steps'] = int(8 / number_of_excavation_layers_per_step)
    params['number_of_excavation_layers_per_step'] = number_of_excavation_layers_per_step
    params['number_of_sub_steps'] = number_of_sub_steps
    params['viscous_damping'] = viscous_damping
    params['time_excavation'] = 180.0
    params['transfer_method'] = "identical"
    params['account_for_water'] = False
    params['dry_run'] = False
    params['output'] = output
    params['logging'] = logging

    timer = Timer()
    timer.Start("Total execution time")
    simulator = BlockExcavationSimulator(params)

    model_virgin = simulator.PrepareVirgin(params)
    model1 = simulator.PrepareSystem(params)
    simulator.RunExcavationAnalysis(model_virgin, model1, params)

    return model1

def test_with_params(number_of_excavation_layers_per_step=8, number_of_sub_steps=1, ref_uy=0.0, tol=1e-10):
    model1 = main(logging=False, output=False, number_of_excavation_layers_per_step=number_of_excavation_layers_per_step, number_of_sub_steps=number_of_sub_steps)

    tolp = 1e-6
    for node in model1.model_part.Nodes:
        if abs(node.X0) < tolp and abs(node.Y0 - 2.0) < tolp:
            mon_node = node

    uy = mon_node.GetSolutionStepValue(DISPLACEMENT_Y)
    print(f"uy ({number_of_excavation_layers_per_step}, {number_of_sub_steps}): %.16e" % (uy))
    assert(abs(uy - ref_uy) < tol)

def test():
    test_with_params(number_of_excavation_layers_per_step=8, number_of_sub_steps=1, ref_uy=4.6147390906348673e-02)
    test_with_params(number_of_excavation_layers_per_step=8, number_of_sub_steps=4, ref_uy=4.7303426581162156e-02)
    test_with_params(number_of_excavation_layers_per_step=4, number_of_sub_steps=1, ref_uy=4.2410766155600098e-02)
    test_with_params(number_of_excavation_layers_per_step=4, number_of_sub_steps=4, ref_uy=4.2397059269436710e-02)
    # test_with_params(number_of_excavation_layers_per_step=2, number_of_sub_steps=1, ref_uy=4.1671666839188409e-02)
    # test_with_params(number_of_excavation_layers_per_step=2, number_of_sub_steps=4, ref_uy=4.1637795774181662e-02)
    test_with_params(number_of_excavation_layers_per_step=1, number_of_sub_steps=1, ref_uy=4.1215978910920575e-02)
    test_with_params(number_of_excavation_layers_per_step=1, number_of_sub_steps=4, ref_uy=4.1149122216637565e-02)

    print("Test passed")

def tag():
    return "CASM"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        # main(logging=True, output=True)
        main(logging=True, output=True, number_of_excavation_layers_per_step=2, number_of_sub_steps=1, viscous_damping=1e-3)
