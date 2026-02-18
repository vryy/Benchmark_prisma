import sys

import simulation_script_block_excavation
from simulation_script_block_excavation import *

def main(logging=True, output=True, number_of_excavation_layers_per_step=8, sub_steps_range=[1, 1], \
    viscous_damping=0.0, abs_tol=1e-13, rel_tol=1e-13, local_error_tolerance=1e-6):
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
    params['sub_steps_range'] = sub_steps_range
    params['viscous_damping'] = viscous_damping
    params['time_excavation'] = 180.0
    params['transfer_method'] = "identical"
    params['account_for_water'] = False
    params['abs_tol'] = abs_tol
    params['rel_tol'] = rel_tol
    params['local_error_tolerance'] = local_error_tolerance
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

def test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=[1, 1], ref_u=[0.0, 0.0], tol=1e-10, viscous_damping=0.0):
    model1 = main(logging=False, output=False, number_of_excavation_layers_per_step=number_of_excavation_layers_per_step, sub_steps_range=sub_steps_range, viscous_damping=viscous_damping)

    tolp = 1e-6
    for node in model1.model_part.Nodes:
        if abs(node.X0 - 2.0) < tolp and abs(node.Y0 - 4.0) < tolp:
            pointA = node

    ux = pointA.GetSolutionStepValue(DISPLACEMENT_X)
    uy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    print(f"u ({number_of_excavation_layers_per_step}, {sub_steps_range}): %.16e, %.16e" % (ux, uy))
    assert(abs(ux - ref_u[0]) < tol)
    assert(abs(uy - ref_u[1]) < tol)

def test():
    test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=[1, 16], ref_u=[-2.9219068214881776e-02, -1.4419627060506536e-02])
    test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=[4, 16], ref_u=[-2.9524305602192845e-02, -1.4800707481889360e-02])
    test_with_params(number_of_excavation_layers_per_step=4, sub_steps_range=[1, 4], ref_u=[-2.9016965365597010e-02, -1.4843600199118174e-02])
    test_with_params(number_of_excavation_layers_per_step=4, sub_steps_range=[4, 4], ref_u=[-2.9118654952549262e-02, -1.4977494796756670e-02])
    # test_with_params(number_of_excavation_layers_per_step=2, sub_steps_range=[1, 8], ref_u=[-2.9118654952549262e-02, -1.4977494796756670e-02], viscous_damping=1e-5)
    # test_with_params(number_of_excavation_layers_per_step=2, sub_steps_range=[4, 4], ref_u=[-2.9118654952549262e-02, -1.4977494796756670e-02])
    # test_with_params(number_of_excavation_layers_per_step=1, sub_steps_range=[1, 2], ref_u=[-2.8586930631622601e-02, -1.4774990705525131e-02])
    test_with_params(number_of_excavation_layers_per_step=1, sub_steps_range=[4, 4], ref_u=[-2.8414126360852200e-02, -1.4718068239814870e-02])

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
        main(logging=True, output=True, number_of_excavation_layers_per_step=2, sub_steps_range=[4, 8], viscous_damping=1e-5, abs_tol=1e-10, rel_tol=1e-10, local_error_tolerance=1e-8)
        # main(logging=True, output=True, number_of_excavation_layers_per_step=1, sub_steps_range=[4, 8], viscous_damping=0.0)
