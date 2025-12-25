import sys

import simulation_script_block_excavation
from simulation_script_block_excavation import *

def main(logging=True, output=True, number_of_excavation_layers_per_step=8, sub_steps_range=[1, 1000]):
    params = {}
    params['path'] = os.getcwd() + '/'
    params['name'] = 'block_excavation'
    params['matfile_virgin'] = 'matfile.dat.plane_strain'
    params['mattype_virgin'] = 'elastoplastic'
    params['matfile'] = 'matfile.dat.mcc_ii_s'
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

def test_with_params(logging=False, output=False, number_of_excavation_layers_per_step=8, sub_steps_range=[1, 100], ref_u=[0.0, 0.0], tol=1e-10):
    model1 = main(logging=logging, output=output, number_of_excavation_layers_per_step=number_of_excavation_layers_per_step, sub_steps_range=sub_steps_range)

    tolp = 1e-6
    for node in model1.model_part.Nodes:
        if abs(node.X0 - 2.0) < tolp and abs(node.Y0 - 4.0) < tolp:
            pointA = node

    ux = pointA.GetSolutionStepValue(DISPLACEMENT_X)
    uy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    print(f"u ({number_of_excavation_layers_per_step}, {sub_steps_range}): %.16e, %.16e" % (ux, uy))
    assert(abs(ux - ref_u[0]) < tol)
    assert(abs(uy - ref_u[1]) < tol)

def test1(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=8, sub_steps_range=[1, 1], ref_u=[-2.9051121934968786e-02, -1.6145267312155581e-02])

def test2(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=8, sub_steps_range=[4, 4], ref_u=[-2.9447891242997924e-02, -1.6755546171158212e-02])

def test3(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=4, sub_steps_range=[1, 1], ref_u=[-2.8551663316528295e-02, -1.6267471607158496e-02])

def test4(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=4, sub_steps_range=[4, 4], ref_u=[-2.8670078763527947e-02, -1.6473494798127628e-02])

def test5(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=2, sub_steps_range=[1, 1], ref_u=[-2.8032541619871839e-02, -1.6149080652976346e-02])

def test6(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=2, sub_steps_range=[4, 4], ref_u=[-2.8180148638021225e-02, -1.6365879577463254e-02])

def test7(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=1, sub_steps_range=[1, 1], ref_u=[-2.7840609238153652e-02, -1.6020760319816710e-02])

def test8(logging=False, output=False):
    test_with_params(logging=logging, output=output, number_of_excavation_layers_per_step=1, sub_steps_range=[4, 4], ref_u=[-2.8016577117629032e-02, -1.6252762317854646e-02])

def test():
    # test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=1, ref_uy=4.6490466391350807e-02)
    # test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=4, ref_uy=4.7326001890955212e-02)
    # test_with_params(number_of_excavation_layers_per_step=4, sub_steps_range=1, ref_uy=4.2529271155825228e-02)
    # test_with_params(number_of_excavation_layers_per_step=4, sub_steps_range=4, ref_uy=4.2509708632310607e-02)
    # test_with_params(number_of_excavation_layers_per_step=2, sub_steps_range=1, ref_uy=4.1671666839188409e-02)
    # test_with_params(number_of_excavation_layers_per_step=2, sub_steps_range=4, ref_uy=4.1637795774181662e-02)
    # test_with_params(number_of_excavation_layers_per_step=1, sub_steps_range=1, ref_uy=4.1346788393986748e-02)
    # test_with_params(number_of_excavation_layers_per_step=1, sub_steps_range=4, ref_uy=4.1302839138054506e-02)

    # with HANDLE_ZERO_INCREMENT_STATE, point (0, 2)
    # test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=[1, 1], ref_uy=4.6490446268247762e-02)
    # test_with_params(number_of_excavation_layers_per_step=8, sub_steps_range=[4, 4], ref_uy=4.7326001888868194e-02)
    # test_with_params(number_of_excavation_layers_per_step=4, sub_steps_range=[1, 1], ref_uy=4.2529271155821211e-02)
    # test_with_params(number_of_excavation_layers_per_step=4, sub_steps_range=[4, 4], ref_uy=4.2509708632312571e-02)
    # test_with_params(number_of_excavation_layers_per_step=2, sub_steps_range=[1, 1], ref_uy=4.1671666839208150e-02)
    # test_with_params(number_of_excavation_layers_per_step=2, sub_steps_range=[4, 4], ref_uy=4.1637795774179608e-02)
    # test_with_params(number_of_excavation_layers_per_step=1, sub_steps_range=[1, 1], ref_uy=4.1346788394070931e-02)
    # test_with_params(number_of_excavation_layers_per_step=1, sub_steps_range=[4, 4], ref_uy=4.1302839138019784e-02)

    # with HANDLE_ZERO_INCREMENT_STATE, point (2, 4)
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    print("Test passed")

def tag():
    return "MCC"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        # main(logging=True, output=False, number_of_excavation_layers_per_step=4, sub_steps_range=[1, 1])
        test7(logging=True, output=False)
