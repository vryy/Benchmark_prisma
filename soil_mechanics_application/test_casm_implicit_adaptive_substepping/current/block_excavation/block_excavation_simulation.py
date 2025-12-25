import sys

import simulation_script_block_excavation
from simulation_script_block_excavation import *

def main(logging=True, output=True, number_of_excavation_layers_per_step=8, viscous_damping=0.0, abs_tol=1e-13, rel_tol=1e-13):
    params = {}
    params['path'] = os.getcwd() + '/'
    params['name'] = 'block_excavation'
    params['matfile_virgin'] = 'matfile.dat.plane_strain'
    params['mattype_virgin'] = 'elastoplastic'
    params['matfile'] = 'matfile.dat.casm_as'
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
    params['sub_steps_range'] = 20 # 10 # this parameter is number of levels in adaptive sub-stepping
    params['viscous_damping'] = viscous_damping
    params['time_excavation'] = 180.0
    params['transfer_method'] = "identical"
    params['account_for_water'] = False
    params['abs_tol'] = abs_tol
    params['rel_tol'] = rel_tol
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

def test_with_params(logging=False, output=False, number_of_excavation_layers_per_step=8, ref_u=[0.0, 0.0], tol=1e-10, viscous_damping=0.0, abs_tol=1e-13, rel_tol=1e-13):
    model1 = main(logging=logging, output=output, number_of_excavation_layers_per_step=number_of_excavation_layers_per_step, viscous_damping=viscous_damping, abs_tol=abs_tol, rel_tol=rel_tol)

    tolp = 1e-6
    for node in model1.model_part.Nodes:
        if abs(node.X0 - 2.0) < tolp and abs(node.Y0 - 4.0) < tolp:
            pointA = node

    ux = pointA.GetSolutionStepValue(DISPLACEMENT_X)
    uy = pointA.GetSolutionStepValue(DISPLACEMENT_Y)
    print(f"u ({number_of_excavation_layers_per_step}): %.16e, %.16e" % (ux, uy))
    assert(abs(ux - ref_u[0]) < tol)
    assert(abs(uy - ref_u[1]) < tol)

def test1(logging=False, output=False, abs_tol=1e-13, rel_tol=1e-13):
    test_with_params(logging=logging, output=output, \
        number_of_excavation_layers_per_step=8, \
        abs_tol=abs_tol, rel_tol=rel_tol, \
        ref_u=[-3.0037319193767124e-02, -1.5219816565195023e-02])

def test2(logging=False, output=False, abs_tol=1e-13, rel_tol=1e-13):
    test_with_params(logging=logging, output=output, \
        number_of_excavation_layers_per_step=4, \
        abs_tol=abs_tol, rel_tol=rel_tol, \
        ref_u=[-2.9478775925052682e-02, -1.5249802029518346e-02])

def test3(logging=False, output=False, abs_tol=1e-13, rel_tol=1e-13):
    test_with_params(logging=logging, output=output, \
        number_of_excavation_layers_per_step=2, \
        abs_tol=abs_tol, rel_tol=rel_tol, \
        viscous_damping=1e-5, \
        ref_u=[-2.8797079690590258e-02, -1.4918505615879419e-02])

def test4(logging=False, output=False, abs_tol=1e-13, rel_tol=1e-13):
    test_with_params(logging=logging, output=output, \
        number_of_excavation_layers_per_step=1, \
        abs_tol=abs_tol, rel_tol=rel_tol, \
        ref_u=[-2.8595405296473569e-02, -1.4785777999613420e-02])

def test():
    test1()
    test2()
    test3()
    test4()
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
        # main(logging=True, output=True, number_of_excavation_layers_per_step=2, viscous_damping=1e-5)
        # test4(logging=True, output=True, abs_tol=1e-13, rel_tol=1e-13)

        test_with_params(logging=True, output=True, \
            number_of_excavation_layers_per_step=2, \
            abs_tol=1e-13, rel_tol=1e-13, \
            viscous_damping=1e-5, \
            ref_u=[-2.8797079690590258e-02, -1.4918505615879419e-02])
