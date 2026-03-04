##################################################################
import sys
import os
##################################################################
current_dir_ = os.path.dirname(os.path.realpath(__file__)) + os.sep
if current_dir_ not in sys.path:
    sys.path.append(current_dir_)
##################################################################
try:
    from KratosMultiphysics import *
    from KratosMultiphysics.StructuralApplication import *
    from KratosMultiphysics.TopOptApplication import *
    all_modules_are_imported_successfully = True
except Exception as e:
    all_modules_are_imported_successfully = False
##################################################################

def main(logging=True, output=True, nsteps=94, parallel_type = 'shared'):
    if(parallel_type == 'shared'):
        import MBB_60_20_20_structural_include
        model = MBB_60_20_20_structural_include.Model('MBB_60_20_20',current_dir_,current_dir_,logging=logging)
    elif(parallel_type == 'distributed'):
        import MBB_60_20_20_structural_parallel_include
        model = MBB_60_20_20_structural_parallel_include.Model('MBB_60_20_20',current_dir_,current_dir_,logging=logging)

    ##################################################################
    ## TOPOLOGY OPTIMIZATION SETUP ###################################
    ##################################################################
    volfrac =          0.5
    rmin =          1.5
    ft = 1
    model.topopt_proc = TopologyUpdateProcess(model.model_part, volfrac, rmin, ft)
    model.solver.solver.attached_processes.append(model.topopt_proc)

    ##################################################################
    ## ADDITIONAL MATERIAL PARAMETERS ################################
    ##################################################################
    model.model_part.Properties[1].SetValue(YOUNG_MODULUS_MIN,        1e-09 )
    model.model_part.Properties[1].SetValue(PENALIZATION_FACTOR,            3 )

    ##################################################################
    ## INITIALIZATION ################################################
    ##################################################################
    model.InitializeModel()

    ##################################################################
    ###  SIMULATION  #################################################
    ##################################################################

    time = 0.0
    delta_time = 1.0

    for step in range(0, nsteps):
        time = time + delta_time
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("***********************************************************************")
        print("*****Step %d, Objective: %.16e, Density change: %.16e" % (step+1, model.topopt_proc.GetObjective(), model.topopt_proc.GetTopologyChange()))
        print("***********************************************************************")

    return model

def test():
    model = main(logging=False, output=False, nsteps=2)

    topopt_proc = model.topopt_proc

    obj = topopt_proc.GetObjective()
    ref_obj = 1.1011948347047752e+07

    assert(abs(obj - ref_obj) / ref_obj < 1e-10)

    print("Test passed")

def tag():
    tags = "topopt"
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
