##################################################################
import sys
import os
##################################################################

def main(logging=True, output=True, nsteps=94, parallel_type = 'shared'):
    if(parallel_type == 'shared'):
        import quad_60_20_structural_include
        model = quad_60_20_structural_include.Model('quad_60_20',os.getcwd()+"/",os.getcwd()+"/", logging=logging)
    elif(parallel_type == 'distributed'):
        import quad_60_20_structural_parallel_include
        model = quad_60_20_structural_parallel_include.Model('quad_60_20',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    ##################################################################
    ###  SIMULATION  #################################################
    ##################################################################

    topopt_proc = model.solver.solver.attached_processes[0]

    time = 0.0
    delta_time = 1.0

    for step in range(0, nsteps):
        time = time + delta_time
        model.Solve(time, 0, 0, 0, 0)
        if output:
            model.WriteOutput(time)
        print("***********************************************************************")
        print("*****Step %d, Objective: %.16e, Density change: %.16e" % (step+1, topopt_proc.GetObjective(), topopt_proc.GetTopologyChange()))
        print("***********************************************************************")

    return model

def test():
    model = main(logging=False, output=False)

    topopt_proc = model.solver.solver.attached_processes[0]

    obj = topopt_proc.GetObjective()
    ref_obj= 2.0319245466199612e+02

    assert(abs(obj - ref_obj) < 1e-10)

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
