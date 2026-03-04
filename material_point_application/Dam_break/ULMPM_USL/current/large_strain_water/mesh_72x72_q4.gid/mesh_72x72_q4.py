##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import mesh_72x72_q4_include
from mesh_72x72_q4_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def main(logging=True, output=True, total_time=5.0):

    model = mesh_72x72_q4_include.Model('mesh_72x72_q4',os.getcwd()+"/",os.getcwd()+"/")
    model.InitializeModel()

    ## initial condition
    for node in model.model_part.Nodes:
        node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
        node.SetSolutionStepValue(VELOCITY_X, 0.0)
        node.SetSolutionStepValue(VELOCITY_Y, 0.0)
        node.SetSolutionStepValue(ACCELERATION_X, 0.0)
        node.SetSolutionStepValue(ACCELERATION_Y, 0.0)

    ## get the dam elements
    dam_elements = ElementsArray()
    for elem_id in model.layer_sets['Dam']:
        elem = model.model_part.Elements[elem_id]
        dam_elements.append(elem)

    ## create particle manager
    particle_manager = SolidParticleManager()
    # query_tool = SpatialGridElementalBinning(0.0, 0.0, 0.0, 1.0e-2, 1.0e-2, 0.0, 1.0e-10)
    query_tool = StructuredGridElementalIndexing2D(0.0, 0.0, 6.0/72, 6.0/72, 72, 72)
    # query_tool = StructuredGridElementalIndexing2D(0.0, 0.0, 6.0/72, 6.0/72, 84, 84)
    query_tool.Initialize(model.model_part.Elements)
    particle_manager.SetQueryTool(query_tool)
    particle_manager.EchoLevel = 1
    sample_particle = SolidParticleLargeStrain2D(0, 0.0, 0.0)
    lastId = 0
    lastId = particle_manager.AddParticles(dam_elements, sample_particle, lastId, model.model_part.ProcessInfo)

    for p in particle_manager.Particles:
        p.Vx = 0.0
        p.Vy = 0.0

    ## create the MP model
    mpm = ULMpmModel(model.model_part, particle_manager)
    mpm.Check()
    print(mpm)

    # boundary condition
    tol = 1.0e-6
    L = 6.0
    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) or (abs(node.X0-L) < tol) or (node.X0 > L):
            node.Fix(DISPLACEMENT_X)
            node.Fix(VELOCITY_X)
            node.Fix(ACCELERATION_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(VELOCITY_X, 0.0)
            node.SetSolutionStepValue(ACCELERATION_X, 0.0)

        if (abs(node.Y0) < tol) or (abs(node.Y0-L) < tol) or (node.Y0 > L):
            node.Fix(DISPLACEMENT_Y)
            node.Fix(VELOCITY_Y)
            node.Fix(ACCELERATION_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(ACCELERATION_Y, 0.0)

    ## create the solver and solve
    import mpm_explicit_strategy
    mpm_solver = mpm_explicit_strategy.MPMExplicitStructuralUSLSolver(mpm)
    mpm_solver.Initialize()

    time = 0.0
    hx = 1.0/72
    kappa = 2.0e6
    rho = 997.5
    delta_time = 0.1*hx/math.sqrt(kappa/rho)
    print("delta_time: " + str(delta_time))

    step = 0
    sample_output = 100
    sample_cond_name = "PointForce2D"
    if logging:
        itime = open("time.txt", "w")
        itime.write("step\ttime\n")

    while time < total_time:
        time = time + delta_time
        print("#########################################################")
        print("######### TIME STEP " + str(time) + " BEGIN #############")

        mpm_solver.Solve(time)

        if output:
            if step % sample_output == 0:
                [list_nodes, list_conds] = particle_manager.ExportParticles(model.model_part, sample_cond_name, 10000, 10000, model.model_part.Properties[2])

                # model.WriteOutput(time*1.0e3)
                # for cond_id in list_conds:
                #     model.model_part.RemoveCondition(cond_id)
                # for node_id in list_nodes:
                #     model.model_part.RemoveNode(node_id)

                ## export to hdf5
                stime = f"{time * 1.0e3:.12g}"
                filename = "mesh_72x72_q4_" + stime + ".h5"
                post_util = MpmHDF5PostUtility(filename)
                post_util.WriteParticles(mpm)
                post_util.WriteParticlesResults(PARTICLE_MASS, mpm)
                post_util.WriteParticlesResults(VELOCITY, mpm)
                post_util.WriteParticlesResults(ACCELERATION, mpm)
                post_util.WriteParticlesResults(STRESSES, mpm)
                itime.write(str(step) + "\t" + str(time) + "\n")

                itime.write(str(step) + "\t" + str(time) + "\n")
                itime.flush()

        step = step + 1

        print("######### TIME STEP " + str(time) + " COMPLETED #########")
        print("#########################################################")

    if logging:
        itime.close()

    return mpm, particle_manager

def test():
    mpm, particle_manager = main(logging = False, output = False, total_time = 0.5)

    monitoring_particle_id = 11
    cx = particle_manager.Particles[monitoring_particle_id].X
    cy = particle_manager.Particles[monitoring_particle_id].Y

    ref_cx = 1.9684075250901925e-01
    ref_cy = 6.1545812425188327e-02

    print("cx: %.16e, cy: %.16e" % (cx, cy))
    assert(abs(cx - ref_cx) < 1e-10)
    assert(abs(cy - ref_cy) < 1e-10)
    print("Test passed")

def tag():
    return "mpm"

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
