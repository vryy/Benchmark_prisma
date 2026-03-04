##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
import quad4_2x2_include
from quad4_2x2_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

def ApplyBC(model):
    ## boundary condition for the background mesh
    tol = 1.0e-6
    for node in model.model_part.Nodes:
        if abs(node.X0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.Fix(VELOCITY_X)
            node.Fix(ACCELERATION_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(VELOCITY_X, 0.0)
            node.SetSolutionStepValue(ACCELERATION_X, 0.0)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.Fix(VELOCITY_Y)
            node.Fix(ACCELERATION_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(ACCELERATION_Y, 0.0)

def main(logging=True, output=True, nsteps=1000, delta_time=1e-4):
    ## create background mesh
    model = quad4_2x2_include.Model('quad4_2x2',os.getcwd()+"/",os.getcwd()+"/")
    model.InitializeModel()

    ## create the particles
    particle_manager = SolidParticleManager()
    query_tool = SpatialGridElementalBinning(0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 1.0e-10)
    query_tool.Initialize(model.model_part.Elements)
    particle_manager.SetQueryTool(query_tool)
    sample_particle = SolidParticleSmallStrain2D(0, 0.0, 0.0)
    lastId = 0
    lastId = particle_manager.AddParticles(model.model_part.Elements, sample_particle, lastId, model.model_part.ProcessInfo)

    monitoring_particle_id = 11
    # print(particle_manager.Particles[monitoring_particle_id])
    ox = particle_manager.Particles[monitoring_particle_id].X
    oy = particle_manager.Particles[monitoring_particle_id].Y

    # sys.exit(0)

    ## create the MP model
    mpm = ULMpmModel(model.model_part, particle_manager)
    mpm.EchoLevel = 3
    print(mpm)

    ApplyBC(model)

    ## create the solver and solve
    import mpm_explicit_strategy
    mpm_solver = mpm_explicit_strategy.MPMExplicitStructuralUSFSolver(mpm)
    mpm_solver.Initialize()

    time = 0.0
    aux_util = MpmAuxiliaryUtility()
    for step in range(0, nsteps):
        time = time + delta_time
        mpm_solver.Solve(time)

        if output:
            last_node_id = aux_util.GetLastNodeId(model.model_part)
            last_cond_id = aux_util.GetLastConditionId(model.model_part)
            print("last_node_id: %d, last_cond_id: %d" % (last_node_id, last_cond_id))
            sample_cond_name = "PointForce2D"
            [list_nodes, list_conds] = particle_manager.ExportParticles(model.model_part, sample_cond_name, last_node_id, last_cond_id, model.model_part.Properties[1])
            model.WriteOutput(time)
            for cond_id in list_conds:
                model.model_part.RemoveCondition(cond_id)
            for node_id in list_nodes:
                model.model_part.RemoveNode(node_id)

        print("######### TIME STEP " + str(time) + " COMPLETED #########")

    # print(particle_manager.Particles[monitoring_particle_id])
    cx = particle_manager.Particles[monitoring_particle_id].X
    cy = particle_manager.Particles[monitoring_particle_id].Y
    print("Particle %d movement: %.16e" % (monitoring_particle_id, math.sqrt((cx-ox)**2 + (cy-oy)**2)))

    return mpm, particle_manager

def test():
    mpm, particle_manager = main(logging = False, output = False, nsteps = 1000, delta_time = 1e-4)

    monitoring_particle_id = 11
    cx = particle_manager.Particles[monitoring_particle_id].X
    cy = particle_manager.Particles[monitoring_particle_id].Y

    ref_cx = 1.0566242403837978e-01
    ref_cy = 8.9433744486633049e-01

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
        main(logging=True, output=True, nsteps=1000, delta_time=1e-4)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
