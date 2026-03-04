import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FiniteCellApplication import *
from KratosMultiphysics.MaterialPointApplication import *

import finite_cell_simulator
from finite_cell_simulator import *

import cbgeo_mpm_writer

class Simulator(FiniteCellSimulator, object):

    def __init__(self, params):
        super(Simulator, self).__init__(params)
        ########Create the B-Rep representing the boundary
        ls1 = CircularLevelSet(0.4, 0.4, 0.2)
        ls2 = CircularLevelSet(0.8, 0.8, 0.2)
        # self.brep = ls2
        self.brep = UnionLevelSet(ls1, ls2)
        #########Some default parameters##################
        if "ignore_invalid_particles" not in self.params:
            self.params["ignore_invalid_particles"] = False
        if "write_output_per_each_step" not in self.params:
            self.params["write_output_per_each_step"] = False
        if "write_cbgeo" not in self.params:
            self.params["write_cbgeo"] = False
        if "hit_angle" not in self.params:
            self.params["hit_angle"] = 0.0
        if "delta_time" not in self.params:
            self.params["delta_time"] = 0.0051
        if "total_time" not in self.params:
            self.params["total_time"] = 2.5
        if "mpm_scheme" not in self.params:
            self.params["mpm_scheme"] = "USL"
        if "add_rigid_wall" not in self.params:
            self.params["add_rigid_wall"] = True
        if "particle_definition" not in self.params:
            self.params["particle_definition"] = "solid_particle_large_strain"
        if "velocity_magnitude" not in self.params:
            self.params["velocity_magnitude"] = 0.2
        if "write_output" not in self.params:
            self.params["write_output"] = False
        if "write_energy" not in self.params:
            self.params["write_energy"] = True
        if "mpm_echo_level" not in self.params:
            self.params["mpm_echo_level"] = 1
        if "search_tolerance" not in self.params:
            self.params["search_tolerance"] = 1e-10

    ###SIMULATION DRIVER#############
    def Initialize(self, model_solid):
        ######## Extract the solid elements
        self.aux_util = FiniteCellAuxiliaryUtility()
        self.mpm_aux_util = MpmAuxiliaryUtility()
        self.solid_elements = ElementsArray()
        for elem_id in model_solid.layer_sets["soil"]:
            if elem_id in model_solid.model_part.Elements:
                elem = model_solid.model_part.Elements[elem_id]
                self.aux_util.AddElement(self.solid_elements, elem)
#        self.aux_util.GetElements(solid_elements, model_solid.model_part, model_solid.layer_sets['solid'])
        print("len(solid_elements):", len(self.solid_elements))

        super(Simulator, self).Initialize(model_solid, self.solid_elements)

    def Run(self, model_solid):
        self.Initialize(model_solid)
        # model_solid.WriteOutput(0.0)
        if self.params["write_cbgeo"]:
            cbgeo_writer = cbgeo_mpm_writer.CbGeoMPMWriter(2)
            cbgeo_writer.WriteMesh(model_solid.model_part, "mesh.txt")
            cbgeo_writer.WriteParticles(model_solid.model_part, "particles.txt")
            cbgeo_writer.WriteParticlesVolumes(model_solid.model_part, "particles-volumes.txt")

            tol = 1.0e-6
            layer_particle_sets = {}
            layer_node_sets = {}
            layer_node_sets[0] = []
            layer_node_sets[1] = []
            layer_node_sets[2] = []
            for node in model_solid.model_part.Nodes:
                if abs(node.X0) < tol:
                    layer_node_sets[0].append(node.Id)
                if abs(node.X0-75.0) < tol:
                    layer_node_sets[1].append(node.Id)
                if abs(node.Y0) < tol:
                    layer_node_sets[2].append(node.Id)
            cbgeo_writer.WriteEntitySets(layer_particle_sets, layer_node_sets, "entity_sets.json")
            sys.exit(0)
        ##################################################################
        ###  SIMULATION  #################################################
        ##################################################################

        ## initialize the particles
        ball1_elements = ElementsArray()
        ball2_elements = ElementsArray()
        for elem in model_solid.model_part.Elements:
            elem_center = [0.0, 0.0]
            for node in elem.GetNodes():
                elem_center[0] = elem_center[0] + node.X0
                elem_center[1] = elem_center[1] + node.Y0
            elem_center[0] = elem_center[0] / len(elem.GetNodes())
            elem_center[1] = elem_center[1] / len(elem.GetNodes())
            if elem_center[0] < 0.6 and elem.Is(ACTIVE):
                ball1_elements.append(elem)
            if elem_center[0] > 0.6 and elem.Is(ACTIVE):
                ball2_elements.append(elem)
        print("len(ball1_elements): " + str(len(ball1_elements)))
        print("len(ball2_elements): " + str(len(ball2_elements)))

        ## create the particle container
        query_tool = NonuniformStructuredGridElementalIndexing2D(self.params['search_tolerance'])
        query_tool.Initialize(model_solid.model_part.Elements)

        if "solid_particle" in self.params["particle_definition"]:
            if self.params["particle_definition"] == "solid_particle_small_strain":
                sample_particle = SolidParticleSmallStrain2D(0, 0.0, 0.0)
            elif self.params["particle_definition"] == "solid_particle_large_strain":
                sample_particle = SolidParticleLargeStrain2D(0, 0.0, 0.0)
            elif self.params["particle_definition"] == "solid_particle_large_strain_t2":
                sample_particle = SolidParticleLargeStrain2DT2(0, 0.0, 0.0)
            else:
                raise Exception("Unknown particle " + str(self.params["particle_definition"]))

            particle_manager1 = SolidParticleManager()
            particle_manager1.SetQueryTool(query_tool)
            particle_manager1.EchoLevel = 1
            lastId = 0
            lastId = particle_manager1.AddParticles(ball1_elements, sample_particle, lastId, model_solid.model_part.ProcessInfo)

            if self.params["number_of_balls"] > 1:
                particle_manager2 = SolidParticleManager()
                particle_manager2.SetQueryTool(query_tool)
                particle_manager2.EchoLevel = 1
                lastId = particle_manager2.AddParticles(ball2_elements, sample_particle, lastId, model_solid.model_part.ProcessInfo)
        elif self.params["particle_definition"] == "solid_cpdi_large_strain":
            sample_particle = SolidCPDILargeStrain2D(0, 0.0, 0.0)

            particle_manager1 = SolidCPDIManager2D()
            particle_manager1.SetQueryTool(query_tool)
            particle_manager1.EchoLevel = 1
            lastId = 0
            lastId = particle_manager1.AddParticles(ball1_elements, model_solid.model_part.Elements, sample_particle, lastId, model_solid.model_part.ProcessInfo)

            if self.params["number_of_balls"] > 1:
                particle_manager2 = SolidCPDIManager2D()
                particle_manager2.SetQueryTool(query_tool)
                particle_manager2.EchoLevel = 1
                lastId = particle_manager2.AddParticles(ball2_elements, model_solid.model_part.Elements, sample_particle, lastId, model_solid.model_part.ProcessInfo)

        if self.params["add_rigid_wall"]:
            particle_manager3 = SolidParticleManager()
            particle_manager3.SetQueryTool(query_tool)
            particle_manager3.EchoLevel = 1
            mesh_size_x = 1.2/self.params['mesh_division_x']
            mesh_size_y = 1.2/self.params['mesh_division_y']
            wall_prop = model_solid.model_part.Properties[3]
            for elem in model_solid.model_part.Elements:
                elem_center = [0.0, 0.0]
                for node in elem.GetNodes():
                    elem_center[0] = elem_center[0] + node.X0
                    elem_center[1] = elem_center[1] + node.Y0
                elem_center[0] = elem_center[0] / len(elem.GetNodes())
                elem_center[1] = elem_center[1] / len(elem.GetNodes())

                if (elem_center[0] < mesh_size_x) or (elem_center[0] > 1.2-mesh_size_x) or (elem_center[1] < mesh_size_y) or (elem_center[1] > 1.2-mesh_size_y):
                    particle = RigidParticle(lastId + 1, elem_center[0], elem_center[1], 0.0)
                    particle_manager3.AddParticle(particle, wall_prop, model_solid.model_part.Elements, model_solid.model_part.ProcessInfo, True, self.params["ignore_invalid_particles"])
                    lastId = lastId + 1

        ## initial condition
        for node in model_solid.model_part.Nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(VELOCITY_X, 0.0)
            node.SetSolutionStepValue(VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(ACCELERATION_X, 0.0)
            node.SetSolutionStepValue(ACCELERATION_Y, 0.0)

        angle = 45.0 + self.params['hit_angle']
        c = math.cos(angle/180.0*math.pi)
        s = math.sin(angle/180.0*math.pi)
        velo_magnitude = self.params['velocity_magnitude']
        for p in particle_manager1.Particles:
            p.Vx = c*velo_magnitude
            p.Vy = s*velo_magnitude
        if self.params["number_of_balls"] > 1:
            for p in particle_manager2.Particles:
                p.Vx = -c*velo_magnitude
                p.Vy = -s*velo_magnitude

        if "_t2" in self.params["particle_definition"]:
            for p in particle_manager1.Particles:
                p.Vx1 = c*velo_magnitude
                p.Vy1 = s*velo_magnitude
            if self.params["number_of_balls"] > 1:
                for p in particle_manager2.Particles:
                    p.Vx1 = -c*velo_magnitude
                    p.Vy1 = -s*velo_magnitude

        ## create the MP model
        mpm = ULMultiMpmModel(model_solid.model_part)
        mpm.IgnoredFlag = self.params["ignore_invalid_particles"]
        mpm.AddParticleManager(1, particle_manager1)
        if self.params["number_of_balls"] > 1:
            mpm.AddParticleManager(2, particle_manager2)
        if self.params["add_rigid_wall"]:
            mpm.AddParticleManager(3, particle_manager3)
        if 'friction_coefficient' in self.params:
            mpm.FrictionCoefficient = self.params['friction_coefficient']
        mpm.Check()
        mpm.EchoLevel = self.params['mpm_echo_level']
        print(mpm)

        ## create the solver and solve
        import mpm_explicit_strategy
        if self.params['mpm_scheme'] == "USL":
            mpm_solver = mpm_explicit_strategy.MPMExplicitContactUSLSolver(mpm)
        elif self.params['mpm_scheme'] == "USF":
            mpm_solver = mpm_explicit_strategy.MPMExplicitContactUSFSolver(mpm)
        elif self.params['mpm_scheme'] == "MUSL":
            mpm_solver = mpm_explicit_strategy.MPMExplicitContactMUSLSolver(mpm)
        mpm_solver.Initialize()

        time = 0.0
        delta_time = self.params['delta_time']
        total_time = self.params['total_time']

        if self.params['write_energy']:
            ifile = open("monitoring_energy.txt", "w")
            ifile.write("step\ttime\tstrain_energy\tkinetic_energy\ttotal_energy\n")

        # activate all background elements for post-processing
        for elem in model_solid.model_part.Elements:
            elem.Set(ACTIVE, True)

        step = 0
        sample_output = self.params['sample_output']
        write_output = self.params['write_output']
        sample_cond_name = "PointForce2D"
        while time < total_time:
            time = time + delta_time
            print("#########################################################")
            print("######### TIME STEP " + str(time) + " BEGIN #############")

            mpm_solver.Solve(time)

            if step % sample_output == 0:
                if write_output:
                    # export to model_part
                    lastNodeId = self.mpm_aux_util.GetLastNodeId(model_solid.model_part)
                    lastCondId = self.mpm_aux_util.GetLastConditionId(model_solid.model_part)
                    [list_nodes1, list_conds1] = particle_manager1.ExportParticles(model_solid.model_part, sample_cond_name, lastNodeId, lastCondId, model_solid.model_part.Properties[2])
                    if self.params["number_of_balls"] > 1:
                        lastNodeId = self.mpm_aux_util.GetLastNodeId(model_solid.model_part)
                        lastCondId = self.mpm_aux_util.GetLastConditionId(model_solid.model_part)
                        [list_nodes2, list_conds2] = particle_manager2.ExportParticles(model_solid.model_part, sample_cond_name, lastNodeId, lastCondId, model_solid.model_part.Properties[3])
                    if self.params["add_rigid_wall"]:
                        lastNodeId = self.mpm_aux_util.GetLastNodeId(model_solid.model_part)
                        lastCondId = self.mpm_aux_util.GetLastConditionId(model_solid.model_part)
                        [list_nodes3, list_conds3] = particle_manager3.ExportParticles(model_solid.model_part, sample_cond_name, lastNodeId, lastCondId, model_solid.model_part.Properties[4])
                    model_solid.WriteOutput(time)
                    for cond_id in list_conds1:
                        model_solid.model_part.RemoveCondition(cond_id)
                    for node_id in list_nodes1:
                        model_solid.model_part.RemoveNode(node_id)
                    if self.params["number_of_balls"] > 1:
                        for cond_id in list_conds2:
                            model_solid.model_part.RemoveCondition(cond_id)
                        for node_id in list_nodes2:
                            model_solid.model_part.RemoveNode(node_id)
                    if self.params["add_rigid_wall"]:
                        for cond_id in list_conds3:
                            model_solid.model_part.RemoveCondition(cond_id)
                        for node_id in list_nodes3:
                            model_solid.model_part.RemoveNode(node_id)

                # compute the energy
                strain_energy = 0.0
                kinetic_energy = 0.0
                for p in particle_manager1.Particles:
                    se = p.GetValue(STRAIN_ENERGY, model_solid.model_part.ProcessInfo)
                    ke = p.GetValue(KINETIC_ENERGY, model_solid.model_part.ProcessInfo)
                    strain_energy = strain_energy + se
                    kinetic_energy = kinetic_energy + ke
                if self.params["number_of_balls"] > 1:
                    for p in particle_manager2.Particles:
                        se = p.GetValue(STRAIN_ENERGY, model_solid.model_part.ProcessInfo)
                        ke = p.GetValue(KINETIC_ENERGY, model_solid.model_part.ProcessInfo)
                        strain_energy = strain_energy + se
                        kinetic_energy = kinetic_energy + ke

                if self.params['write_energy']:
                    ifile.write(str(step) + '\t' + str(time) + '\t' + str(strain_energy) + '\t' + str(kinetic_energy) + '\t' + str(strain_energy+kinetic_energy) + '\n')
                    ifile.flush()

                # ## export to hdf5
                # filename = "mesh_12x12_q4_" + str(time) + ".h5"
                # post_util = MpmHDF5PostUtility(filename)
                # post_util.WriteParticles(mpm)
                # post_util.WriteParticlesResults(PARTICLE_MASS, mpm)
                # post_util.WriteParticlesResults(VELOCITY, mpm)
                # post_util.WriteParticlesResults(ACCELERATION, mpm)
                # post_util.WriteParticlesResults(THREED_STRESSES, mpm)
                # itime.write(str(step) + "\t" + str(time) + "\n")

            step = step + 1

            print("######### TIME STEP " + str(time) + " COMPLETED #########")
            print("#########################################################")

        if self.params['write_energy']:
            ifile.close()

        print("ANALYSIS COMPLETED")
        # print(Timer)

        particle_manager_list = []
        particle_manager_list.append(particle_manager1)
        if self.params["number_of_balls"] > 1:
            particle_manager_list.append(particle_manager2)

        energies = []
        energies.append(strain_energy)
        energies.append(kinetic_energy)

        return particle_manager_list, energies
