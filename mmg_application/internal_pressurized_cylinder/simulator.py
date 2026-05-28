import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MultithreadedSolversApplication import *

import mmg_simulator

###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
def Compute_L2_error(model, solution, P):
    nom = 0.0
    denom = 0.0
    for element in model.model_part.Elements:
        u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        er = []
        for i in range(0, len(u)):
            ana_u = solution.get_displacement(P, Q[i][0], Q[i][1])
            d = (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2))
            er.append(math.sqrt(d))
            nom = nom + d * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
        element.SetValuesOnIntegrationPoints(L2_ERROR, er, model.model_part.ProcessInfo)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

###COMPUTE GLOBAL DISPLACEMENT (H1) ERROR###
def Compute_H1_error(model, solution, P):
    nom = 0.0
    denom = 0.0
    for element in model.model_part.Elements:
        o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        er = []
        for i in range(0, len(o)):
            ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1])
            d = (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2)))
            er.append(math.sqrt(d))
            nom = nom + d * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        element.SetValuesOnIntegrationPoints(H1_ERROR, er, model.model_part.ProcessInfo)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

###COMPUTE GLOBAL DISPLACEMENT (H1) ERROR FOR RECOVERY STRESS###
def Compute_H1_error_rstress(model, solution, P):
    nom = 0.0
    denom = 0.0
    nu = float(model.model_part.Properties[1].GetValue(POISSON_RATIO))
    for element in model.model_part.Elements:
        o = element.GetValuesOnIntegrationPoints(RECOVERY_STRESSES, model.model_part.ProcessInfo)
        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        er = []
        for i in range(0, len(o)):
            o_zz = nu * (o[i][0] + o[i][1])
            ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1])
            d = (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o_zz - ana_o[2], 2) + 2.0*(pow(o[i][2] - ana_o[3], 2) + pow(0.0 - ana_o[4], 2) + pow(0.0 - ana_o[5], 2)))
            er.append(math.sqrt(d))
            nom = nom + d * W[i][0] * J0[i][0]
            denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        element.SetValuesOnIntegrationPoints(H1_ERROR, er, model.model_part.ProcessInfo)
    if denom == 0.0:
        if nom == 0.0:
            return 0.0
        else:
            return float('nan');
    else:
        return math.sqrt(abs(nom / denom))

class Simulator:

    def __init__(self, params):
        self.params = params
        #########Some default parameters##################
        if "write_output_per_each_step" not in self.params:
            self.params["write_output_per_each_step"] = False
        if "compute_condition_number" not in self.params:
            self.params["compute_condition_number"] = False

    ###SIMULATION DRIVER#############
    def Initialize(self, model):
        # register the boundary
        self.a = 100.0
        self.b = 200.0
        tol = 1.0e-4
        ########Create the B-Rep representing the boundary
        ils = CircularLevelSet(0.0, 0.0, self.a)
        ols = CircularLevelSet(0.0, 0.0, self.b)
        ils.Tolerance = tol*self.a
        ols.Tolerance = tol*self.b
        self.params['boundaries'] = {}
        self.params['boundaries'][1] = ils
        self.params['boundaries'][2] = ols

        # create mmg model
        self.mmg_model = mmg_simulator.CreateMmgModel(model.model_part, self.params)

    def Update(self, model):
        # create the model_part
        [mp, layer_sets, layer_cond_sets, node_groups, element_assignments] = mmg_simulator.ConstructSystemModelPart(self.mmg_model, self.params)
        model.SetModelPart(mp)
        model.layer_sets = layer_sets
        model.layer_cond_sets = layer_cond_sets
        model.node_groups = node_groups
        model.InitializeModel(stress_recovery_type = self.params["stress_recovery_type"], neighbour_expansion_level=self.params['neighbour_expansion_level'])

        # boundary condition
        tol = 1.0e-6
        for node in model.model_part.Nodes:
            if abs(node.X0 - 0.0) < tol:
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            if abs(node.Y0 - 0.0) < tol:
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)

        model.solver.solver.MoveMeshFlag = False

    def Run(self, model, time=1.0, logging=True):
        ##################################################################
        ###  SIMULATION  #################################################
        ##################################################################
        analytical_solution = self.params["analytical_solution"]
        report_conv_name = self.params["report_convergence_name"]
        report_disp_b_name = self.params["report_disp_b_name"]
        compute_condition_number = self.params["compute_condition_number"]
        stress_recovery_type = self.params["stress_recovery_type"]

        if compute_condition_number:
            model.solver.solver.ReformDofSetAtEachStep = False
            eigen_solver = SpectraEigenvaluesSolver()

        if stress_recovery_type > 0:
            FindElementalNeighboursProcess(model.model_part, 2, 6).Execute()

        ## Loading path
        E = float(model.model_part.Properties[1].GetValue(YOUNG_MODULUS))
        nu = float(model.model_part.Properties[1].GetValue(POISSON_RATIO))
        oy = float(model.model_part.Properties[1].GetValue(TENSILE_STRENGTH))
        print("E:", E)
        print("oy:", oy)

        self.elastic_sol = analytical_solution.ElasticSolution(E, nu, oy, self.a, self.b)

        P0 = self.elastic_sol.computePO()
        Plim = self.elastic_sol.computePlim()
        print("P0:", P0)
        print("Plim:", Plim)

        self.elastic_loading_path = [0.5*P0]

        # model.WriteOutput(0.0)
        # sys.exit(0)

        if logging:
            ###############Reporting
            ifile_disp_y_0_b = open(report_disp_b_name, "w")
            ifile_report = open(report_conv_name, "w")

            #ifile_disp_y_0_a.write("Load\t\tdisp-y\t\tana_disp-y\t\trel_error\t\tana_plastic_front\n")
            ifile_disp_y_0_b.write("Load\t\tdisp-y\t\tana_disp-y\t\trel_error\t\tana_plastic_front\n")
            ifile_report.write("Load\t\tDelta_Load\t\tl2_error\t\th1_error\n")

        Pold = 0.0
        disp_y_0_b = []
        disp_y_0_a = []

        tol = 1.0e-6
        node_0_b_0_existed = False
        #node_0_a_0_existed = False
        for node in model.model_part.Nodes:
            if abs(node.X0 - 0.0) < tol and abs(node.Y0 - self.b) < tol and abs(node.Z0 - 0.0) < tol:
                node_0_b_0 = node
                node_0_b_0_existed = True
        #    if abs(node.X0 - 0.0) < tol and abs(node.Y0 - self.a) < tol and abs(node.Z0 - 0.0) < tol:
        #        node_0_a_0 = node
        #        node_0_a_0_existed = True
        print("node_0_b_0_existed:", node_0_b_0_existed)
        #print("node_0_a_0_existed:", node_0_a_0_existed)
        #sys.exit(0)

        analysis_time = 0.0

        compute_error = True

        for P in self.elastic_loading_path:
            start_time = time_module.time()

            print("######LOADING P = " + str(P) + "##################")
            ## Boundary load
            if self.params['condition_name'] == "LinePressure":
                for cond_id in model.layer_cond_sets["inner"]:
                    condition = model.model_part.Conditions[cond_id]
                    condition.SetValue(PRESSURE, P)
            elif self.params['condition_name'] == "LineForce":
                for node_id in model.node_groups['inner']:
                    node = model.model_part.Nodes[node_id]
                    x = node.X0
                    y = node.Y0
                    r = math.sqrt(x*x + y*y)
                    t = math.acos(x/r)
                    node.SetSolutionStepValue(FACE_LOAD, [P * math.cos(t), P * math.sin(t), 0.0])
            print("setup loading P = " + str(P) + " completed")
            ############################################

            model.Solve(time, 0, 0, 0, 0)

            end_time = time_module.time()
            analysis_time = analysis_time + (end_time - start_time)

            if compute_condition_number:
                max_evalues = eigen_solver.SolveLargestSym(model.solver.solver.A, 3, 5)
                print("largest eigenvalues:", max_evalues)
                min_evalues = eigen_solver.SolveSmallestSPD(model.solver.solver.A, SuperLUSolver(), 3, 5)
                print("smallest eigenvalues:", min_evalues)
                print("condition number:", max_evalues[0]/min_evalues[2])

            if node_0_b_0_existed:
                node = node_0_b_0
                disp_y_nume = node.GetSolutionStepValue(DISPLACEMENT_Y)
                disp_y_0_b.append(disp_y_nume)
                disp_0_b_ana = self.elastic_sol.get_displacement(P, node.X0, node.Y0, node.Z0)
                norm_disp_0_b_ana = math.sqrt(math.pow(disp_0_b_ana[0], 2) + math.pow(disp_0_b_ana[1], 2) + math.pow(disp_0_b_ana[2], 2))
                if norm_disp_0_b_ana == 0.0:
                    norm_disp_0_b_ana = 1.0
                print("displacement at (0, b): ", node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(DISPLACEMENT_Z))
                print("analytical displacement at (0, b): ", disp_0_b_ana[0], disp_0_b_ana[1], disp_0_b_ana[2])
                if logging:
                    ifile_disp_y_0_b.write(str(P) + "\t" + str(disp_y_nume) + "\t" + str(disp_0_b_ana[1]) + "\t" + str(disp_0_b_ana[2]) + "\t" + str(abs(disp_y_nume-disp_0_b_ana[1])/norm_disp_0_b_ana) + "\t" + str(self.a) + "\n")

        #    if node_0_a_0_existed:
        #        node = node_0_a_0
        #        disp_y_nume = node.GetSolutionStepValue(DISPLACEMENT_Y)
        #        disp_y_0_a.append(disp_y_nume)
        #        disp_0_a_ana = self.elastic_sol.get_displacement(P, node.X0, node.Y0, node.Z0)
        #        norm_disp_0_a_ana = math.sqrt(math.pow(disp_0_a_ana[0], 2) + math.pow(disp_0_a_ana[1], 2) + math.pow(disp_0_a_ana[2], 2))
        #        if norm_disp_0_a_ana == 0.0:
        #            norm_disp_0_a_ana = 1.0
        #        print("displacement at (0, a): ", node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(DISPLACEMENT_Z))
        #        print("analytical displacement at (0, a): ", disp_0_a_ana[0], disp_0_a_ana[1], disp_0_a_ana[2])
        #        ifile_disp_y_0_a.write(str(P) + "\t" + str(disp_y_nume) + "\t" + str(disp_0_a_ana[1]) + "\t" + str(disp_0_a_ana[2]) + "\t" + str(abs(disp_y_nume-disp_0_a_ana[1])/norm_disp_0_a_ana) + "\t" + str(self.a) + "\n")

            if compute_error:
                ## compute the L2 error
                if P == 0.0:
                    l2_error = 0.0
                    h1_error = 0.0
                    h1_error_rs = 0.0
                else:
                    l2_error = Compute_L2_error(model, self.elastic_sol, P)
                    h1_error = Compute_H1_error(model, self.elastic_sol, P)
                    if stress_recovery_type > 0:
                        h1_error_rs = Compute_H1_error_rstress(model, self.elastic_sol, P)
                    else:
                        h1_error_rs = 0.0

                    print("Global displacement (L2) error: %.16e" % l2_error)
                    print("Global displacement (H1) error: %.16e" % h1_error)
                    print("Global displacement (H1) error (recovery stress): %.16e" % h1_error_rs)

                ## reporting
                if logging:
                    ifile_report.write(str(P) + "\t" + str(P - Pold) + "\t" + str(l2_error) + "\t" + str(h1_error) + "\n")

            if self.params["write_output_per_each_step"]:
                model.WriteOutput(time)

            Pold = P

        print("Analysis time: " + str(analysis_time) + " s")

        if logging:
            #ifile_disp_y_0_a.close()
            ifile_disp_y_0_b.close()
            ifile_report.close()

        model.l2_error = l2_error
        model.h1_error = h1_error
        model.h1_error_rs = h1_error_rs
