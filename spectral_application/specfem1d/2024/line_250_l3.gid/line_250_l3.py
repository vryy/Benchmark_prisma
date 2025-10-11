##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Mo 5. Jul 23:39:15 CEST 2021 
##################################################################
import sys
import os
import math
import numpy as np
import time as time_module
##################################################################
##################################################################
sys.path.append('./line_250_l3.gid')
import line_250_l3_include
from line_250_l3_include import *
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

# =====================
# | USER SCRIPT FOR CALCULATION OF EKATE.GID |
# vvvvvvvvvvvvvvvvvvvvv

# ### DEBUGGING: observe the stiffness and mass matrix
# mp_utils = ModelPartUtilities()
# mp_utils.CalculateLocalSystem(model.model_part.Elements[1], model.model_part.ProcessInfo, 1)
# mp_utils.CalculateMassMatrix(model.model_part.Elements[1], model.model_part.ProcessInfo, 1)
# print("----------------")
# mp_utils.CalculateLocalSystem(model.model_part.Elements[2], model.model_part.ProcessInfo, 1)
# mp_utils.CalculateMassMatrix(model.model_part.Elements[2], model.model_part.ProcessInfo, 1)
# print("----------------")
# mp_utils.CalculateLocalSystem(model.model_part.Elements[3], model.model_part.ProcessInfo, 1)
# mp_utils.CalculateMassMatrix(model.model_part.Elements[3], model.model_part.ProcessInfo, 1)
# sys.exit(0)
# ########

class Source:
    """Contains the source properties"""

    def __init__(self, delta_time):
        """Init"""
        self.ampl = 1.0e7
        self.hdur = 100.0 * delta_time  # Duration of the source (s)
        self.decayRate = 2.628
        self.alpha = self.decayRate / self.hdur

    def __getitem__(self, t):
        # """What happens when we do source[t]"""
        # t -= self.hdur
        # return -2 * self.ampl * self.alpha**3 / math.sqrt(math.pi) * (
        #     t * math.exp(-(self.alpha * t)**2))
        """What happens when we do source[t]"""
        t -= self.hdur
        return -2 * self.ampl * self.alpha**3 / np.sqrt(np.pi) * (
            t * np.exp(-(self.alpha * t)**2))

    def plotSource(self):
        """Plot the source"""
        import matplotlib.pyplot as plt
        t = np.linspace(0, self.hdur, 1000)
        plt.figure()
        plt.plot(t, self[t], 'b')
        plt.title('Source(t)')
        plt.grid(True)
        plt.show()

def main(logging=True, output=True, plot=True, nsteps=10000):
    model = line_250_l3_include.Model('line_250_l3',os.getcwd()+"/",os.getcwd()+"/",logging=logging)
    model.InitializeModel()

    nelements = len(model.model_part.Elements)
    cfl = 0.45 # Courant number
    order = 2
    h = 3000.0/nelements/order # mesh size
    E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
    nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
    mu = 0.5*E/(1+nu)
    rho = model.model_part.Properties[1].GetValue(DENSITY)

    vMax = math.sqrt(mu/rho)
    delta_time = cfl*h/vMax

    print("delta_time: " + str(delta_time))

    source = Source(delta_time)
    source_cond = model.model_part.Conditions[1]
    source_node = source_cond.GetNodes()[0]
    # source.plotSource()
    # sys.exit(0)

    if plot:
        grid = []
        for node in model.model_part.Nodes:
            grid.append(node.X0)

        import matplotlib.pyplot as plt
        plt.ion()
        fig = plt.figure()
        # plt.hold(True)

        dplot = 10

    time = 0.0
    for step in range(0, nsteps):

        # adding source term
        source_node.SetSolutionStepValue(SPECTRAL_EXCITATION_Y, source[time])
        time = time + delta_time

        model.SolveModel(time)

        if output:
            model.WriteOutput(time)

        if plot:
            if step % dplot == 0:
                u = []
                for node in model.model_part.Nodes:
                    u.append(node.GetSolutionStepValue(DISPLACEMENT_Y))

                plt.clf()
                plt.plot(grid, u)
                    # plt.ylim([-0.10, 0.10])
                # print("Update plot at it = " + str(it))
                plt.title("it : {}".format(step))
                # plt.draw()
                fig.canvas.draw()
                if step > 0:
                    plt.savefig("plot_"+str(step)+".png")

    if logging:
        ifile = open("results.txt", "w")
        ifile.write("node_x\ta\tv\tu\n")
        for node in model.model_part.Nodes:
            a = node.GetSolutionStepValue(ACCELERATION_Y)
            v = node.GetSolutionStepValue(DISPLACEMENT_DT_Y)
            u = node.GetSolutionStepValue(DISPLACEMENT_Y)
            ifile.write(str(node.X0) + "\t" + str(a) + "\t" + str(v) + "\t" + str(u) + "\n")
        ifile.close()

    return model

def test():
    model = main(logging=False, output=False, plot=False, nsteps=100)

    monitoring_node = model.model_part.Nodes[499]
    a = monitoring_node.GetSolutionStepValue(ACCELERATION_Y)
    v = monitoring_node.GetSolutionStepValue(DISPLACEMENT_DT_Y)
    u = monitoring_node.GetSolutionStepValue(DISPLACEMENT_Y)

    ref_a = -47363.99084454844
    ref_v = 208.10877001448884
    ref_u = 21.50164817032548

    assert(abs(a - ref_a) / abs(ref_a) < 1e-12)
    assert(abs(v - ref_v) / abs(ref_v) < 1e-12)
    assert(abs(u - ref_u) / abs(ref_u) < 1e-12)
    print("Test passed")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True, output=False, plot=True, nsteps=100)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
