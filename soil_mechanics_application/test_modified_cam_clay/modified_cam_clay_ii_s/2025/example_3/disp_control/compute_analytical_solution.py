import sys
import math

import analytical_solution_disp_control_axial

_M = 1.2
_lambda = 0.066
_kappa = 0.0077
_nu = 0.3
_r = 3.0*(1-2*_nu) / (2*(1+_nu))

load = 100000.0
loadq = 0.0

ana_sol = analytical_solution_disp_control_axial.proportional_loading_drained_solution(_M, _lambda, _kappa, _nu)
N = 1.788 # intial void ratio

de_a_list = []
for i in range(0, 67):
    de_a_list.append(1.0e-4)
for i in range(0, 9):
    de_a_list.append(1.0e-5)
for i in range(0, 20):
    de_a_list.append(1.0e-5)
for i in range(0, 300):
    de_a_list.append(1.0e-4)
for i in range(0, 300):
    de_a_list.append(1.0e-3)


ifile = open("analytical_loading_path.txt", "w")
ifile.write("p\tq\tpc\tev\teq\tea\ter\tv\n")

p_n = load
q_n = loadq
pc_n = 500000.0
e_q = 0.0
k = 3.0

e = N - _lambda*math.log(pc_n) - _kappa*math.log(p_n/pc_n)
print("consistent void ratio: " + str(e))
v_n = 1+e # initial specific volume

e_v = 0.0

ifile.write(str(p_n))
ifile.write("\t" + str(q_n))
ifile.write("\t" + str(pc_n))
ifile.write("\t" + str(e_v))
ifile.write("\t" + str(e_q))
ifile.write("\t" + str(e_v/3+e_q))
ifile.write("\t" + str(e_v/3-e_q/2))
ifile.write("\t" + str(v_n))
ifile.write("\n")
ifile.flush()

for step in range(0, len(de_a_list)):
    de_a = de_a_list[step]

    # compute the trial elastic solution
    p = ana_sol.get_elas_p(v_n, p_n, q_n, k, de_a)
    q = q_n + k*(p - p_n)
    pc = pc_n

    # check the yield state
    is_yield = ana_sol.check_yield(p, q, pc)

    if is_yield == False:
        de_v = ana_sol.get_elas_de_e_v(v_n, p_n, q_n, k, de_a)
        de_q = ana_sol.get_elas_de_e_q(v_n, p_n, q_n, k, de_a)
    else:
        # compute for yield state
        eta = ana_sol.compute_eta(v_n, p_n, q_n, k, de_a)
        p = ana_sol.get_p(v_n, p_n, q_n, k, eta)
        q = q_n + k*(p - p_n)
        pc = ana_sol.get_pc(v_n, p_n, q_n, k, eta)
        de_v = ana_sol.get_de_v(v_n, p_n, q_n, k, eta)
        de_q = ana_sol.get_de_q(v_n, p_n, q_n, k, eta)

    e_v = e_v + de_v
    e_q = e_q + de_q
    v = v_n - de_v*v_n # Eq. (16) Peric et al, 2016, here we make it minus because the sample is compressed

    ifile.write(str(p))
    ifile.write("\t" + str(q))
    ifile.write("\t" + str(pc))
    ifile.write("\t" + str(e_v))
    ifile.write("\t" + str(e_q))
    ifile.write("\t" + str(e_v/3+e_q))
    ifile.write("\t" + str(e_v/3-e_q/2))
    ifile.write("\t" + str(v))
    ifile.write("\n")
    ifile.flush()

    # update values
    p_n = p
    q_n = q
    v_n = v
    pc_n = pc

ifile.close()
