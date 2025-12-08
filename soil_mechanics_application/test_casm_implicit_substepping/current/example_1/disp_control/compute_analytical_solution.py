import sys
import math

import analytical_solution_disp_control_axial as analytical_solution

_M = 1.2
_N = 3.0
_R = 2.0
_lambda = 0.066
_kappa = 0.0077
_nu = 0.3
_r = 3.0*(1-2*_nu) / (2*(1+_nu))

load = 200000.0
loadq = 0.0

ana_sol = analytical_solution.proportional_loading_drained_solution(_M, _N, _R, _lambda, _kappa, _nu)
N = 1.788 # intial void ratio

de_a_list = []
for i in range(0, 400):
    de_a_list.append(1.0e-3)

ifile = open("analytical_loading_path.txt", "w")
ifile.write("p q pc ev eq ea er v dphi\n")

p_n = load
q_n = loadq
pc_n = p_n
e_q = 0.0
k = 3.0

e = N - _lambda*math.log(pc_n) - _kappa*math.log(p_n/pc_n)
print("consistent void ratio: " + str(e))
v_n = 1+e # initial specific volume

e_v = 0.0

ifile.write(str(p_n))
ifile.write(" " + str(q_n))
ifile.write(" " + str(pc_n))
ifile.write(" " + str(e_v))
ifile.write(" " + str(e_q))
ifile.write(" " + str(e_v/3+e_q))
ifile.write(" " + str(e_v/3-e_q/2))
ifile.write(" " + str(v_n))
ifile.write(" " + str(0.0))
ifile.write("\n")
ifile.flush()

e_a = 0.0
for step in range(0, len(de_a_list)):
    de_a = de_a_list[step]
    e_a = e_a + de_a

    # compute the trial elastic solution
    p = ana_sol.get_elas_p(v_n, p_n, q_n, k, de_a)
    q = q_n + k*(p - p_n)
    pc = pc_n

    # check the yield state
    is_yield = ana_sol.check_yield(p, q, pc)

    if is_yield == False:
        de_v = ana_sol.get_elas_de_e_v(v_n, p_n, q_n, k, de_a)
        de_q = ana_sol.get_elas_de_e_q(v_n, p_n, q_n, k, de_a)
        dphi = 0.0
    else:
        # compute for yield state
        eta = ana_sol.compute_eta(v_n, p_n, q_n, k, de_a)
        p = ana_sol.get_p(v_n, p_n, q_n, k, eta)
        q = q_n + k*(p - p_n)
        pc = ana_sol.get_pc(v_n, p_n, q_n, k, eta)
        pc_ = ana_sol.get__pc(v_n, p_n, q_n, k, eta, pc_n)
        # print("diff pc: " + str(pc - pc_))
        # print("yield: " + str(ana_sol.compute_f(p, q, pc)))
        de_v = ana_sol.get_de_v(v_n, p_n, q_n, k, eta)
        de_q = ana_sol.get_de_q(v_n, p_n, q_n, k, eta)
        dphi = ana_sol.get_dphi(v_n, p_n, q_n, k, eta)

    e_v = e_v + de_v
    e_q = e_q + de_q
    v = v_n - de_v*v_n # Eq. (16) Peric et al, 2016, here we make it minus because the sample is compressed

    ifile.write(str(p))
    ifile.write(" " + str(q))
    ifile.write(" " + str(pc))
    ifile.write(" " + str(e_v))
    ifile.write(" " + str(e_q))
    ifile.write(" " + str(e_v/3+e_q))
    # print("diff e_a: " + str(e_a - (e_v/3+e_q)))
    ifile.write(" " + str(e_v/3-e_q/2))
    ifile.write(" " + str(v))
    ifile.write(" " + str(dphi))
    ifile.write("\n")
    ifile.flush()

    # update values
    p_n = p
    q_n = q
    v_n = v
    pc_n = pc

ifile.close()
