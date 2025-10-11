import sys
import math

import analytical_solution

_M = 1.2
_lambda = 0.066
_kappa = 0.0077
_nu = 0.3

load = 200000.0
loadq = 0.0

ana_sol = analytical_solution.proportional_loading_drained_solution(_M, _lambda, _kappa, _nu)
N = 1.788 # intial void ratio

delta_load_list = []
delta_loadq_list = []
for i in range(0, 266):
    delta_load_list.append(5.0e2)
    delta_loadq_list.append(15.0e2)
# for i in range(0, 200):
#     delta_load_list.append(5.0e2/10)
#     delta_loadq_list.append(15.0e2/10)
# for i in range(0, 350):
#     delta_load_list.append(5.0e2/20)
#     delta_loadq_list.append(15.0e2/20)
# for i in range(0, 400):
#     delta_load_list.append(5.0e2/40)
#     delta_loadq_list.append(15.0e2/40)

ifile = open("analytical_loading_path.txt", "w")
ifile.write("p\tq\tpc\tev\teq\tea\ter\tv\n")

p_n = load
q_n = loadq
pc_n = load
e_q = 0.0
e_v = 0.0

e = N - _lambda*math.log(pc_n) - _kappa*math.log(p_n/pc_n)
print("consistent void ratio: " + str(e))
v_n = 1+e # initial specific volume

ifile.write("\t" + str(p_n))
ifile.write("\t" + str(q_n))
ifile.write("\t" + str(pc_n))
ifile.write("\t" + str(e_v))
ifile.write("\t" + str(e_q))
ifile.write("\t" + str(e_v/3+e_q))
ifile.write("\t" + str(e_v/3-e_q/2))
ifile.write("\t" + str(v_n))
ifile.write("\n")
ifile.flush()

for step in range(0, len(delta_load_list)):
    load = load + delta_load_list[step]
    loadq = loadq + delta_loadq_list[step]

    p = load
    q = loadq
    pc = ana_sol.get_pc(p, q)

    de_e_v = ana_sol.get_de_e_v(v_n, p, p_n, q, q_n)
    de_p_v = ana_sol.get_de_p_v(v_n, p, p_n, q, q_n)
    de_e_q = ana_sol.get_de_e_q(v_n, p, p_n, q, q_n)
    de_p_q = ana_sol.get_de_p_q(v_n, p, p_n, q, q_n)
    e_v = e_v + de_e_v + de_p_v
    e_q = e_q + de_e_q + de_p_q
    v = v_n - (de_e_v + de_p_v)*v_n # Eq. (16) Peric et al, 2016, here we make it minus because the sample is compressed

    ifile.write("\t" + str(p))
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
