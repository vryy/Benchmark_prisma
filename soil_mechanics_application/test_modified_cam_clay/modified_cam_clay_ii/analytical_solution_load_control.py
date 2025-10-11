import math

##
## Analytical solution for drained solution and proportional loading, modified Cam-Clay
## Assumption:
##  +   constant Poisson ratio
##  +   the specific volume is constant during a step increment
##  +   the stress path is a line, i.e. \dot{q} = k \dot{p}
## Reference:
##  +   Peric, Analytical solutions for a three-invariant Cam clay model subjected to drained loading histories
##  +   Hoang-Giang Bui, CC note

class proportional_loading_drained_solution:
    def __init__(self, M, lambda_, kappa, poisson_ratio):
        self._M = M
        self._lambda  = lambda_
        self._kappa  = kappa
        self._mu = poisson_ratio
        self._r = 3.0*(1-2*self._mu) / (2*(1+self._mu))

    # Check the yield state of the stress point
    def check_yield(self, p, q, pc, tol = 0):
        return (((q/self._M)**2 + p*(p-pc)) > tol)

    # Compute the preconsolidation pressure in the plastic domain
    # Input:
    #   p       current hydrostatic pressure
    #   q       current deviatoric pressure
    def get_pc(self, p, q):
        pc = p + (q/self._M)**2/p
        return pc

    ######################################################
    #################### LOAD CONTROL ####################
    ######################################################

    # Compute the increment of the elastic volumetric strain in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p       current hydrostatic pressure
    #   p_n     hydrostatic pressure of the previous step
    #   q       current deviatoric pressure
    #   q_n     deviatoric pressure of the previous step
    def get_elas_de_e_v(self, v_n, p, p_n, q, q_n):
        return self.get_de_e_v(v_n, p, p_n, q, q_n)

    # Compute the increment of the elastic deviatoric strain in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p       current hydrostatic pressure
    #   p_n     hydrostatic pressure of the previous step
    #   q       current deviatoric pressure
    #   q_n     deviatoric pressure of the previous step
    def get_elas_de_e_q(self, v_n, p, p_n, q, q_n):
        return self.get_de_e_q(v_n, p, p_n, q, q_n)

    # Compute the increment of the elastic volumetric strain in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p       current hydrostatic pressure
    #   p_n     hydrostatic pressure of the previous step
    #   q       current deviatoric pressure
    #   q_n     deviatoric pressure of the previous step
    def get_de_e_v(self, v_n, p, p_n, q, q_n):
        de_e_v = (self._kappa/v_n)*math.log(p/p_n)
        return de_e_v

    # Compute the increment of the plastic volumetric strain in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p       current hydrostatic pressure
    #   p_n     hydrostatic pressure of the previous step
    #   q       current deviatoric pressure
    #   q_n     deviatoric pressure of the previous step
    def get_de_p_v(self, v_n, p, p_n, q, q_n):
        eta = q/p
        eta_n = q_n/p_n
        c = self._lambda - self._kappa
        de_p_v = (c/v_n)*math.log((p/p_n)*(self._M**2+eta**2)/(self._M**2+eta_n**2))
        return de_p_v

    # Compute the increment of the plastic volumetric strain based on preconsolidation pressure in the plastic domain
    # It is noted that this is only applicable when the loading path is purely plastic
    # Input:
    #   v_n     specific volume of the previous step
    #   pc      current preconsolidation pressure
    #   pc_n    preconsolidation pressure of the previous step
    def get2_de_p_v(self, v_n, pc, pc_n):
        c = self._lambda - self._kappa
        de_p_v = (c/v_n)*math.log(pc/pc_n)
        return de_p_v

    # Compute the increment of the elastic deviatoric strain in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p       current hydrostatic pressure
    #   p_n     hydrostatic pressure of the previous step
    #   q       current deviatoric pressure
    #   q_n     deviatoric pressure of the previous step
    # Output:
    #  de_p_q      increment of the plastic deviatoric strain
    def get_de_e_q(self, v_n, p, p_n, q, q_n):
        k = (q-q_n)/(p-p_n)
        de_e_q = (self._kappa*k/(3*self._r*v_n))*math.log(p/p_n)
        return de_e_q

    # Compute the increment of the plastic deviatoric strain in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p       current hydrostatic pressure
    #   p_n     hydrostatic pressure of the previous step
    #   q       current deviatoric pressure
    #   q_n     deviatoric pressure of the previous step
    # Output:
    #  de_p_q      increment of the plastic deviatoric strain
    def get_de_p_q(self, v_n, p, p_n, q, q_n):
        #   eta     current load ratio (eta=q/p)
        #   eta_n   load ratio of the previous step
        #   k       slope of the loading path (k=dq/dp)
        eta = q/p
        eta_n = q_n/p_n
        k = (q-q_n)/(p-p_n)
        theta = v_n/(self._lambda - self._kappa)
        M = self._M
        A1 = math.log((M - eta) / (M - eta_n))
        A2 = math.log((M + eta) / (M + eta_n))
        A3 = math.log((k - eta) / (k - eta_n))
        A4 = math.atan(eta/M) - math.atan(eta_n/M)
        de_p_q = (1/theta)*(-k/(M*(k-M))*A1 + k/(M*(k+M))*A2 + 2*k/(k**2-M**2)*A3 - 2/M*A4)
        return de_p_q
