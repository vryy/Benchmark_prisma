import math

##
## Analytical solution for drained solution and proportional loading, CASM
## Assumption:
##  +   constant Poisson ratio
##  +   the specific volume is constant during a step increment
##  +   the stress path is a line, i.e. \dot{q} = k \dot{p}
## Reference:
##  +   Hoang-Giang Bui, CASM note

## compute integral \int x^N/(1-x), x=x1..x2
## N must be positive integer
def F(N, x1, x2):
    if N == 0:
        return math.log((1-x1)/(1-x2))
    else:
        return F(N-1, x1, x2) + (x1**N - x2**N)/N

class proportional_loading_drained_solution:
    def __init__(self, M, N, R, lambda_, kappa, poisson_ratio):
        self._M = M
        self._N = N
        self._R = R
        self._lambda  = lambda_
        self._kappa  = kappa
        self._nu = poisson_ratio
        self._r = 3.0*(1-2*self._nu) / (2*(1+self._nu))

    # Check the yield state of the stress point
    def check_yield(self, p, q, pc, tol = 0):
        return ((q/(self._M*p))**self._N + math.log(p/pc) / math.log(self._R) > tol)

    # Compute the preconsolidation pressure in the plastic domain
    # Input:
    #   p       current hydrostatic pressure
    #   q       current deviatoric pressure
    def get_pc(self, p, q):
        pc = p * math.exp(math.log(self._R)*(q/(self._M*p))**self._N)
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
        k = (q-q_n)/(p-p_n)
        N = self._N
        M = self._M
        R = self._R
        theta = v_n/(self._lambda - self._kappa)
        FN1 = F(N-1, eta_n/k, eta/k)
        FN = F(N, eta_n/k, eta/k)
        de_p_v = (1.0/theta)*(math.log(p/p_n) + ((k/M)**N)*N*math.log(R)*(FN1 - FN))
        return de_p_v

    # Compute the increment of the plastic volumetric strain based on preconsolidation pressure in the plastic domain
    # It is noted that this is only applicable when the loading path is purely plastic
    # Input:
    #   v_n     specific volume of the previous step
    #   pc      current preconsolidation pressure
    #   pc_n    preconsolidation pressure of the previous step
    def get2_de_p_v(self, v_n, pc, pc_n):
        theta = v_n/(self._lambda - self._kappa)
        de_p_v = (1.0/theta)*math.log(pc/pc_n)
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
        de_e_q = (self._kappa*k/(3.0*self._r*v_n))*math.log(p/p_n)
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
        N = self._N
        R = self._R
        A = (2.0*(M**2) - 3.0*(3.0+M))/(M-k)
        B = (3.0*(3.0+M) - 2.0*k*M)/(M-k)
        FN1 = F(N-1, eta_n/M, eta/M)
        FN = F(N, eta_n/M, eta/M)
        de_p_v = self.get_de_p_v(v_n, p, p_n, q, q_n)
        de_p_q = A/(9.0*theta)*(math.log((M-eta_n)/(M-eta)) + N*math.log(R)*(k/M*FN1 - FN)) + B/9.0*de_p_v
        return de_p_q
