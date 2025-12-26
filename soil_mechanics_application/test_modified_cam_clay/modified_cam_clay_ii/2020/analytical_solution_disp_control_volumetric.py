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

    ##############################################################
    #################### DISPLACEMENT CONTROL ####################
    ##############################################################

    # Compute the increment of the deviatoric strain in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_elas_de_e_q(self, v_n, p_n, q_n, k, de_v):
        de_q = k/(3.0*self._r) * de_v
        return de_q

    # Compute the hydrostatic pressure in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_elas_p(self, v_n, p_n, q_n, k, de_v):
        p = p_n*math.exp(v_n/self._kappa * de_v)
        return p

    # Compute the eta in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def compute_eta(self, v_n, p_n, q_n, k, de_v):
        eta_n = q_n/p_n
        theta = v_n/(self._lambda - self._kappa)
        # print(eta_n)
        C = de_v + 1/theta*math.log(eta_n**2+self._M**2) - (self._kappa/v_n+1/theta)*math.log(k-eta_n)
        eta = eta_n
        max_iters = 30
        it = 1
        converged = False
        while (converged == False) and (it < max_iters):
            l = 1/theta*math.log(eta**2+self._M**2) - (self._kappa/v_n+1/theta)*math.log(k-eta) - C
            if abs(l) < 1.0e-6:
                break
            dl = 1/theta*(2.0*eta)/(eta**2+self._M**2) + (self._kappa/v_n+1/theta)/(k-eta)
            deta = l/dl
            eta = eta - deta
            if abs(deta) < 1.0e-8:
                converged = True
            it = it + 1
        return eta

    # Compute the hydrostatic pressure in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_p(self, v_n, p_n, q_n, k, eta):
        eta_n = q_n/p_n
        p = p_n*(k-eta_n)/(k-eta)
        return p

    # Compute the preconsolidation pressure in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_pc(self, v_n, p_n, q_n, k, eta):
        p = self.get_p(v_n, p_n, q_n, k, eta)
        pc = p*(1 + (eta/self._M)**2)
        return pc

    # Compute the increment of deviatoric strain in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_de_q(self, v_n, p_n, q_n, k, eta):
        eta_n = q_n/p_n
        theta = v_n/(self._lambda - self._kappa)
        M = self._M
        A1 = math.log((M - eta) / (M - eta_n))
        A2 = math.log((M + eta) / (M + eta_n))
        A3 = math.log((k - eta) / (k - eta_n))
        A4 = math.atan(eta/M) - math.atan(eta_n/M)
        de_e_q = self._kappa*k/(3.0*self._r*v_n)*math.log((k-eta_n)/(k-eta))
        de_p_q = (1/theta)*(-k/(M*(k-M))*A1 + k/(M*(k+M))*A2 + 2*k/(k**2-M**2)*A3 - 2/M*A4)
        return (de_e_q + de_p_q)
