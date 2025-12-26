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
        self._nu = poisson_ratio
        self._r = 3.0*(1-2*self._nu) / (2*(1+self._nu))

    # Check the yield state of the stress point
    def check_yield(self, p, q, pc, tol = 0):
        return (((q/self._M)**2 + p*(p-pc)) > tol)

    ##############################################################
    #################### DISPLACEMENT CONTROL ####################
    ##############################################################

    # Compute the increment of the volumetric strain in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_elas_de_e_v(self, v_n, p_n, q_n, k, de_a):
        de_v = (3.0*self._r)/(k + self._r) * de_a
        return de_v

    # Compute the increment of the deviatoric strain in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_elas_de_e_q(self, v_n, p_n, q_n, k, de_a):
        de_q = k/(k + self._r) * de_a
        return de_q

    # Compute the hydrostatic pressure in the elastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_elas_p(self, v_n, p_n, q_n, k, de_a):
        de_v = self.get_elas_de_e_v(v_n, p_n, q_n, k, de_a)
        p = p_n*math.exp(v_n/self._kappa * de_v)
        return p

    # Compute the eta in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def compute_eta(self, v_n, p_n, q_n, k, de_a):
        eta_n = q_n/p_n
        theta = v_n/(self._lambda - self._kappa)
        M = self._M
        # print(k)
        # print(M)
        # print(eta_n)
        a1 = -(self._kappa/v_n + 1.0/theta)/3 - self._kappa*k/(3.0*self._r*v_n) + 2.0*k/theta/(k**2-M**2)
        a2 = -k/(theta*M*(k-M))
        a3 = k/(theta*M*(k+M))
        a4 = 1.0/(3.0*theta)
        a5 = -2.0/(theta*M)
        eta = eta_n
        max_iters = 30
        it = 1
        converged = False
        while (converged == False) and (it < max_iters):
            l = a1*math.log((k-eta)/(k-eta_n)) + a2*math.log((M-eta)/(M-eta_n)) + a3*math.log((M+eta)/(M+eta_n)) + a4*math.log((M**2+eta**2)/(M**2+eta_n**2)) + a5*(math.atan(eta/M) - math.atan(eta_n/M)) - de_a
            if abs(l) < 1.0e-10:
                break
            dl = -a1/(k-eta) - a2/(M-eta) + a3/(M+eta) + 2.0*a4*eta/(M**2+eta**2) + M*a5/(M**2+eta**2)
            deta = l/dl
            eta = eta - deta
            if abs(deta) < 1.0e-10:
                converged = True
            it = it + 1
        if ((it == max_iters) and (converged==False)):
            raise Exception("The loop to compute eta does not converge in " + str(max_iters) + " iterations")
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

    # Compute the increment of hydrostatic strain in the plastic domain
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get_de_v(self, v_n, p_n, q_n, k, eta):
        eta_n = q_n/p_n
        theta = v_n/(self._lambda - self._kappa)
        M = self._M
        de_e_v = self._kappa/v_n*math.log((k-eta_n)/(k-eta))
        de_p_v = 1.0/theta*(math.log((k-eta_n)/(k-eta)) + math.log((eta**2+M**2)/(eta_n**2+M**2)))
        return (de_e_v + de_p_v)

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
