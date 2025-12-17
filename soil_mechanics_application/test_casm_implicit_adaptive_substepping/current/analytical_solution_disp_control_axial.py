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

    # Compute the yield function
    def compute_f(self, p, q, pc):
        f = (q/(self._M*p))**self._N + math.log(p/pc) / math.log(self._R)
        return f

    # Check the yield state of the stress point
    def check_yield(self, p, q, pc, tol = 0):
        f = self.compute_f(p, q, pc)
        return (f > tol)

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
        N = self._N
        R = self._R
        # print(k)
        # print(M)
        # print(eta_n)
        A = (2.0*(M**2) - 3.0*(3.0+M))/(M-k)
        B = (3.0*(3.0+M) - 2.0*k*M)/(M-k)
        a1 = -(self._kappa/v_n + (1.0+B/3)/theta)/3 - k*self._kappa/(3.0*self._r*v_n)
        a2 = -A/(9.0*theta)
        a3 = ((k/M)**N)*N*math.log(R)*(1.0+B/3)/(3.0*theta)
        a4 = (k/M)*N*math.log(R)*A/(9.0*theta)
        a5 = -a3
        a6 = -N*math.log(R)*A/(9.0*theta)
        eta = eta_n
        max_iters = 30
        it = 1
        converged = False
        # print("eta_n: " + str(eta))
        while (converged == False) and (it < max_iters):
            FN1_k = F(N-1, eta_n/k, eta/k)
            FN1_M = F(N-1, eta_n/M, eta/M)
            FN_k = F(N, eta_n/k, eta/k)
            FN_M = F(N, eta_n/M, eta/M)
            l = a1*math.log((k-eta)/(k-eta_n)) + a2*math.log((M-eta)/(M-eta_n)) + a3*FN1_k + a4*FN1_M + a5*FN_k + a6*FN_M - de_a
            # diff_l = self.get_de_v(v_n, p_n, q_n, k, eta)/3 + self.get_de_q(v_n, p_n, q_n, k, eta) - (a1*math.log((k-eta)/(k-eta_n)) + a2*math.log((M-eta)/(M-eta_n)) + a3*FN1_k + a4*FN1_M + a5*FN_k + a6*FN_M)
            # print("diff_l: " + str(diff_l))
            if abs(l) < 1.0e-10:
                break
            dl = -a1/(k-eta) - a2/(M-eta) + a3*((eta/k)**(N-1))/(k-eta) + a4*((eta/M)**(N-1))/(M-eta) + a5*((eta/k)**N)/(k-eta) + a6*((eta/M)**N)/(M-eta)
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
        pc = p * math.exp(math.log(self._R)*(eta/self._M)**self._N)
        return pc

    # Compute the preconsolidation pressure in the plastic domain (second method)
    # Input:
    #   v_n     specific volume of the previous step
    #   p_n     hydrostatic pressure of the previous step
    #   q_n     deviatoric pressure of the previous step
    #   k       proportional loading ratio
    #   de_v    increment of volumetric strain
    def get__pc(self, v_n, p_n, q_n, k, eta, pc_n):
        eta_n = q_n/p_n
        theta = v_n/(self._lambda - self._kappa)
        M = self._M
        N = self._N
        R = self._R
        FN1 = F(N-1, eta_n/k, eta/k)
        FN = F(N, eta_n/k, eta/k)
        de_p_v = 1.0/theta*(math.log((k-eta_n)/(k-eta)) + ((k/M)**N)*N*math.log(R)*(FN1 - FN))
        pc = pc_n*math.exp(theta*de_p_v)
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
        N = self._N
        R = self._R
        de_e_v = self._kappa/v_n*math.log((k-eta_n)/(k-eta))
        FN1 = F(N-1, eta_n/k, eta/k)
        FN = F(N, eta_n/k, eta/k)
        de_p_v = 1.0/theta*(math.log((k-eta_n)/(k-eta)) + ((k/M)**N)*N*math.log(R)*(FN1 - FN))
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
        N = self._N
        R = self._R
        A = (2.0*(M**2) - 3.0*(3.0+M))/(M-k)
        B = (3.0*(3.0+M) - 2.0*k*M)/(M-k)
        FN1_k = F(N-1, eta_n/k, eta/k)
        FN_k = F(N, eta_n/k, eta/k)
        FN1_M = F(N-1, eta_n/M, eta/M)
        FN_M = F(N, eta_n/M, eta/M)
        p = self.get_p(v_n, p_n, q_n, k, eta)
        de_e_q = (self._kappa*k/(3*self._r*v_n))*math.log(p/p_n)
        de_p_v = (1.0/theta)*(math.log(p/p_n) + ((k/M)**N)*N*math.log(R)*(FN1_k - FN_k))
        de_p_q = A/(9.0*theta)*(math.log((M-eta_n)/(M-eta)) + N*math.log(R)*(k/M*FN1_M - FN_M)) + B/9.0*de_p_v
        return (de_e_q + de_p_q)

    def get_dphi(self, v_n, p_n, q_n, k, eta):
        eta_n = q_n/p_n
        theta = v_n/(self._lambda - self._kappa)
        M = self._M
        N = self._N
        R = self._R
        FN1 = F(N-1, eta_n/M, eta/M)
        FN = F(N, eta_n/M, eta/M)
        p = self.get_p(v_n, p_n, q_n, k, eta)
        q = q_n + k*(p-p_n)
        de_p_v = (1.0/theta)*(math.log(p/p_n) + ((k/M)**N)*N*math.log(R)*(FN1 - FN))
        dgdp = 3.0*(3.0+2.0*M)/(3.0*p+2.0*q) - 3.0*(3.0-M)/(3.0*p-q)
        dphi = de_p_v / dgdp
        return dphi
