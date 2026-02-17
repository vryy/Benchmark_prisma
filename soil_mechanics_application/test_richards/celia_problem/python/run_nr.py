import numpy as np

# -------------------------------------------------
# Mesh
# -------------------------------------------------
L = 40.0
ne = 40
nn = ne + 1
dz = L / ne
z = np.linspace(0, L, nn)

# -------------------------------------------------
# Time
# -------------------------------------------------
dt = 1.0
nt = 360

# -------------------------------------------------
# Initial & boundary conditions
# -------------------------------------------------
psi_top = -20.7
psi_bot = -61.5
psi = np.full(nn, psi_bot)

# -------------------------------------------------
# Celia (1990) parameters
# -------------------------------------------------
alpha = 1.611e6
theta_s = 0.287
theta_r = 0.075
beta = 3.96

Ks = 0.00944
Ac = 1.175e6
gamma = 4.74

# -------------------------------------------------
# Constitutive laws
# -------------------------------------------------
def theta(p):
    return theta_r + alpha * (theta_s - theta_r) / (alpha + np.abs(p)**beta)

def K(p):
    return Ks * Ac / (Ac + np.abs(p)**gamma)

def capacity(p):
    pabs = np.abs(p)
    return (-alpha * (theta_s - theta_r) * beta *
            pabs**(beta - 1) / (alpha + pabs**beta)**2
            * np.sign(p))

def dK_dpsi(p):
    pabs = np.abs(p)
    return (-Ks * Ac * gamma *
            pabs**(gamma - 1) / (Ac + pabs**gamma)**2
            * np.sign(p))

# -------------------------------------------------
# Time stepping
# -------------------------------------------------
for nstep in range(nt):
    psi_old = psi.copy()
    psi_k = psi.copy()
    psi_k[0] = psi_bot
    psi_k[-1] = psi_top

    for it in range(30):
        J = np.zeros((nn, nn))
        R = np.zeros(nn)

        Ck = capacity(psi_k)
        Kk = K(psi_k)
        dKk = dK_dpsi(psi_k)
        th_k = theta(psi_k)
        th_n = theta(psi_old)

        for e in range(ne):
            i, j = e, e + 1

            if e == ne-1:
                print(psi_k[i])
                print(psi_k[j])
                print(Kk[i])
                print(Kk[j])

            Ke = 0.5 * (Kk[i] + Kk[j])
            dKe_i = 0.5 * dKk[i]
            dKe_j = 0.5 * dKk[j]

            ## use the mid-point value of psi to compute K
            ## it will give different solution
            # psi_mid = 0.5 * (psi_k[i] + psi_k[j])
            # Ke = K(psi_mid)
            # dKe = dK_dpsi(psi_mid)
            # dKe_i = 0.5 * dKe
            # dKe_j = 0.5 * dKe

            grad_psi = (psi_k[j] - psi_k[i]) / dz

            # -------- Residual --------
            R[i] += (dz/2) * (th_k[i] - th_n[i]) / dt
            R[j] += (dz/2) * (th_k[j] - th_n[j]) / dt

            R[i] += -Ke * grad_psi - Ke
            R[j] +=  Ke * grad_psi + Ke

            # -------- Jacobian --------
            # Mass
            J[i,i] += (dz/2) * Ck[i] / dt
            J[j,j] += (dz/2) * Ck[j] / dt

            # Diffusion + gravity (FULL Newton)
            J[i,i] +=  Ke/dz - dKe_i*(grad_psi + 1)
            J[i,j] += -Ke/dz - dKe_j*(grad_psi + 1)

            J[j,i] += -Ke/dz + dKe_i*(grad_psi + 1)
            J[j,j] +=  Ke/dz + dKe_j*(grad_psi + 1)

            if e == ne-1:
                print(R[e:e+2])
                print(J[e:e+2,e:e+2])
                print('------------')

        # Dirichlet BCs
        J[0,:] = 0.0
        J[0,0] = 1.0
        R[0] = psi_k[0] - psi_bot

        J[-1,:] = 0.0
        J[-1,-1] = 1.0
        R[-1] = psi_k[-1] - psi_top

        dpsi = np.linalg.solve(J, -R)
        psi_k += dpsi

        if np.linalg.norm(dpsi, np.inf) < 1e-6:
            break

    psi = psi_k
    print(f"Step {nstep+1:3d}, Newton iterations = {it+1}")

# -------------------------------------------------
# Output
# -------------------------------------------------
print("\nFinal pressure head:")
for zi, pi in zip(z, psi):
    print(f"z = {zi:5.1f} cm : psi = {pi:9.4f} cm")

ifile = open("pressure_head.txt", "w")
ifile.write("z\tp\n")

for zi, pi in zip(z, psi):
    ifile.write("%.16e\t%.16e\n" % (zi, pi))

ifile.close()
