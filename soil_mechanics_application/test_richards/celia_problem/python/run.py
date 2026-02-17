import numpy as np

# -------------------------------------------------
# Mesh (Celia)
# -------------------------------------------------
L = 40.0           # cm
ne = 40
nn = ne + 1
dz = L / ne
z = np.linspace(0, L, nn)

# -------------------------------------------------
# Time
# -------------------------------------------------
dt = 1.
nt = 1 # 360

# -------------------------------------------------
# Initial & boundary conditions
# -------------------------------------------------
psi_top = -20.7
psi_bot = -61.5
psi = np.full(nn, psi_bot)

# -------------------------------------------------
# Celia (1990) Haverkamp parameters
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
    # Derivative of the Haverkamp theta function w.r.t p
    denom = (alpha + np.abs(p)**beta)
    return (alpha * (theta_s - theta_r) * beta * np.abs(p)**(beta - 1)) / (denom**2)

# # -------------------------------------------------
# # Consistent mass matrix
# # -------------------------------------------------
# M = np.zeros((nn, nn))
# for e in range(ne):
#     i, j = e, e + 1
#     Me = dz / 6 * np.array([[2, 1],
#                             [1, 2]])
#     M[i:i+2, i:i+2] += Me

# -------------------------------------------------
# Time stepping
# -------------------------------------------------
for nstep in range(nt):
    psi_old = psi.copy()
    psi_k = psi.copy()

    for it in range(30):
        Kk = K(psi_k)
        Ck = capacity(psi_k) # You need to define this
        thetak = theta(psi_k)
        theta_old = theta(psi_old)

        # Initialize A and rhs for this iteration
        A = np.zeros((nn, nn))
        rhs = np.zeros(nn)

        for e in range(ne):
            i, j = e, e + 1

            # Use Lumped Mass for both A and RHS for better stability
            Ce_lumped = (dz / 2) * np.array([Ck[i], Ck[j]])

            # Stiffness/Diffusion
            Ke = 0.5 * (Kk[i] + Kk[j])
            Se = (Ke / dz) * np.array([[1, -1], [-1, 1]])

            if e == 0:
                print(psi_k)
                print(Kk[i])
                print(Kk[j])
                print(Ke)
                print(Se)
                print(Ce_lumped)

            # Gravity (Assuming z is positive upward)
            Ge = Ke * np.array([-1, 1])

            if e == 0:
                print(Ge)
                print('--------------')

            # Assembly
            A[i, i] += (Ce_lumped[0] / dt) + Se[0, 0]
            A[i, j] += Se[0, 1]
            A[j, i] += Se[1, 0]
            A[j, j] += (Ce_lumped[1] / dt) + Se[1, 1]

            rhs[i] += (Ce_lumped[0] / dt) * (psi_k[i] - (thetak[i] - theta_old[i])/Ck[i]) + Ge[0]
            rhs[j] += (Ce_lumped[1] / dt) * (psi_k[j] - (thetak[j] - theta_old[j])/Ck[j]) + Ge[1]

        # Dirichlet BCs
        A[0, :] = 0.0
        A[0, 0] = 1.0
        rhs[0] = psi_top

        A[-1, :] = 0.0
        A[-1, -1] = 1.0
        rhs[-1] = psi_bot

        # print(A)
        # print(rhs)
        psi_new = np.linalg.solve(A, rhs)
        # print("diff: %.6e" % (np.linalg.norm(psi_new - psi_k, np.inf)))
        # print(psi_k)
        # print(psi_new)

        if np.linalg.norm(psi_new - psi_k, np.inf) < 1e-6:
            break

        psi_k = psi_new
        # print(psi_k)
        # print(len(psi_k))

    psi = psi_new.copy()
    print(f"Step {nstep+1:3d}, Picard iterations = {it+1}")

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
