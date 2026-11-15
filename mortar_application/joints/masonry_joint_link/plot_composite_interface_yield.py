"""
Composite interface model (Lourenco, 1996 PhD thesis; Lourenco & Rots, 1997)
2D formulation: traction vector t = [sigma, tau]
    sigma : normal stress (tension positive)
    tau   : shear stress

Yield surfaces
    f1 (tension cut-off) : f1 = sigma - sig1(k1)
    f2 (Coulomb friction): f2 = |tau| + sigma*tan(phi) - sig2(k2)
    f3 (compression cap) : f3 = Cnn*sigma^2 + Css*tau^2 + Cn*sigma - sig3(k3)^2

Softening / hardening laws
    sig1 = ft * exp(-ft/GfI * k1)
    sig2 = c  * exp(-c/GfII * k2)
    tan(phi) = tan(phi0) + (tan(phir) - tan(phi0)) * (c - sig2)/c
    sig3 : parabolic hardening -> parabolic softening -> exponential residual

Material values below are illustrative; replace them with your own.
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- parameters
P = dict(
    ft=0.25,       # tensile strength            [MPa]
    GfI=0.018,     # mode I fracture energy      [N/mm]
    c=0.35,        # cohesion                    [MPa]
    tanphi0=0.75,  # initial friction coefficient
    tanphir=0.75,  # residual friction coefficient
    GfII=0.125,    # mode II fracture energy     [N/mm]
    fc=10.5,       # compressive strength        [MPa]
    kp=0.09,       # equivalent plastic relative displacement at peak (cap)
    km=0.49,       # equivalent plastic relative displacement at residual (cap)
    Cnn=1.0,       # cap: normal stress coefficient
    Css=9.0,       # cap: shear stress contribution
    Cn=0.0,        # cap: linear term (centre shift)
)


# ---------------------------------------------------------- hardening laws
def sig1(k1, p=P):
    return p["ft"] * np.exp(-p["ft"] / p["GfI"] * k1)


def sig2(k2, p=P):
    return p["c"] * np.exp(-p["c"] / p["GfII"] * k2)


def tanphi(k2, p=P):
    return p["tanphi0"] + (p["tanphir"] - p["tanphi0"]) * (p["c"] - sig2(k2, p)) / p["c"]


def cap_params(p=P):
    fc = p["fc"]
    si, sp, sm, sr = fc / 3.0, fc, fc / 2.0, fc / 7.0
    kp = p["kp"]
    km = p["km"]
    m = 2.0 * (sm - sp) / (km - kp)          # slope continuity at km
    return si, sp, sm, sr, kp, km, m


def sig3(k3, p=P):
    si, sp, sm, sr, kp, km, m = cap_params(p)
    k = np.asarray(k3, dtype=float)
    a = si + (sp - si) * np.sqrt(np.clip(2 * k / kp - (k / kp) ** 2, 0, None))
    b = sp + (sm - sp) * ((k - kp) / (km - kp)) ** 2
    c = sr + (sm - sr) * np.exp(m * (k - km) / (sm - sr))
    return np.where(k <= kp, a, np.where(k <= km, b, c))


# ---------------------------------------------------------- yield functions
def yield_functions(s, t, k1=0.0, k2=0.0, k3=0.0, p=P):
    f1 = s - sig1(k1, p)
    f2 = np.abs(t) + s * tanphi(k2, p) - sig2(k2, p)
    f3 = p["Cnn"] * s**2 + p["Css"] * t**2 + p["Cn"] * s - sig3(k3, p) ** 2
    return f1, f2, f3


def composite(s, t, k1=0.0, k2=0.0, k3=0.0, p=P):
    """Normalised max of the three functions: zero level = composite surface."""
    f1, f2, f3 = yield_functions(s, t, k1, k2, k3, p)
    return np.maximum.reduce([f1 / p["ft"], f2 / p["c"], f3 / p["fc"] ** 2])


# ---------------------------------------------------------- states to plot
_, _, _, _, KP, KM, _ = cap_params()
STATES = [
    # dict(k1=0.00, k2=0.00, k3=0.0,    label="initial ($\\kappa_i=0$)",        color="k"),
    # dict(k1=0.00, k2=0.00, k3=KP,     label="cap at peak ($\\kappa_3=\\kappa_p$)", color="tab:blue"),
    # dict(k1=0.02, k2=0.10, k3=KM,     label="partly softened",                color="tab:orange"),
    # dict(k1=1.00, k2=1.00, k3=2 * KM, label="residual",                       color="tab:red"),
    # dict(k1=0.000638622, k2=0.255202, k3=0, label="j4d",                       color="tab:red"),
    dict(k1=0.0, k2=2.3143786754206164e-01, k3=0, label="j4d",                       color="tab:red"),
]

# ---------------------------------------------------------- stress points
# Each point is checked against the surface of STATES[state]
STRESS_POINTS = [
    # dict(sigma=-2.0, tau=1.0, state=0, label="A"),
    # dict(sigma=0.10, tau=0.20, state=0, label="B"),
    # dict(sigma=0.30, tau=0.10, state=0, label="C"),
    dict(sigma=0.652023, tau=0.564063, state=0, label="X"),
]

SURFACE_NAMES = ["tension cut-off", "Coulomb friction", "cap"]

def check_point(pt, p=P):
    """Return the values f1, f2, f3 and the list of violated surfaces."""
    st = STATES[pt["state"]]
    f = yield_functions(pt["sigma"], pt["tau"], st["k1"], st["k2"], st["k3"], p)
    f = [float(fi) for fi in f]
    violated = [name for name, fi in zip(SURFACE_NAMES, f) if fi > 0.0]
    return f, violated

def plot_points(ax, s_lim, t_lim):
    for pt in STRESS_POINTS:
        s, t = pt["sigma"], pt["tau"]
        if not (s_lim[0] <= s <= s_lim[1] and t_lim[0] <= t <= t_lim[1]):
            continue
        _, violated = check_point(pt)
        if violated:
            ax.plot(s, t, "X", color="red", ms=10, mec="k", zorder=5)
        else:
            ax.plot(s, t, "o", color="limegreen", ms=8, mec="k", zorder=5)
        ax.annotate(pt["label"], (s, t), xytext=(6, 6),
                    textcoords="offset points", fontweight="bold")

def report_points():
    print("Stress point check")
    for pt in STRESS_POINTS:
        f, violated = check_point(pt)
        st = STATES[pt["state"]]
        status = "inadmissible: " + ", ".join(violated) if violated else "elastic"
        print(f"  {pt['label']}: sigma = {pt['sigma']:7.3f}, tau = {pt['tau']:7.3f} "
              f"| state '{st['label']}' | f1 = {f[0]:8.4f}, f2 = {f[1]:8.4f}, "
              f"f3 = {f[2]:9.4f} -> {status}")

def plot_surfaces(ax, s_lim, t_lim, n=2000, show_individual=True):
    s = np.linspace(*s_lim, n)
    t = np.linspace(*t_lim, n)
    S, T = np.meshgrid(s, t)

    # admissible elastic domain of the initial state
    F0 = composite(S, T, **{k: STATES[0][k] for k in ("k1", "k2", "k3")})
    ax.contourf(S, T, F0, levels=[-1e9, 0], colors=["0.9"])

    # individual surfaces (full, dashed) for the initial state
    if show_individual:
        f1, f2, f3 = yield_functions(S, T)
        for f, col, name in [(f1, "tab:green", "$f_1$ tension cut-off (k1=0)"),
                             (f2, "tab:purple", "$f_2$ Coulomb friction (k2=0)"),
                             (f3, "tab:brown", "$f_3$ cap (k3=0)")]:
            ax.contour(S, T, f, levels=[0], colors=col, linestyles="--", linewidths=0.9)
            ax.plot([], [], ls="--", color=col, lw=0.9, label=name)

    # composite surfaces for each state
    for st in STATES:
        F = composite(S, T, st["k1"], st["k2"], st["k3"])
        ax.contour(S, T, F, levels=[0], colors=st["color"], linewidths=2)
        ax.plot([], [], color=st["color"], lw=2, label=st["label"])

    plot_points(ax, s_lim, t_lim)
    ax.plot([], [], "o", color="limegreen", mec="k", label="stress point (elastic)")
    ax.plot([], [], "X", color="red", mec="k", label="stress point (outside)")

    ax.axhline(0, color="0.5", lw=0.5)
    ax.axvline(0, color="0.5", lw=0.5)
    ax.set_xlabel("$\\sigma$ [MPa]  (tension +)")
    ax.set_ylabel("$\\tau$ [MPa]")
    ax.set_xlim(s_lim)
    ax.set_ylim(t_lim)


# ---------------------------------------------------------- figure 1: surfaces
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

fc = P["fc"]
plot_surfaces(ax1, (-1.15 * fc, 0.15 * fc), (-0.45 * fc, 0.45 * fc))
ax1.set_title("Composite interface yield surface (full view)")
ax1.legend(loc="upper left", fontsize=8)
ax1.set_aspect('equal', 'box')

plot_surfaces(ax2, (-1.5, 0.5), (-1.6, 1.6), show_individual=True)
ax2.set_title("Zoom: tension / shear corner")
ax2.set_aspect('equal', 'box')

fig.tight_layout()

# # ---------------------------------------------------------- figure 2: laws
# fig2, axs = plt.subplots(1, 3, figsize=(14, 4))

# k = np.linspace(0, 0.3, 400)
# axs[0].plot(k, sig1(k))
# axs[0].set_title("Tension softening $\\bar\\sigma_1(\\kappa_1)$")
# axs[0].set_xlabel("$\\kappa_1$ [mm]")

# k = np.linspace(0, 3.0, 400)
# axs[1].plot(k, sig2(k), label="$\\bar\\sigma_2$ (cohesion)")
# axs[1].plot(k, tanphi(k), label="tan $\\phi$")
# axs[1].set_title("Shear softening")
# axs[1].set_xlabel("$\\kappa_2$ [mm]")
# axs[1].legend()

# k = np.linspace(0, 3 * KM, 600)
# axs[2].plot(k, sig3(k))
# for kv, name in [(KP, "$\\kappa_p$"), (KM, "$\\kappa_m$")]:
#     axs[2].axvline(kv, color="0.6", ls=":")
#     axs[2].text(kv, 0.05 * fc, name, ha="left")
# axs[2].set_title("Cap hardening/softening $\\bar\\sigma_3(\\kappa_3)$")
# axs[2].set_xlabel("$\\kappa_3$ [mm]")

# for a in axs:
#     a.set_ylabel("[MPa]")
#     a.grid(alpha=0.3)

# fig2.tight_layout()

plt.show()
