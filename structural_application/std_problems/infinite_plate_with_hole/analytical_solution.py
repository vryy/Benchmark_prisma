import math
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

#REF: Schillinger thesis, G.R.Liu book (example 6.11)
def get_displacement(x, y, P, ri, G, kappa):
    r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    theta = math.acos(x/r)
    c = math.cos(theta)
    s = math.sin(theta)
    c2 = math.cos(2.0*theta)
    s2 = math.sin(2.0*theta)
    aux2 = math.pow(ri/r, 2)
    aux4 = math.pow(ri/r, 4)

    ur = P*r/(4*G) * ( (kappa-1.0)/2 + c2 + aux2 * (1.0+(1.0+kappa)*c2) - aux4*c2 )
    ut = P*r/(4*G) * ( (1.0-kappa) * aux2 - 1.0 - aux4) * s2

    ux = c*ur - s*ut;
    uy = s*ur + c*ut;

    return [ux, uy]

def get_stress_2d(x, y, P, ri, G, kappa):
    r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    theta = math.acos(x/r)
    c2 = math.cos(2.0*theta)
    s2 = math.sin(2.0*theta)
    c4 = math.cos(4.0*theta)
    s4 = math.sin(4.0*theta)
    aux2 = math.pow(ri/r, 2)
    aux4 = math.pow(ri/r, 4)

    o_xx = P * (1.0 - aux2 * (1.5*c2 + c4) + 1.5*aux4*c4)
    o_yy = P * (-aux2 * (0.5*c2 - c4) - 1.5*aux4*c4)
    o_xy = P * (-aux2 * (0.5*s2 + s4) + 1.5*aux4*s4)

    return [o_xx, o_yy, o_xy]

def get_strain_2d(x, y, P, ri, G, kappa):
    r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    theta = math.acos(x/r)
    c = math.cos(theta)
    s = math.sin(theta)
    c2 = math.cos(2.0*theta)
    s2 = math.sin(2.0*theta)
    aux2 = math.pow(ri/r, 2)
    aux4 = math.pow(ri/r, 4)

    ur = P*r/(4*G) * ( (kappa-1.0)/2 + c2 + aux2 * (1.0+(1.0+kappa)*c2) - aux4*c2 )
    d_ur_dr = P/(4*G) * ( (kappa-1.0)/2 + c2 + aux2 * (1.0+(1.0+kappa)*c2) - aux4*c2 ) + P*r/(4*G) * ( -2.0*aux2/r * (1.0+(1.0+kappa)*c2) + 4*aux4*c2/r )
    d_ur_dt = P*r/(4*G) * ( -2.0*s2 + aux2 * (-2.0*(1.0+kappa)*s2) + 2.0*aux4*s2 )
    ut = P*r/(4*G) * ( (1.0-kappa) * aux2 - 1.0 - aux4) * s2
    d_ut_dr = P/(4*G) * ( (1.0-kappa) * aux2 - 1.0 - aux4) * s2 + P*r/(4*G) * ( (1.0-kappa) * (-2.0*aux2/r) + 4.0*aux4/r) * s2
    d_ut_dt = P*r/(4*G) * ( (1.0-kappa) * aux2 - 1.0 - aux4) * 2.0*c2

    # ref for conversion: http://www.public.iastate.edu/~e_m.424/Strain-Cylindrical.pdf
    err = d_ur_dr
    ett = d_ut_dt/r + ur/r
    ert = 0.5 * (d_ur_dt/r + d_ut_dr - ut/r)
    exx = err*math.pow(c, 2) + ett*math.pow(s, 2) - 2.0*ert*s*c
    eyy = err*math.pow(s, 2) + ett*math.pow(c, 2) + 2.0*ert*s*c
    exy = (err - ett)*s*c + ert*(math.pow(c, 2) - math.pow(s, 2))
    return [exx, eyy, 2.0*exy]

def get_stress_3d(x, y, P, ri, G, kappa): # be careful, this is only for plane stress case
    o = get_stress_2d(x, y, P, ri, G, kappa)
    return [o[0], o[1], 0.0, o[2], 0.0, 0.0]

def get_strain_3d(x, y, P, ri, G, kappa, nu):
    e = get_strain_2d(x, y, P, ri, G, kappa)
    ezz = -nu*(e[0] + e[1]) # plane stress solution, set nu = 0 to obtain the plane strain solution
    return [e[0], e[1], ezz, e[2], 0.0, 0.0]

def get_strain_energy_function(P, ri, G, kappa):
    f_x = MonomialFunctionR3R1X()
    f_x2 = MonomialFunctionR3R1X2()
    f_y2 = MonomialFunctionR3R1Y2()
    f_r = PowFunctionR3R1(SumFunctionR3R1(f_x2, f_y2), 0.5)
    f_t = AcosFunctionR3R1(ProductFunctionR3R1(f_x, InverseFunctionR3R1(f_r)))
    f_2t = ScaleFunctionR3R1(2.0, f_t)
    f_4t = ScaleFunctionR3R1(4.0, f_t)
    f_cost = CosFunctionR3R1(f_t)
    f_cos2t = CosFunctionR3R1(f_2t)
    f_cos4t = CosFunctionR3R1(f_4t)
    f_sint = SinFunctionR3R1(f_t)
    f_sin2t = SinFunctionR3R1(f_2t)
    f_sin4t = SinFunctionR3R1(f_4t)
    f_ri_div_r = ScaleFunctionR3R1(ri, InverseFunctionR3R1(f_r))
    f_ri_div_r_pow2 = PowFunctionR3R1(2, f_ri_div_r)
    f_ri_div_r_pow4 = PowFunctionR3R1(4, f_ri_div_r)

    f_aux1 = ScalarFunctionR3R1(0.5*(kappa-1))
    f_aux2 = SumFunctionR3R1(f_aux1, f_cos2t)
    f_aux3 = SumFunctionR3R1(f_aux2, ProductFunctionR3R1(f_ri_div_r_pow2, SumFunctionR3R1( ScalarFunctionR3R1(1.0), ScaleFunctionR3R1(1.0+kappa, f_cos2t) ) ) )
    f_aux4 = SumFunctionR3R1(f_aux3, NegateFunctionR3R1( ProductFunctionR3R1( f_ri_div_r_pow4, f_cos2t) ) )
    f_aux5 = ScaleFunctionR3R1(P/(4*G), f_r)

    f_ur = ProductFunctionR3R1(f_aux5, f_aux4)

    f_aux6 = ScaleFunctionR3R1(1.0 - kappa, f_ri_div_r_pow2)
    f_aux7 = SumFunctionR3R1(f_aux6, ScalarFunctionR3R1(-1.0))
    f_aux8 = SumFunctionR3R1(f_aux7, NegateFunctionR3R1(f_ri_div_r_pow4))
    f_aux81 = ProductFunctionR3R1(f_aux8, f_sin2t)

    f_ut = ProductFunctionR3R1(f_aux5, f_aux81)

    f_ux = SumFunctionR3R1(ProductFunctionR3R1(f_cost, f_ur), NegateFunctionR3R1(ProductFunctionR3R1(f_sint, f_ut)))
    f_uy = SumFunctionR3R1(ProductFunctionR3R1(f_sint, f_ur), ProductFunctionR3R1(f_cost, f_ut))

    f_exx = f_ux.GetDiffFunction(0)
    f_exy = ScaleFunctionR3R1(0.5, SumFunctionR3R1(f_ux.GetDiffFunction(1), f_uy.GetDiffFunction(0)))
    f_eyy = f_uy.GetDiffFunction(1)

    f_aux9 = NegateFunctionR3R1(ProductFunctionR3R1(f_ri_div_r_pow2, SumFunctionR3R1(ScaleFunctionR3R1(1.5, f_cos2t), f_cos4t))) #- (ri/r)^2 * (1.5*cos(2*theta) + cos(4*theta))
    f_aux10 = NegateFunctionR3R1(ProductFunctionR3R1(f_ri_div_r_pow2, SumFunctionR3R1(ScaleFunctionR3R1(0.5, f_cos2t), NegateFunctionR3R1(f_cos4t)))) #-(ri/r)^2 * (0.5*cos(2*theta) - cos(4*theta))
    f_aux11 = NegateFunctionR3R1(ProductFunctionR3R1(f_ri_div_r_pow2, SumFunctionR3R1(ScaleFunctionR3R1(0.5, f_sin2t), f_sin4t))) #-(ri/r)^2 * (0.5*sin(2*theta) + sin(4*theta))
    f_aux12 = ScaleFunctionR3R1(1.5, ProductFunctionR3R1(f_ri_div_r_pow4, f_cos4t))
    f_aux13 = ScaleFunctionR3R1(1.5, ProductFunctionR3R1(f_ri_div_r_pow4, f_sin4t))

    f_oxx = ScaleFunctionR3R1(P, SumFunctionR3R1(ScalarFunctionR3R1(1.0), SumFunctionR3R1(f_aux9, f_aux12)))
    f_oyy = ScaleFunctionR3R1(P, SumFunctionR3R1(f_aux10, NegateFunctionR3R1(f_aux12)))
    f_oxy = ScaleFunctionR3R1(P, SumFunctionR3R1(f_aux11, f_aux13))

    f_aux14 = ProductFunctionR3R1(f_exx, f_oxx)
    f_aux15 = SumFunctionR3R1(f_aux14, ProductFunctionR3R1(f_eyy, f_oyy))
    f_aux16 = SumFunctionR3R1(f_aux15, ScaleFunctionR3R1(2.0, ProductFunctionR3R1(f_exy, f_oxy)))
    f_se = ScaleFunctionR3R1(0.5, f_aux16)
    return f_se

class PlaneStressSolution:
    def __init__(self, P, ri, G, nu):
        self.P = P
        self.ri = ri
        self.G = G
        self.nu = nu
        self.kappa = (3.0-nu)/(1.0+nu)

    def get_displacement(self, x, y, z):
        u = get_displacement(x, y, self.P, self.ri, self.G, self.kappa)
        return [u[0], u[1], 0.0]

    def get_stress_3d(self, x, y, z):
        o = get_stress_2d(x, y, self.P, self.ri, self.G, self.kappa)
        return [o[0], o[1], 0.0, o[2], 0.0, 0.0]

    def get_strain_3d(self, x, y, z):
        e = get_strain_2d(x, y, self.P, self.ri, self.G, self.kappa)
        ezz = -nu*(e[0] + e[1])
        return [e[0], e[1], ezz, e[2], 0.0, 0.0]

class PlaneStrainSolution:
    def __init__(self, P, ri, G, nu):
        self.P = P
        self.ri = ri
        self.G = G
        self.nu = nu
        self.kappa = 3.0-4.0*nu

    def get_displacement(self, x, y, z):
        u = get_displacement(x, y, self.P, self.ri, self.G, self.kappa)
        return [u[0], u[1], 0.0]

    def get_stress_3d(self, x, y, z):
        o = get_stress_2d(x, y, self.P, self.ri, self.G, self.kappa)
        ozz = self.nu*(o[0] + o[1])
        return [o[0], o[1], ozz, o[2], 0.0, 0.0]

    def get_strain_3d(self, x, y, z):
        e = get_strain_2d(x, y, self.P, self.ri, self.G, self.kappa)
        return [e[0], e[1], 0.0, e[2], 0.0, 0.0]

