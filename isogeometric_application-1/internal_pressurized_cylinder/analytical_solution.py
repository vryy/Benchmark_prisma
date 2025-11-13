import math

## REF: http://solidmechanics.org/Text/Chapter4_1/Chapter4_1.php#Sect4_1_9, sec 1.4.9
def get_displacement(x, y, a, b, pa, pb, E, nu):
    r = math.sqrt(x*x + y*y)
    t = math.acos(x/r)
    ur = (1.0+nu)*math.pow(a*b,2)/E/(math.pow(b,2) - math.pow(a,2)) * ((pa-pb)/r + (1.0-2*nu)*(pa*math.pow(a,2)-pb*math.pow(b,2))*r/math.pow(a*b,2))
    ut = 0.0
    ux = ur*math.cos(t) - ut*math.sin(t)
    uy = ur*math.sin(t) + ut*math.cos(t)
    return [ux, uy]

def get_stress_3d(x, y, a, b, pa, pb, E, nu):
    r = math.sqrt(x*x + y*y)
    t = math.acos(x/r)
    c = x/r
    s = y/r
    orr = (pa*math.pow(a, 2) - pb*math.pow(b, 2)) / (math.pow(b, 2) - math.pow(a, 2)) - math.pow(a*b/r, 2) / (math.pow(b, 2) - math.pow(a, 2)) * (pa - pb)
    ott = (pa*math.pow(a, 2) - pb*math.pow(b, 2)) / (math.pow(b, 2) - math.pow(a, 2)) + math.pow(a*b/r, 2) / (math.pow(b, 2) - math.pow(a, 2)) * (pa - pb)
    ort = 0.0
    ozz = 2.0*nu * (pa*math.pow(a, 2) - pb*math.pow(b, 2)) / (math.pow(b, 2) - math.pow(a, 2))
    oxx = orr*math.pow(c, 2) + ott*math.pow(s, 2) - 2.0*ort*s*c
    oyy = orr*math.pow(s, 2) + ott*math.pow(c, 2) + 2.0*ort*s*c
    oxy = (orr - ott)*s*c + ort*(math.pow(c, 2) - math.pow(s, 2))
    return [oxx, oyy, ozz, oxy, 0.0, 0.0]

def get_strain_3d(x, y, a, b, pa, pb, E, nu):
    r = math.sqrt(x*x + y*y)
    t = math.acos(x/r)
    c = x/r
    s = y/r
    ezz = 0.0 # plane strain solution
    err = (1.0+nu)*pow(a*b, 2)/E * (-(pa-pb)/pow(r, 2) + (1.0-2.0*nu)*(pa*pow(a, 2) - pb*pow(b, 2))/pow(a*b, 2)) / (math.pow(b, 2) - math.pow(a, 2)) - nu*ezz
    ett = (1.0+nu)*pow(a*b, 2)/E * ((pa-pb)/pow(r, 2) + (1.0-2.0*nu)*(pa*pow(a, 2) - pb*pow(b, 2))/pow(a*b, 2)) / (math.pow(b, 2) - math.pow(a, 2)) - nu*ezz
    ert = 0.0
    exx = err*math.pow(c, 2) + ett*math.pow(s, 2) - 2.0*ert*s*c
    eyy = err*math.pow(s, 2) + ett*math.pow(c, 2) + 2.0*ert*s*c
    exy = (err - ett)*s*c + ert*(math.pow(c, 2) - math.pow(s, 2))
    return [exx, eyy, ezz, 2.0*exy, 0.0, 0.0]

class Solution:
    def __init__(self, a, b, pa, pb, E, nu):
        self.a = a
        self.b = b
        self.pa = pa
        self.pb = pb
        self.E = E
        self.nu = nu

    def get_displacement(self, x, y, z):
        u = get_displacement(x, y, self.a, self.b, self.pa, self.pb, self.E, self.nu)
        return [u[0], u[1], 0.0]

    def get_stress_3d(self, x, y, z):
        return get_stress_3d(x, y, self.a, self.b, self.pa, self.pb, self.E, self.nu)

    def get_strain_3d(self, x, y, z):
        return get_strain_3d(x, y, self.a, self.b, self.pa, self.pb, self.E, self.nu)


####COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
#nom = 0.0
#denom = 0.0
#for element in model.model_part.Elements:
#    if element.GetValue(IS_INACTIVE) == False:
#        u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model.model_part.ProcessInfo)
#        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
#        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
#        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
#        for i in range(0, len(u)):
#            ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1], a, b, P, 0.0, E, nu)
#            nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2)) * W[i][0] * J0[i][0]
#            denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
#print("Global displacement (L2) error:", math.sqrt(nom / denom))

####COMPUTE GLOBAL STRESS (H1) ERROR###
#nom = 0.0
#denom = 0.0
#for element in model.model_part.Elements:
#    if element.GetValue(IS_INACTIVE) == False:
#        o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model.model_part.ProcessInfo)
#        J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
#        Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
#        W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
#        for i in range(0, len(o)):
#            ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], a, b, P, 0.0, E, nu)
#            nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
#            denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
#print("Global stress (H1) error:", math.sqrt(nom / denom))

