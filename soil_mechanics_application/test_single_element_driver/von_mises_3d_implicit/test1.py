################################################
### Single element test for constitutive law ###
### Von Mises
### Hoang-Giang Bui, 2018, Ruhr University Bochum
################################################

import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SoilMechanicsApplication import *

def main(logging=True):

    driver = SingleElementTestDriver()
    con = VonMises3dImplicit()
    econ = LinearElastic3D()
    #print(con)

    ## material parameters
    E        = 1.0e6    # [kPa] Young's modulus
    nu       = 0.25     # [--]  Poisson's coefficient
    c_0      = 10.0     # [kPa] cohesion, initial
    c_inf    = 100.0    # [kPa] cohesion, limit condition at eps_s^p = inf
    rho      = 50.0     # [--] hardening rate
    tol_f    = 1.0e-5*c_0   # [kPa] tolerance on yield condition check
    max_ksub = 1000         # max no. of substeps allowed
    err_tol  = 1.0e-4       # normalized error tolerance for adaptive integration

    material_properties = Properties(1)
    material_properties.SetValue(YOUNG_MODULUS, E)
    material_properties.SetValue(POISSON_RATIO, nu)
    hlaw = ExponentialHardeningLaw(2*c_0, 2*c_inf, rho)
    driver.SetValue(ISOTROPIC_HARDENING_LAW, hlaw, material_properties)

    process_info = ProcessInfo()

    driver.InitializeMaterial(con, material_properties)

    ## define the stress path
    nspb = 1 # number of stress paths
    nsteps = [100]
    DX = [0.1/100]
    strain = Vector(6)
    for i in range(0, 6):
        strain[i] = 0.0

    recorded_values = {}
    recorded_values['EQUIVALENT_DEVIATORIC_STRAIN'] = []
    recorded_values['EQUIVALENT_VOLUMETRIC_STRAIN'] = []
    recorded_values['PRESSURE_P'] = []
    recorded_values['PRESSURE_Q'] = []

    ## initial stress state
    initial_stress = Vector(6)
    initial_stress[0] = 100.0
    initial_stress[1] = 100.0
    initial_stress[2] = 100.0
    initial_stress[3] = 0.0
    initial_stress[4] = 0.0
    initial_stress[5] = 0.0
    initial_elastic_strain = Vector(6)
    econ.ApplyInverse(initial_stress, initial_elastic_strain, E, nu)
    if logging:
        print("initial_elastic_strain:" + str(initial_elastic_strain))
    con.SetValue(STRESSES, initial_stress, process_info)
    con.SetValue(ELASTIC_STRAIN_VECTOR, initial_elastic_strain, process_info)

    ## loop over the stress path
    for ipath in range(0, nspb):
        n = nsteps[ipath]
        dload = DX[ipath]

        ## loop over the steps in the path
        for step in range(0, n):

            # compute load
            strain[0] = strain[0] - 0.5*dload
            strain[1] = strain[1] - 0.5*dload
            strain[2] = strain[2] + dload

            # compute the stress
            driver.CalculateMaterialResponse(con, strain, process_info, material_properties)
            driver.FinalizeSolutionStep(con, process_info, material_properties)
            stress = con.GetValue(THREED_STRESSES)
            if logging:
                print("Stress: " + str(stress))
                print("Step " + str(step+1) + ", path " + str(ipath+1) + " completed")

            recorded_values['EQUIVALENT_VOLUMETRIC_STRAIN'].append(con.GetValue(EQUIVALENT_VOLUMETRIC_STRAIN))
            recorded_values['EQUIVALENT_DEVIATORIC_STRAIN'].append(con.GetValue(EQUIVALENT_DEVIATORIC_STRAIN))
            recorded_values['PRESSURE_P'].append(con.GetValue(PRESSURE_P))
            recorded_values['PRESSURE_Q'].append(con.GetValue(PRESSURE_Q))

    if logging:
        ifile = open("output.txt", "w")
        ifile.write("step\tepv\teps\tp\tq\n")
        for i in range(0, len(recorded_values['EQUIVALENT_DEVIATORIC_STRAIN'])):
            ifile.write(str(i+1) + "\t" + str(recorded_values['EQUIVALENT_VOLUMETRIC_STRAIN'][i]))
            ifile.write("\t" + str(recorded_values['EQUIVALENT_DEVIATORIC_STRAIN'][i]))
            ifile.write("\t" + str(recorded_values['PRESSURE_P'][i]))
            ifile.write("\t" + str(recorded_values['PRESSURE_Q'][i]))
            ifile.write("\n")
        ifile.close()

    return recorded_values

def test():
    recorded_values = main(logging=False)

    last_ep = recorded_values['EQUIVALENT_VOLUMETRIC_STRAIN'][-1]
    last_eq = recorded_values['EQUIVALENT_DEVIATORIC_STRAIN'][-1]
    last_p = recorded_values['PRESSURE_P'][-1]
    last_q = recorded_values['PRESSURE_Q'][-1]
    print("%.15e %.15e %.15e %.15e" % (last_ep, last_eq, last_p, last_q))

    tol = 1e-10
    ref_last_ep = 0.0
    ref_last_eq = 1e-1
    ref_last_p = -1e2
    ref_last_q = 1.987770827054755e+02
    assert(abs(last_ep - ref_last_ep) < tol)
    assert(abs(last_eq - ref_last_eq) < tol)
    assert(abs(last_p - ref_last_p) < tol)
    assert(abs(last_q - ref_last_q) < tol)
    print("Test passed")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main(logging=True)
