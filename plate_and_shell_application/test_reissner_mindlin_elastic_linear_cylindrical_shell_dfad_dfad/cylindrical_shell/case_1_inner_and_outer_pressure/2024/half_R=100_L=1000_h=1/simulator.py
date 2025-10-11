import math

from KratosMultiphysics import *
from KratosMultiphysics.PlateAndShellApplication import *

def compute_L2_error(model, p, E, nu, h, R):
    sigma = nu / (1.0 - nu)
    Rb = R / h
    pb = p*h / (E/(2*(1.0+nu)))
    for elem in model.model_part.Elements:
        u = elem.CalculateOnIntegrationPoints(LOCAL_PROJECTED_DISPLACEMENT, model.model_part.ProcessInfo)
        # print(u)
        q = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model.model_part.ProcessInfo)
        w = elem.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
        j0 = elem.CalculateOnIntegrationPoints(JACOBIAN_0, model.model_part.ProcessInfo)
        # print(q)
        nom = 0.0
        denom = 0.0
        # print("elem %d --------" % elem.Id)
        for i in range(0, len(q)):
            r = math.sqrt(q[i][1]**2 + q[i][2]**2)
            # uref = p*(r**2)*(1-nu**2)/(E*h)
            # uref = p*(r**2)*(1-nu**2)/(E*h) * (1+(h/r)/(2*(1-nu)))
            # uref = p*(r**2)*(1-nu**2)/(E*h) * (1+(h/r)*(1-2*nu)/(2*(1-nu)))
            # uref = p*(r**2)/(2.0*(1+sigma)) * (1.0 + (1.0-sigma)/(2*r))
            uref = pb*Rb**2/(2*(1.0+sigma))
            # print(" u computed: %.10e" % u[i][2])
            # print(" u ref: %.10e" % uref)
            nom += ((u[i][2] - uref)**2)*w[i][0]*j0[i][0]
            denom += (uref**2)*w[i][0]*j0[i][0]
    return nom/denom

def run(model, logging=True, output=True, p=1e2, compute_l2_error = True):
    ## boundary condition
    #     u1 = psi1 = 0 at x1 = 0, L.
    #     u2 = 0, psi2 = 0 at x2 = 0, S.
    tol = 1e-4
    l = 1000.0
    r = 100.0

    for node in model.model_part.Nodes:
        if (abs(node.X0) < tol) or (abs(node.X0-l) < tol):
            node.Fix(DISPLACEMENT_X)
            node.Fix(ROTATION_X)

        if (abs(node.Z0) < tol):
            # node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            # node.Fix(ROTATION_X)
            node.Fix(ROTATION_Y)

    # # debugging
    # for node in model.model_part.Nodes:
    #     node.Fix(ROTATION_X)
    #     node.Fix(DISPLACEMENT_X)
    # # end of debugging

    # # debugging
    # for node in model.model_part.Nodes:
    #     node.Fix(ROTATION_X)
    #     node.Fix(ROTATION_Y)
    #     node.Fix(ROTATION_Z)
    # # end of debugging

    ## DEBUGGING
    # for node in model.model_part.Nodes:
    #     node.Fix(ROTATION_X)
    #     node.Fix(ROTATION_Y)
    ## END OF DEBUGGING

    # ## body force
    # body_force = Vector(3)
    # body_force[0] = 0
    # body_force[1] = 0
    # body_force[2] = -1.e2
    # model.model_part.Properties[1].SetValue(BODY_FORCE, body_force)

    ## pressure load
    for element in model.model_part.Elements:
        element.SetValue(POSITIVE_FACE_PRESSURE, -0.5*p)
        element.SetValue(NEGATIVE_FACE_PRESSURE, 0.5*p)
        element.SetValue(RADIUS, r)

    ## analysis
    time = 1.0
    model.SolveModel(time)
    if output:
        model.WriteOutput(time)

    # for node in model.model_part.Nodes:
    #     print(node.GetSolutionStepValue(DISPLACEMENT))

    if compute_l2_error:
        E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
        nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
        h = model.model_part.Properties[1].GetValue(THICKNESS)
        R = 100.0
        l2_error = compute_L2_error(model, p, E, nu, h, R)
        print("L2 error: %.10e" % (l2_error))
        return model, l2_error
    else:
        return model
