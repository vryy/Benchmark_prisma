import math

from KratosMultiphysics import *
from KratosMultiphysics.PlateAndShellApplication import *

def compute_L2_error(model, p, E, nu, h, R):
    sigma = nu / (1.0 - nu)
    # uref1 = p*R**2/(2*(1+sigma)) - p/16*(1.0+3.0*sigma)/(1.0+sigma) + p/(64.0*R**2)
    Rb = R / h
    pb = p*h / (E/(2*(1.0+nu)))
    uref1 = pb/(4*Rb) * ((Rb**2-0.25)**2/Rb + (1.0-2.0*nu)*(Rb**2+0.25)*Rb)
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
            rb = r/h
            # uref = p*(r**2)*(1-nu**2)/(E*h)
            # uref = p*(r**2)*(1-nu**2)/(E*h) * (1+(h/r)/(2*(1-nu)))
            # uref = p*(r**2)*(1-nu**2)/(E*h) * (1+(h/r)*(1-2*nu)/(2*(1-nu)))
            # uref = uref1
            # uref = pb/(4*Rb) * ((Rb**2-0.25)**2/rb + (1.0-2.0*nu)*(Rb**2+0.25)*rb)
            uref = pb*Rb**2/(2*(1.0+sigma))
            # print(" u computed: %.10e" % u[i][2])
            # print(" u ref: %.10e" % uref)
            nom += ((u[i][2] - uref)**2)*w[i][0]*j0[i][0]
            denom += (uref**2)*w[i][0]*j0[i][0]
    return nom/denom

def run(model, logging=True, output=True, p=1e2, drill_stiff=1e2, compute_l2_error=False):
    model.model_part.Properties[1].SetValue(DRILLING_STIFFNESS,            drill_stiff )

    # ## boundary condition
    # #     The middle surface in the cylindrical coordinates r, \theta , x is described by
    # # r = R, -\pi < \theta < \pi, 0 < x < L
    # # At x=0 and x=L the shell edges are free. We impose the following boundary condition on the edges \theta =\pm \pi
    # # u_\theta =0, \psi_\theta =0, \psi_x=0, no restriction for u_x and u_r
    # # The solution is:
    # #       u_1=u_2=0, \psi_1=\psi_2=0, u=p R^2(1-\nu^2)/(Eh),
    # #       \gamma_{11}=0, \gamma_{22}=p R^2(1-\nu^2)/(Eh), \gamma_{12}=0,
    # #       \rho_{11}=\rho_{22}=\rho_{12}=0.
    tol = 1e-4
    # a = 100.0
    # b = 100.0
    l = 1000.0
    # ## BC type 1
    # for node in model.model_part.Nodes:
    #     if (abs(node.Z0) < tol) and (abs(node.Y0 - a) < tol):
    #         node.Fix(DISPLACEMENT_X)
    #         node.Fix(DISPLACEMENT_Z)
    #         node.Fix(ROTATION_X)
    #         node.Fix(ROTATION_Y)
    #         # node.Fix(ROTATION_Z)
    #     if (abs(node.Y0) < tol) and (abs(node.Z0 - b) < tol):
    #         node.Fix(DISPLACEMENT_X)
    #         node.Fix(DISPLACEMENT_Y)
    #         node.Fix(ROTATION_X)
    #         node.Fix(ROTATION_Y)
    #         # node.Fix(ROTATION_Z)

    #     if (abs(node.X0) < tol) or (abs(node.X0-l) < tol):
    #         node.Fix(DISPLACEMENT_X)
    #         # node.Fix(ROTATION_X)
    #         node.Fix(ROTATION_Y)
    #         node.Fix(ROTATION_Z)

    ## BC type 2
    #     u1 = psi1 = 0 at x1 = 0, L.
    #     u2 = 0, psi2 = 0 at x2 = 0, S/2.
    for node in model.model_part.Nodes:
        if abs(node.Z0) < tol:
            node.Fix(DISPLACEMENT_Z)
            node.Fix(ROTATION_Z)
        if abs(node.Y0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.Fix(ROTATION_Y)

        if (abs(node.X0) < tol) or (abs(node.X0 - l) < tol):
            node.Fix(DISPLACEMENT_X)
            node.Fix(ROTATION_X)

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
    #     # node.Fix(DISPLACEMENT_X)
    # # end of debugging

    # # DEBUGGING
    # for node in model.model_part.Nodes:
    #     node.Fix(ROTATION_X)
    #     node.Fix(ROTATION_Y)
    #     node.Fix(ROTATION_Z)
    # # END OF DEBUGGING

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

    if compute_l2_error:
        return model, l2_error
    else:
        return model
