try:
    from KratosMultiphysics import *
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.P4estApplication import *
    kratos_all_modules_are_imported_successfully = True
except:
    kratos_all_modules_are_imported_successfully = False

def main():
    # specify the order of the discretization
    p4est_order = P8estSecondOrder()

    # create a sample quad for the p4est tree
    p4est_quad_nodal = P8estQuadData(p4est_order.Value)

    p4est_quad_int = P8estQuadData(p4est_order.Value)

    p4est_quad = P8estQuad(p4est_quad_nodal, p4est_quad_int)

    # create the P8est model
    p4est_util = P8estUtilities()
    p4est_model = p4est_util.CreateModel(p4est_order, mpi.world, p4est_quad)

    p4est_model.BeginModelPart()
    p4est_model.EndModelPart()

    mp = p4est_model.GetModelPart()
    # print(mp.GetCommunicator())

    a = 10*(mpi.rank+1)
    sa = mp.GetCommunicator().SumAllInt(a)
    assert(sa == 5*mpi.size*(mpi.size+1))
    b = 1.5*(mpi.rank+1)
    sb = mp.GetCommunicator().SumAllScalar(b)
    assert(sb == 0.75*mpi.size*(mpi.size+1))

    ma = mp.GetCommunicator().MinAllInt(a)
    assert(ma == 10)
    mb = mp.GetCommunicator().MinAllScalar(b)
    assert(mb == 1.5)

    Ma = mp.GetCommunicator().MaxAllInt(a)
    assert(Ma == 10*mpi.size)
    Mb = mp.GetCommunicator().MaxAllScalar(b)
    assert(Mb == 1.5*mpi.size)

    #############

    p4est_quad_nodal = ComplexP8estQuadData(p4est_order.Value)

    p4est_quad_int = ComplexP8estQuadData(p4est_order.Value)

    p4est_quad = ComplexP8estQuad(p4est_quad_nodal, p4est_quad_int)

    # create the complex P8est model
    p4est_util = ComplexP8estUtilities()
    p4est_complex_model = p4est_util.CreateModel(p4est_order, mpi.world, p4est_quad)

    p4est_complex_model.BeginModelPart()
    p4est_complex_model.EndModelPart()

    mp = p4est_complex_model.GetModelPart()
    # print(mp)

    b = complex(1.5,1.0)*(mpi.rank+1)
    sb = mp.GetCommunicator().SumAllScalar(b)
    assert(sb == complex(0.75,0.5)*mpi.size*(mpi.size+1))

def test():
    main()
    print("Test passed")

def tag():
    if kratos_all_modules_are_imported_successfully:
        return "p4est"
    else:
        return "untested"

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
