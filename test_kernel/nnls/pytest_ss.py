import KratosMultiphysics
from KratosMultiphysics import *

def test():
    G_ = [[0.44,   0.50,   0.22,   0.36],
        [0.59,   0.01,   0.23,   0.05],
        [0.76,   0.63,   0.21,   0.91],
        [0.23,   0.98,   0.73,   0.94],
        [0.70,   0.62,   0.26,   0.59],
        [0.67,   0.43,   0.94,   0.85],
        [0.03,   0.09,   0.36,   0.42],
        [0.35,   0.47,   0.99,   0.16],
        [0.57,   0.07,   0.62,   0.93],
        [0.21,   0.76,   0.37,   0.50]]

    b_ = [0.19,
        0.99,
        0.08,
        0.26,
        0.80,
        0.41,
        0.40,
        0.71,
        0.70,
        0.86]

    G = Matrix.FromDoubleList(G_)
    b = Vector.FromList(b_)

    x = Vector(1)

    nnlsSolver = SsNnlsSolver()
    nnlsSolver.Solve(G, x, b, 1e-10, 30, 0)

    print("%.16e %.16e %.16e %.16e" % (x[0], x[1], x[2], x[3]))
    x_ref = [5.9290449022138347e-01, 0.0, 4.2564260703088602e-01, 0.]
    for i in range(4):
        dif = x[i] - x_ref[i]
        print("%.6e" % (dif))
        assert(abs(dif) < 1e-8)

    print("Test passed")

def tag():
    tags = "nnls"
    if not hasattr(KratosMultiphysics, "SsNnlsSolver"):
        tags += ",untested"
    return tags

def print_tag():
    print("Tag(s): " + tag())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        test()
