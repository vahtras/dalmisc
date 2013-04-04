import os
from util.full import matrix, unit, init
from dalmisc import oli
from dalmisc.test.test_common import assert_

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"

def test_e2n_S():
    refe2 = init([
 [ 0.44049280,     0.04640262,    -0.01404939,    -0.00000000,     0.02644672,    0.02107409],
 [ 0.04640262,     1.02997762,    -0.03772361,     0.02644672,    -0.44005673,   -0.03075434],
 [-0.01404939,    -0.03772361,     0.47089100,     0.02107409,    -0.03075434,    0.00000000],
 [-0.00000000,     0.02644672,     0.02107409,     0.44049280,     0.04640262,   -0.01404939],
 [ 0.02644672,    -0.44005673,    -0.03075434,     0.04640262,     1.02997762,   -0.03772361],
 [ 0.02107409,    -0.03075434,    -0.00000000,    -0.01404939,    -0.03772361,    0.47089100]
        ])
    e2 = [oli.e2n(n, tmpdir=suppdir) for n in unit(6)]
    assert_(e2, refe2)

def test_s2n_S():
    refs2 = matrix.diag([1., 2., 1., -1., -2., -1.])
    s2 = [oli.s2n(n, tmpdir=suppdir) for n in unit(6)]
    assert_(s2, refs2)

if __name__ == "__main__":
    setup()

