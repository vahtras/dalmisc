import os
from util.full import matrix, unit, init
from dalmisc import oli
from .common_tests import assert_

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"

def test_e2n_S():
    refe2 = init([
       [0.91587038, -0.00000000 ],
       [0.00000000,  0.91587038 ]
        ])
    e2 = init([oli.e2n(n, tmpdir=suppdir) for n in unit(2)])
    assert_(e2, refe2)

def test_s2n_S():
    refs2 = matrix.diag([1, -1])
    s2 = init([oli.s2n(n, tmpdir=suppdir) for n in unit(2)])
    assert_(s2, refs2)

if __name__ == "__main__":
    setup()

