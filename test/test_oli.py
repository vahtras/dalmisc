import os
from util.full import unit, init, matrix
from dalmisc import oli
from dalmisc.test.test_common import assert_

def setup():
    global suppdir
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_oli.d')

def test_e2n_S():
    refe2 = init([
        [0.78643356, -0.50624296],
        [-0.50624296, 0.78643356]
        ])
    e2 = [oli.e2n(n, tmpdir=suppdir) for n in unit(2)]
    assert_(e2, refe2)

def test_s2n_S():
    refs2 = matrix.diag([2., -2.])
    s2 = [oli.s2n(n, tmpdir=suppdir) for n in unit(2)]
    assert_(s2, refs2)

if __name__ == "__main__":
    setup()
