import os
import numpy as np
from daltools import one, dens
from dalmisc import two


def setup():
    global h1, S, EN, aooneint, aotwoint, na, nb
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_rohf.d')
    aotwoint = os.path.join(suppdir, "AOTWOINT")

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)

def test_basinfo():
    info = two.info(aotwoint)
    assert_(info["nsym"], 1)
    assert_(info["nbas"], [7,0,0,0,0,0,0,0])
    assert_(info["lbuf"], 600)
    assert_(info["nibuf"], 1)
    assert_(info["nbits"], 8)


def test_first_integral():
    for ig, g in two.list_integrals(aotwoint):
        break
    assert_(g, 4.78506540471)
    assert_(ig, [1,1,1,1])
