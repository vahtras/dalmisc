import os
import numpy as np
from daltools import one, dens
from dalmisc import twoso


def setup():
    global ao2soint
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_twoso.d')
    ao2soint = os.path.join(suppdir, "AO2SOINT")

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)


def test_first_integral():
    for c, ig, g in twoso.list_integrals(ao2soint):
        break
    assert c == 'x'
    assert_(g, 2.88719908251)
