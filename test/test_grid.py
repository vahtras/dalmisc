import os
import numpy as np
from dalmisc import grid

def assert_(this,ref):
    print this
    print ref
    assert np.allclose(this, ref)

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"
    os.chdir(suppdir)

def teardown():
    pass

def test_readfirst():
    wRref = [
        -8.575690271798617, 
        -4.223255559902539, 
        -2.8653033453872707, 
        3.3506109330112563,
        ]
    wR = grid.readfirst('DALTON.QUAD')
    assert_(wR, wRref)


