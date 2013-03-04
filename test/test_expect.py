import os
import numpy as np
from daltools import prop
from dalmisc.expect import value


def setup():
    thisdir  = os.path.dirname(__file__)
    global suppdir
    suppdir = os.path.join(thisdir, 'test_expect.d')

def test_z(): 
    zdiplen = value('ZDIPLEN', suppdir)
    assert np.allclose(zdiplen, -1.91204774)
    

if __name__ == "__main__":
    setup()
    test_z()
