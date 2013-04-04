import os
import numpy as np
from daltools import prop
from dalmisc.expect import value

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"

def test_z(): 
    zdiplen, = value('ZDIPLEN', tmpdir=suppdir)
    assert np.allclose(zdiplen, -1.91204774)
    

if __name__ == "__main__":
    setup()
    test_z()
