#
# Finite field strength
#
delta = 0.001
#
# Top part of script: setup
#

setup = """
import os
import shutil
import numpy as np
from dalmisc.findif import *
from mol import inp

def assert_(this,ref):
    print this, ref
    assert np.allclose(this, ref, rtol=1e-4)
        

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"
    if os.path.isdir(suppdir):
        shutil.rmtree(suppdir)
    os.mkdir(suppdir)
    os.chdir(suppdir)

def teardown():
    #shutil.rmtree(suppdir)
    pass
"""

