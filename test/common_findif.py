#
# Finite field strength
#
delta = 0.001
rtol=1e-3
#
# Top part of script: setup
#

setup = """
import os
import shutil
import numpy as np
from dalmisc.findif import *
from mol import inp

def assert_(num,ana):
    atol = 1e-8
    print "Numerical ", num
    print "Analytical", ana
    print "Difference", abs(num-ana)
    print "Target difference", atol + %f*abs(ana)
    assert np.allclose(num, ana, rtol=%f)
        

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
""" % (rtol, rtol)

#
# Bottom part of script: main (to invoke inividual tests)
#

main = """
if __name__ == "__main__":
    import sys
    setup()
    eval("test_findif_%s()"%sys.argv[1])
    teardown()
"""
