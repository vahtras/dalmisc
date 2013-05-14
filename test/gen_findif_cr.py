#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_cr.py A B C D file_of_functionals

Checks d<<A; B, C>>/dx(D) = <<A; B, C, D>>
"""

import sys
A, B, C, D, file_of_functionals = sys.argv[1:6]

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

#
# Bottom part of script: main (to invoke inividual tests)
#

main = """
if __name__ == "__main__":
    import sys
    setup()
    eval("test_LRx_QR_%s()"%sys.argv[1])
    teardown()
"""

#
# Template for functional test calling findif module 
#

template = {}

template["closed_singlet"] = """
def test_LRx_QR_%s():
    wf='%s'
    qr = FinDif(QuadResp('%s', '%s', '%s', wf=wf, mol=inp["h2o"], field='%s', delta=0.001)).first() 
    cr = CubResp('%s', '%s', '%s', '%s', wf=wf, mol=inp["h2o"]).exe()
    assert_(qr, cr)
""" % ("%s", "%s", A, B, C, D, A, B, C, D)


functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

for runtype in template:
    with open("test_findif_cr_" + runtype + ".py", 'w') as runfile:
        runfile.write(setup)
        runfile.write( template[runtype]%('HF', 'HF'))
        for f in functionals:
            validfname = f.replace('-', '_').replace('/', '_')
            wf = 'DFT\\n%s'%f
            runfile.write(template[runtype]%(validfname, wf))
        runfile.write(main)





