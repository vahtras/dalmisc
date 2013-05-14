#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_qr.py A B C file_of_functionals

Checks d<<A; B;>/dx(C) = <<A; B, C>>
"""

import sys
from common_findif import setup, delta

A, B, C, file_of_functionals = sys.argv[1:5]


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
    lr = FinDif(LinResp('%s', '%s', wf=wf, mol=inp["h2o"], field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s', '%s', wf=wf, mol=inp["h2o"]).exe()
    assert_(lr, qr)
""" % ("%s", "%s", A, B, C, delta, A, B, C)

template["closed_triplet"] = """
def test_LRx_QR_%s():
    wf='%s'
    lr = FinDif(LinResp('%s', '%s', wf=wf, mol=inp["h2o"], triplet=True, field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s', '%s', wf=wf, mol=inp["h2o"], triplet=True, aux=".ISPABC\\n 1 1 0").exe()
    assert_(lr, qr)
""" % ("%s", "%s", A, B, C, delta, A, B, C)

template["open_singlet"] = """
def test_LRx_QR_%s():
    wf='%s'
    lr = FinDif(LinResp('%s', '%s', wf=wf, mol=inp["h2o+"], field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s', '%s', wf=wf, mol=inp["h2o+"], parallel=False).exe()
    assert_(lr, qr)
""" % ("%s", "%s", A, B, C, delta, A, B, C)

template["open_triplet"] = """
def test_LRx_QR_%s():
    wf='%s'
    lr = FinDif(LinResp('%s', '%s 1', wf=wf, mol=inp["h2o+"], field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s 1', '%s', wf=wf, mol=inp["h2o+"], parallel=False).exe()
    assert_(lr, qr)
""" % ("%s", "%s", A, B, C, delta, A, B, C)

functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

for runtype in template:
    with open("test_findif_qr_" + runtype + ".py", 'w') as runfile:
        runfile.write(setup)
        runfile.write( template[runtype]%('HF', 'HF'))
        for f in functionals:
            validfname = f.replace('-', '_').replace('/', '_')
            wf = 'DFT\\n%s'%f
            runfile.write(template[runtype]%(validfname, wf))
        runfile.write(main)





