#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_qr.py A B C file_of_functionals

Checks d<<A; B;>/dx(C) = <<A; B, C>>
"""

import sys
from common_findif import setup, delta, main

A, B, C, file_of_functionals = sys.argv[1:5]


#
# Template for functional test calling findif module 
#

template = {}

template["closed_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    lr = FinDif(LinResp('%s', '%s', wf=wf, dal=dal, mol=inp["h2o"], field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o"]).exe()
    assert_(lr, qr)
""" % ("%s", "%s", "%s", A, B, C, delta, A, B, C)

template["closed_triplet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    lr = FinDif(LinResp('%s', '%s', wf=wf, dal=dal, mol=inp["h2o"], triplet=True, field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o"], triplet=True, aux=".ISPABC\\n 1 1 0").exe()
    assert_(lr, qr)
""" % ("%s", "%s", "%s", A, B, C, delta, A, B, C)

template["open_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    lr = FinDif(LinResp('%s', '%s', wf=wf, dal=dal, mol=inp["h2o+"], field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o+"], parallel=False).exe()
    assert_(lr, qr)
""" % ("%s", "%s", "%s", A, B, C, delta, A, B, C)

template["open_triplet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    lr = FinDif(LinResp('%s', '%s 1', wf=wf, dal=dal, mol=inp["h2o+"], field='%s', delta=%f)).first() 
    qr = QuadResp('%s', '%s 1', '%s', wf=wf, dal=dal, mol=inp["h2o+"], parallel=False).exe()
    assert_(lr, qr)
""" % ("%s", "%s", "%s", A, B, C, delta, A, B, C)

functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

for runtype in template:
    with open("test_findif_qr_" + runtype + ".py", 'w') as runfile:
        runfile.write(setup)
        runfile.write( template[runtype]%('HF', 'HF', 'hf'))
        for f in functionals:
            validfname = f.replace('-', '_').replace('/', '_').replace(' ', '_').replace('*', '')
            dal=validfname.lower()
            if '*' in f: 
                wf = 'DFT\\nGGAKey hf=.5 %s=.5' % f.replace('*', '')
            else:
                wf = 'DFT\\n%s'%f
            runfile.write(template[runtype]%(validfname, wf, dal))
        runfile.write(main)





