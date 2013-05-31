#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_lr.py A B file_of_functionals

Checks d<A>/dx(B) = <<A; B>>
"""

import sys
from common_findif import setup, delta, main

A, B, file_of_functionals = sys.argv[1:4]


#
# Template for functional test calling findif module 
#

template = {}

template["closed_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    ev = FinDif(ExpVal('%s', wf=wf, dal=dal, mol=inp["h2o"], field='%s', delta=%f)).first() 
    lr = LinResp('%s', '%s', wf=wf, dal=dal, mol=inp["h2o"]).exe()
    assert_(ev, lr)
""" % ("%s", "%s", "%s", A, B, delta, A, B)

template["open_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    ev = FinDif(ExpVal('%s', wf=wf, dal=dal, mol=inp["h2o+"], field='%s', delta=%f)).first() 
    lr = LinResp('%s', '%s', wf=wf, dal=dal, mol=inp["h2o+"]).exe()
    assert_(ev, lr)
""" % ("%s", "%s", "%s", A, B, delta, A, B)

template["open_triplet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    ev = FinDif(ExpVal('%s', wf=wf, dal=dal, mol=inp["h2o+"], triplet=True, field='%s', delta=%f)).first() 
    lr = LinResp('%s 1', '%s', wf=wf, dal=dal, mol=inp["h2o+"], triplet=False).exe()
    assert_(ev, lr)
""" % ("%s", "%s", "%s", A, B, delta, A, B)

functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

for runtype in template:
    with open("test_findif_lr_" + runtype + ".py", 'w') as runfile:
        runfile.write(setup)
        runfile.write( template[runtype]%('HF', 'HF', 'hf'))
        for f in functionals:
            validfname = f.replace('-', '_').replace('/', '_').replace(' ', '_').replace('*', '')
            dal=validfname.lower()
            if '*' in f: 
                wf = 'DFT\\nGGAKey hf=.1 %s=.9' % f.replace('*', '')
            else:
                wf = 'DFT\\n%s'%f
            runfile.write(template[runtype]%(validfname, wf, dal))
        runfile.write(main)





