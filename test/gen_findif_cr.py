#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_cr.py A B C D file_of_functionals

Checks d<<A; B, C>>/dx(D) = <<A; B, C, D>>
"""

import sys
from common_findif import setup, delta, main

A, B, C, D, file_of_functionals = sys.argv[1:6]


#
# Template for functional test calling findif module 
#

template = {}

template["closed_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    qr = FinDif(QuadResp('%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o"], field='%s', delta=%f)).first() 
    cr = CubResp('%s', '%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o"]).exe()
    assert_(qr, cr)
""" % ("%s", "%s", "%s", A, B, C, D, delta, A, B, C, D)


functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

for runtype in template:
    with open("test_findif_cr_" + runtype + ".py", 'w') as runfile:
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





