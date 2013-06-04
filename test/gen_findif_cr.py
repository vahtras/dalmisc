#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_cr.py A B C X file_of_functionals

Checks d<<A; B, C>>/dx(X) = <<A; B, C, D>>
"""

import sys
from common_findif import setup, delta, main, hfweight

file_of_functionals = sys.argv.pop()
A, B, C, X = sys.argv[1:]
targs = ("%s", "%s", "%s", A, B, C, X, delta, A, B, C, X)


#
# Template for functional test calling findif module 
#

template = {}

template["cr_closed_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    qr = FinDif(RspCalc('%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o"], field='%s', delta=%f)).first() 
    cr = RspCalc('%s', '%s', '%s', '%s', wf=wf, dal=dal, mol=inp["h2o"]).exe()
    assert_(qr, cr)
""" % ("%s", "%s", "%s", A, B, C, X, delta, A, B, C, X)


functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

from common_findif import process
process(template, functionals)



