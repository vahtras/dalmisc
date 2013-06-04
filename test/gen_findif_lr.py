#!/usr/bin/env python
"""Generate a set of test scipts to pass through nosetests

Tests finite field tests of expectation values

Usage:
./gen_findif_lr.py A X file_of_functionals

Checks d<A>/dx(X) = <<A; X>>
"""

import sys
from common_findif import setup, delta, main, hfweight

file_of_functionals = sys.argv.pop()
A, X = sys.argv[1:]


#
# Template for functional test calling findif module 
#

template = {}

template["lr_closed_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    ev = FinDif(RspCalc('%s', wf=wf, dal=dal, mol=inp["h2o"], field='%s', delta=%f)).first() 
    lr = RspCalc('%s', '%s', wf=wf, dal=dal, mol=inp["h2o"]).exe()
    assert_(ev, lr)
""" % ("%s", "%s", "%s", A, X, delta, A, X)

template["lr_open_singlet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    ev = FinDif(RspCalc('%s', wf=wf, dal=dal, mol=inp["h2o+"], field='%s', delta=%f)).first() 
    lr = RspCalc('%s', '%s', wf=wf, dal=dal, mol=inp["h2o+"]).exe()
    assert_(ev, lr)
""" % ("%s", "%s", "%s", A, X, delta, A, X)

template["lr_open_triplet"] = """
def test_findif_%s():
    wf='%s'
    dal='%s'
    ev = FinDif(RspCalc('%s', wf=wf, dal=dal, mol=inp["h2o+"], triplet=True, field='%s', delta=%f)).first() 
    lr = RspCalc('%s 1', '%s', wf=wf, dal=dal, mol=inp["h2o+"], triplet=False).exe()
    assert_(ev, lr)
""" % ("%s", "%s", "%s", A, X, delta, A, X)

functionals = [ line.strip() for line in open(file_of_functionals) ] 

#
# Process all runtypes and functionals defined in input file
#

from common_findif import process
process(template, functionals)
