import os
import numpy as np
from daltools import one, dens
from util.full import init
from dalmisc import twoso



def setup():
    global ao2soint
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_twoso.d')
    ao2soint = os.path.join(suppdir, "AO2SOINT")

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)


def test_first_integral():
    for c, ig, g in twoso.list_integrals(ao2soint, label="AO2SOINT"):
        break
    assert c == 0
    assert_(g, 2.88719908251)

dc = init([
       [ 2.10662148,    -0.44832042,     0.00000000,     0.10863721,     0.00000000,     -0.02702924,    -0.02702924],
       [-0.44832042,     1.97509195,    -0.00000000,    -0.61381008,     0.00000000,     -0.03765992,    -0.03765992],
       [ 0.00000000,    -0.00000000,     0.73437339,     0.00000000,     0.00000000,      0.54049192,    -0.54049192],
       [ 0.10863721,    -0.61381008,     0.00000000,     1.22749198,     0.00000000,      0.47809924,     0.47809924],
       [ 0.00000000,     0.00000000,     0.00000000,     0.00000000,     2.00000000,      0.00000000,     0.00000000],
       [-0.02702924,    -0.03765992,     0.54049192,     0.47809924,     0.00000000,      0.60856871,    -0.18702525],
       [-0.02702924,    -0.03765992,    -0.54049192,     0.47809924,     0.00000000,     -0.18702525,     0.60856871] 
    ])


fcref = -init([
        [ 0.00000000,     0.00000000,     0.37936213,    -0.00000000,     0.00000000,      0.00253951,    -0.00253951],
        [-0.00000000,     0.00000000,     0.10878631,    -0.00000000,     0.00000000,     -0.01699688,     0.01699688],
        [-0.37936213,    -0.10878631,     0.00000000,   -11.60434639,     0.00000000,     -0.98565692,    -0.98565692],
        [ 0.00000000,     0.00000000,    11.60434639,     0.00000000,     0.00000000,      1.20385366,    -1.20385366],
        [ 0.00000000,     0.00000000,     0.00000000,     0.00000000,     0.00000000,      0.00000000,     0.00000000],
        [-0.00253951,     0.01699688,     0.98565692,    -1.20385366,     0.00000000,      0.00000000,    -0.36840521],
        [ 0.00253951,    -0.01699688,     0.98565692,     1.20385366,     0.00000000,      0.36840521,     0.00000000] 
    ])

def test_fc():
    fc = twoso.fock(dc, 'z' , filename=ao2soint, label='AO2SOINT')
    assert_(fc, fcref)

def test_fab():
    da = db = 0.5*dc
    fa, fb  = twoso.fockab(da, db, 'z' , filename=ao2soint, label='AO2SOINT')
    fc = 0.5*(fa - fb)
    assert_(fc, fcref)
