import os
from util.full import unit, init
from dalmisc import qr
from dalmisc.test.test_common import assert_
from daltools import rspvec, sirifc, dens, one, prop

def setup():
    global suppdir
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_qr.d')

    global A, B, C
    A = "XDIPLEN"
    B = "YDIPLEN"
    C = "XDIPLEN"

    AOONEINT = os.path.join(suppdir, 'AOONEINT')
    S = one.read(label="OVERLAP", filename=AOONEINT).unblock().unpack()

    global D
    SIRIFC = os.path.join(suppdir, 'SIRIFC')
    ifc = sirifc.sirifc(SIRIFC)
    cmo = ifc.cmo.unblock()
    dc, do = dens.ifc(ifc_=ifc)
    d = dc + do
    D = cmo.T*S*d*S*cmo


    RSPVEC = os.path.join(suppdir, 'RSPVEC')
    global NA, NB, NC
    NA = rspvec.read(A, RSPVEC)
    NB = rspvec.read(B, RSPVEC)
    NC = rspvec.read(C, RSPVEC)
    global kA, kB, kC
    kA = rspvec.tomat(NA, ifc, tmpdir = suppdir)
    kB = rspvec.tomat(NB, ifc, tmpdir = suppdir).T
    kC = rspvec.tomat(NC, ifc, tmpdir = suppdir).T

    global a, b, c
    AOPROPER = os.path.join(suppdir, 'AOPROPER')
    a = cmo.T*prop.read(A, AOPROPER).unpack()*cmo
    b = cmo.T*prop.read(B, AOPROPER).unpack()*cmo
    c = cmo.T*prop.read(C, AOPROPER).unpack()*cmo
    

def test_e3():
    ref = -1.53248530         
    this = -NA&qr.E3(B, C, tmpdir=suppdir)
    assert_(this, ref)

def test_b2c():
    ref = 7.13089781          
    this = -(kA^(kC^b))&D
    assert_(this, ref)

def test_c2b():
    ref = 6.00528627         
    this = -(kA^(kB^c))&D
    assert_(this, ref)

def test_a2b():
    ref = 3.00264314         
    this = .5*(kC^(kB^a))&D
    assert_(this, ref)

def test_b2a():
    ref = 3.00264314         
    this = .5*(kB^(kC^a))&D
    assert_(this, ref)

if __name__ == "__main__":
    setup()

#   E3  CONTRIBUTION TO HYPVAL         -1.53248530         -1.53248530
# + B2C CONTRIBUTION TO HYPVAL          7.13089781          5.59841250
# + C2B CONTRIBUTION TO HYPVAL          6.00528627         11.60369878
# + A2B CONTRIBUTION TO HYPVAL          3.00264314         14.60634192
# + A2C CONTRIBUTION TO HYPVAL          3.00264314         17.60898505

