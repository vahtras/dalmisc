import os
from util.full import unit, init
from daltools import rspvec, sirifc, dens, one, prop
from .. import qr
from .common_tests import assert_

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"

    global A, B, C
    A = "XDIPLEN"
    B = "YDIPLEN"
    C = "XDIPLEN"

    AOONEINT = os.path.join(suppdir, 'AOONEINT')
    S = one.read(label="OVERLAP", filename=AOONEINT).unblock().unpack()

    global D, ifc
    SIRIFC = os.path.join(suppdir, 'SIRIFC')
    ifc = sirifc.sirifc(SIRIFC)
    cmo = ifc.cmo.unblock()
    da, db = dens.Dab(ifc_=ifc)
    D = cmo.T*S*(da+db)*S*cmo


    RSPVEC = os.path.join(suppdir, 'RSPVEC')
    global NA, NB, NC
    NA, NB, NC = rspvec.read(A, B, C, propfile=RSPVEC)[0]

    global kA, kB, kC
    kA = rspvec.tomat(NA, ifc, tmpdir = suppdir)
    kB = rspvec.tomat(NB, ifc, tmpdir = suppdir).T
    kC = rspvec.tomat(NC, ifc, tmpdir = suppdir).T

    global a, b, c
    AOPROPER = os.path.join(suppdir, 'AOPROPER')
    #a, b, c = [cmo.T*x*cmo for x in prop.read(A, B, C, filename=AOPROPER, unpack=True)]
    global pmat
    pmat = prop.read(A, B, C, filename=AOPROPER, unpack=True)
    a, b, c = [cmo.T*x*cmo for x in pmat]
    

def test_e3():
    ref = -1.53248530         
    pB = {"kappa": kB}
    pC = {"kappa": kC}
    this = -NA&qr.E3(pB, pC, ifc, tmpdir=suppdir)
    assert_(this, ref)

def test_b2c():
    ref = 7.13089781          
    this = -(kA^(kC^b))&D
    assert_(this, ref)

def test_c2b():
    ref = 6.00528627         
    this = -(kA^(kB^c))&D
    assert_(this, ref)

def test_alt_b2c():
    pB = {"kappa":kB, "matrix":pmat[1]}
    pC = {"kappa":kC, "matrix":pmat[2]}
    ref = 7.13089781 + 6.00528627
    this = -NA&qr.B2C(pB, pC, ifc, tmpdir=suppdir)
    assert_(this, ref)

def test_a2b():
    ref = 3.00264314         
    this = .5*(kC^(kB^a))&D
    assert_(this, ref)

def test_a2c():
    ref = 3.00264314         
    this = .5*(kB^(kC^a))&D
    assert_(this, ref)

def test_alt_a2b():
    pA = {"kappa":kA, "matrix":pmat[0]}
    pB = {"kappa":kB, "matrix":pmat[1]}
    pC = {"kappa":kC, "matrix":pmat[2]}
    pA, pB, pC = [{"kappa":k, "matrix":p} for k, p in zip((kA, kB, kC), pmat)]
    ref = 3.00264314 * 2
    this = (
        -(NB&qr.A2B(pA, pC, ifc, tmpdir=suppdir)) 
        -(NC&qr.A2B(pA, pB, ifc, tmpdir=suppdir))
        )/2
    assert_(this, ref)

if __name__ == "__main__":
    setup()

#   E3  CONTRIBUTION TO HYPVAL         -1.53248530         -1.53248530
# + B2C CONTRIBUTION TO HYPVAL          7.13089781          5.59841250
# + C2B CONTRIBUTION TO HYPVAL          6.00528627         11.60369878
# + A2B CONTRIBUTION TO HYPVAL          3.00264314         14.60634192
# + A2C CONTRIBUTION TO HYPVAL          3.00264314         17.60898505

