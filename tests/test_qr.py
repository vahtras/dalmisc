import os
from util.full import unit, init
from daltools import rspvec, sirifc, dens, one, prop
from dalmisc import qr

import unittest

class TestQR(unittest.TestCase):

    def setUp(self):
        n, e = os.path.splitext(__file__)
        self.suppdir = n + ".d"

        self.A = "XDIPLEN"
        self.B = "YDIPLEN"
        self.C = "XDIPLEN"

        AOONEINT = os.path.join(self.suppdir, 'AOONEINT')
        S = one.read(label="OVERLAP", filename=AOONEINT).unblock().unpack()

        SIRIFC = os.path.join(self.suppdir, 'SIRIFC')
        self.ifc = sirifc.sirifc(SIRIFC)
        cmo = self.ifc.cmo.unblock()
        da, db = dens.Dab(ifc_=self.ifc)
        self.D = cmo.T*S*(da+db)*S*cmo


        RSPVEC = os.path.join(self.suppdir, 'RSPVEC')
        rspvecs = rspvec.read(self.A, self.B, self.C, propfile=RSPVEC)
        self.NA = rspvecs[(self.A, 0)]
        self.NB = rspvecs[(self.B, 0)]
        self.NC = rspvecs[(self.C, 0)]

        self.kA = rspvec.tomat(self.NA, self.ifc, tmpdir = self.suppdir)
        self.kB = rspvec.tomat(self.NB, self.ifc, tmpdir = self.suppdir).T
        self.kC = rspvec.tomat(self.NC, self.ifc, tmpdir = self.suppdir).T

        AOPROPER = os.path.join(self.suppdir, 'AOPROPER')
        #a, b, c = [cmo.T*x*cmo for x in prop.read(A, B, C, filename=AOPROPER, unpack=True)]
        global pmat
        pmat = prop.read(self.A, self.B, self.C, filename=AOPROPER, unpack=True)
        self.a, self.b, self.c = [cmo.T*x*cmo for x in pmat]
        

    def test_e3(self):
        ref = -1.53248530         
        pB = {"kappa": self.kB}
        pC = {"kappa": self.kC}
        this = -self.NA&qr.E3(pB, pC, self.ifc, tmpdir=self.suppdir)
        self.assertAlmostEqual(this, ref)

    def test_b2c(self):
        ref = 7.13089781          
        this = -(self.kA^(self.kC^self.b))&self.D
        self.assertAlmostEqual(this, ref)

    def test_c2b(self):
        ref = 6.00528627         
        this = -(self.kA^(self.kB^self.c))&self.D
        self.assertAlmostEqual(this, ref)

    def test_alt_b2c(self):
        pB = {"kappa":self.kB, "matrix":pmat[1]}
        pC = {"kappa":self.kC, "matrix":pmat[2]}
        ref = 7.13089781 + 6.00528627
        this = -self.NA&qr.B2C(pB, pC, self.ifc, tmpdir=self.suppdir)
        self.assertAlmostEqual(this, ref)

    def test_a2b(self):
        ref = 3.00264314         
        this = .5*(self.kC^(self.kB^self.a))&self.D
        self.assertAlmostEqual(this, ref)

    def test_a2c(self):
        ref = 3.00264314         
        this = .5*(self.kB^(self.kC^self.a))&self.D
        self.assertAlmostEqual(this, ref)

    def test_alt_a2b(self):
        pA = {"kappa":self.kA, "matrix":pmat[0]}
        pB = {"kappa":self.kB, "matrix":pmat[1]}
        pC = {"kappa":self.kC, "matrix":pmat[2]}
        pA, pB, pC = [{"kappa":k, "matrix":p} for k, p in zip((self.kA, self.kB, self.kC), pmat)]
        ref = 3.00264314 * 2
        this = (
            -(self.NB&qr.A2B(pA, pC, self.ifc, tmpdir=self.suppdir)) 
            -(self.NC&qr.A2B(pA, pB, self.ifc, tmpdir=self.suppdir))
            )/2
        self.assertAlmostEqual(this, ref)

if __name__ == "__main__":
    unittest.main()

#   E3  CONTRIBUTION TO HYPVAL         -1.53248530         -1.53248530
# + B2C CONTRIBUTION TO HYPVAL          7.13089781          5.59841250
# + C2B CONTRIBUTION TO HYPVAL          6.00528627         11.60369878
# + A2B CONTRIBUTION TO HYPVAL          3.00264314         14.60634192
# + A2C CONTRIBUTION TO HYPVAL          3.00264314         17.60898505

