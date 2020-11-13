import os

import numpy as np
import numpy.testing as npt
from pytest import approx
from daltools import one, dens
import two

from dalmisc import rohf
from . import ref_rohf as ref


class TestROHF:

    def setup(self):
        global suppdir
        n, e = os.path.splitext(__file__)
        suppdir = n + ".d"

        self.aooneint = os.path.join(suppdir, "AOONEINT")
        self.aotwoint = os.path.join(suppdir, "AOTWOINT")
        self.h1 = one.read("ONEHAMIL", self.aooneint).unpack().unblock()
        self.S = one.read("OVERLAP", self.aooneint).unpack().unblock()
        self.EN = one.readhead(self.aooneint)["potnuc"]
        self.na = 5
        self.nb = 4

    def test_potnuc(self):
        assert self.EN == approx(9.263515863231)

    def test_h1diag_initial_mo(self):

        cmo = dens.cmo(self.h1, self.S)
        npt.assert_allclose(cmo, ref.C1)

    def test_h1diag_initial_dens(self):

        Cmo = dens.cmo(self.h1, self.S)
        Ca = Cmo[:, :self.na]
        Cb = Cmo[:, :self.nb]
        Da = Ca @ Ca.T
        Db = Cb @ Cb.T
        npt.assert_allclose(Da, ref.Da)
        npt.assert_allclose(Db, ref.Db)

    def test_h1diag_initial_energy(self):
        Cmo = dens.cmo(self.h1, self.S)
        Ca = Cmo[:, :self.na]
        Cb = Cmo[:, :self.nb]
        Da = Ca * Ca.T
        Db = Cb * Cb.T
        (Fa, Fb), = two.fockab((Da, Db), filename=self.aotwoint)
        E = self.EN + rohf.energy(Da, Db, self.h1, Fa, Fb)
        Eref = -73.4538472272
        assert E == approx(Eref)
