import os

import numpy as np
import numpy.testing as npt
import pytest
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
        self.nb = 5

    def test_potnuc(self):
        assert np.allclose(self.EN, 9.055004525638)

    def test_h1diag_initial_mo(self):

        Cmo = dens.cmo(self.h1, self.S)
        Ca = Cmo[:, :self.na]
        Cb = Cmo[:, :self.nb]
        npt.assert_allclose(Ca, ref.Ca)
        npt.assert_allclose(Cb, ref.Cb)

    def test_h1diag_initial_dens(self):

        Cmo = dens.cmo(self.h1, self.S)
        Ca = Cmo[:, :self.na]
        Cb = Cmo[:, :self.nb]
        Da = Ca @ Ca.T
        Db = Cb @ Cb.T
        assert np.allclose(Da, ref.Da)
        assert np.allclose(Db, ref.Db)

    def test_h1diag_initial_fock(self):

        Cmo = dens.cmo(self.h1, self.S)
        Ca = Cmo[:, :self.na]
        Cb = Cmo[:, :self.nb]
        Da = Ca @ Ca.T
        Db = Cb @ Cb.T
        np.allclose(Da, ref.Da)
        np.allclose(Db, ref.Db)

    def test_h1diag_initial_energy(self):
        Cmo = dens.cmo(self.h1, self.S)
        Ca = Cmo[:, :self.na]
        Cb = Cmo[:, :self.nb]
        Da = Ca * Ca.T
        Db = Cb * Cb.T
        (Fa, Fb), = two.fockab((Da, Db), filename=self.aotwoint)
        E = self.EN + rohf.energy(Da, Db, self.h1, Fa, Fb)
        Eref = -73.240064311328
        np.allclose(E, Eref)
