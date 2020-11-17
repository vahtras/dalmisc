import pathlib

import numpy as np
import numpy.testing as npt
from pytest import approx
from daltools import one, dens
import two

from dalmisc import rohf
from . import ref_rohf as ref


class TestROHF:

    def setup(self):
        suppdir = pathlib.Path(__file__).with_suffix('.d')

        self.aooneint = suppdir/"AOONEINT"
        self.aotwoint = suppdir/"AOTWOINT"
        self.h1 = one.read("ONEHAMIL", self.aooneint).unpack().unblock()
        self.S = one.read("OVERLAP", self.aooneint).unpack().unblock()
        self.EN = one.readhead(self.aooneint)["potnuc"]
        self.na = 5
        self.nb = 4

    def test_potnuc(self):
        assert self.EN == approx(9.263515863231)

    def cmo(self):
        return dens.cmo(self.h1, self.S)

    def cab(self):
        cmo = self.cmo()
        ca = cmo[:, :self.na]
        cb = cmo[:, :self.nb]
        return ca, cb

    def dab(self):
        ca, cb = self.cab()
        da = ca @ ca.T
        db = cb @ cb.T
        return da, db

    def fab(self):
        da, db = self.dab()
        (fa, fb), = two.fockab((da, db), filename=self.aotwoint)
        return fa, fb

    def test_h1diag_initial_mo(self):

        cmo = dens.cmo(self.h1, self.S)
        cmo = self.cmo()
        npt.assert_allclose(cmo, ref.C1)

    def test_h1diag_initial_dens(self):

        Da, Db = self.dab()
        npt.assert_allclose(Da, ref.Da)
        npt.assert_allclose(Db, ref.Db)

    def test_h1diag_initial_energy(self):

        Da, Db = self.dab()
        (Fa, Fb), = two.fockab((Da, Db), filename=self.aotwoint)
        E = self.EN + rohf.energy(Da, Db, self.h1, Fa, Fb)
        Eref = -73.4538472272
        assert E == approx(Eref)

    def dco(self):
        da, db = self.dab()
        dc = 2*db
        do = da - db
        return dc, do

    def fco(self):
        Dc, Do = self.dco()
        Fc = two.fock(Dc + Do, filename=self.aotwoint)
        Fo = two.fock(Do, hfc=0, filename=self.aotwoint) + Fc
        return Fc, Fo

    def test_fo(self):
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()

        npt.assert_allclose(
            np.array(Fo - Fc),
            .5*np.array((Fa - Fb)),
        )

    def test_fco(self):
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()

        npt.assert_allclose(
             np.array(Fc),
             .5*np.array((Fa + Fb)),
        )
