import math
import pathlib

import numpy as np
import numpy.testing as npt
import pytest
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

    @pytest.mark.skip
    def test_h1diag_initial_mo(self):

        cmo = dens.cmo(self.h1, self.S)
        cmo = self.cmo()
        npt.assert_allclose(cmo, ref.C1)

    @pytest.mark.skip
    def test_h1diag_initial_dens(self):

        Da, Db = self.dab()
        npt.assert_allclose(Da, ref.Da)
        npt.assert_allclose(Db, ref.Db)

    @pytest.mark.skip
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
            atol=1e-8
        )

    def test_fco(self):
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()

        npt.assert_allclose(
            np.array(Fc),
            .5*np.array((Fa + Fb)),
            atol=1e-8
        )

    def test_grad_co_ab(self):
        C = self.cmo()
        Dc, Do = self.dco()
        Da, Db = self.dab()
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()
        gco = rohf.grad(self.S, C, Dc, Do, Fc, Fo)
        gab = rohf.grad(self.S, C, Da, Db, Fa, Fb)
        gco = np.array(gco)
        gab = np.array(gab)
        npt.assert_allclose(gco, gab, atol=1e-8)

    def test_gradao_co_ab(self):
        Dc, Do = self.dco()
        Da, Db = self.dab()
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()
        gco = np.array(rohf.gradao(self.S, Dc, Do, Fc, Fo))
        gab = np.array(rohf.gradao(self.S, Da, Db, Fa, Fb))
        npt.assert_allclose(gco, gab, atol=1e-8)

    def test_gradnorm_co_ab(self):
        C = self.cmo()
        Dc, Do = self.dco()
        Da, Db = self.dab()
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()
        gnref = 2.75678
        gnco = rohf.gradnorm(self.S, C, Dc, Do, Fc, Fo, 4, 1)
        gnab = rohf.gradnorm(self.S, C, Da, Db, Fa, Fb, 4, 1)
        assert gnco == approx(gnref)
        assert gnab == approx(gnref)

    def test_gradao_norm_co_ab(self):

        Dc, Do = self.dco()
        Da, Db = self.dab()
        Fc, Fo = self.fco()
        Fa, Fb = self.fab()
        gnref = 2.75678
        gco = rohf.gradao(self.S, Dc, Do, Fc, Fo)
        gab = rohf.gradao(self.S, Da, Db, Fa, Fb)
        g2co = gco & (self.S.I@gco@self.S.I) / 2
        g2ab = gab & (self.S.I@gab@self.S.I) / 2
        assert math.sqrt(g2co) == approx(gnref)
        assert math.sqrt(g2ab) == approx(gnref)

    def test_feffs(self):
        Dc, Do = self.dco()
        Fc, Fo = self.fco()
        Da, Db = self.dab()
        Fa, Fb = self.fab()
        S = self.S
        Feff_co = rohf.jensen(S, Dc, Do, Fc, Fo)
        Feff_ab = S@rohf.Feff(Da@S, Db@S, S.I@Fa, S.I@Fb)
        npt.assert_allclose(
            np.array(Feff_co),
            np.array(Feff_ab),
        )
