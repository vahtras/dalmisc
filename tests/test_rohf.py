import math
import pathlib

import numpy as np
import numpy.testing as npt
import pytest
from pytest import approx, fixture
from daltools import one, dens
import two

from dalmisc import rohf
from . import ref_rohf as ref

@fixture
def suppdir():
    return pathlib.Path(__file__).with_suffix('.d')

@fixture
def aooneint(suppdir):
    return suppdir / "AOONEINT"

@fixture
def aotwoint(suppdir):
    return suppdir / "AOTWOINT"

@fixture
def h1(aooneint):
    return one.read("ONEHAMIL", aooneint).unpack().unblock()

@fixture
def S(aooneint):
    return one.read("OVERLAP", aooneint).unpack().unblock()

class TestROHF:

    na = 5
    nb = 4

        #self.na = 5
        #self.nb = 4

    def test_potnuc(self, aooneint):
        assert one.readhead(aooneint)["potnuc"] == approx(9.263515863231)

    def cmo(self, h1, S):
        return dens.cmo(h1, S)

    def cab(self, h1, S):
        cmo = self.cmo(h1, S)
        ca = cmo[:, :self.na]
        cb = cmo[:, :self.nb]
        return ca, cb

    def dab(self, h1, S):
        ca, cb = self.cab(h1, S)
        da = ca @ ca.T
        db = cb @ cb.T
        return da, db

    def fab(self, h1, S, aotwoint):
        da, db = self.dab(h1, S)
        (fa, fb), = two.fockab((da, db), filename=aotwoint)
        return fa, fb

    @pytest.mark.skip
    def test_h1diag_initial_mo(self, h1, S):

        cmo = dens.cmo(h1, S)
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

    def dco(self, h1, S):
        da, db = self.dab(h1, S)
        dc = 2*db
        do = da - db
        return dc, do

    def fco(self, h1, S, aotwoint):
        Dc, Do = self.dco(h1, S)
        Fc = two.fock(Dc + Do, filename=aotwoint)
        Fo = two.fock(Do, hfc=0, filename=aotwoint) + Fc
        return Fc, Fo

    def test_fo(self, h1, S, aotwoint):
        Fc, Fo = self.fco(h1, S, aotwoint)
        Fa, Fb = self.fab(h1, S, aotwoint)

        npt.assert_allclose(
            np.array(Fo - Fc),
            .5*np.array((Fa - Fb)),
            atol=1e-8
        )

    def test_fco(self, h1, S, aotwoint):
        Fc, Fo = self.fco(h1, S, aotwoint)
        Fa, Fb = self.fab(h1, S, aotwoint)

        npt.assert_allclose(
            np.array(Fc),
            .5*np.array((Fa + Fb)),
            atol=1e-8
        )

    def test_grad_co_ab(self, h1, S, aotwoint):
        C = self.cmo(h1, S)
        Dc, Do = self.dco(h1, S)
        Da, Db = self.dab(h1, S)
        Fc, Fo = self.fco(h1, S, aotwoint)
        Fa, Fb = self.fab(h1, S, aotwoint)
        gco = rohf.grad(S, C, Dc, Do, Fc, Fo)
        gab = rohf.grad(S, C, Da, Db, Fa, Fb)
        gco = np.array(gco)
        gab = np.array(gab)
        npt.assert_allclose(gco, gab, atol=1e-8)

    def test_gradao_co_ab(self, h1, S, aotwoint):
        Dc, Do = self.dco(h1, S)
        Da, Db = self.dab(h1, S)
        Fc, Fo = self.fco(h1, S, aotwoint)
        Fa, Fb = self.fab(h1, S, aotwoint)
        gco = np.array(rohf.gradao(S, Dc, Do, Fc, Fo))
        gab = np.array(rohf.gradao(S, Da, Db, Fa, Fb))
        npt.assert_allclose(gco, gab, atol=1e-8)

    def test_gradnorm_co_ab(self, h1, S, aotwoint):
        C = self.cmo(h1, S)
        Dc, Do = self.dco(h1, S)
        Da, Db = self.dab(h1, S)
        Fc, Fo = self.fco(h1, S, aotwoint)
        Fa, Fb = self.fab(h1, S, aotwoint)
        gnref = 2.75678
        gnco = rohf.gradnorm(S, C, Dc, Do, Fc, Fo, 4, 1)
        gnab = rohf.gradnorm(S, C, Da, Db, Fa, Fb, 4, 1)
        assert gnco == approx(gnref)
        assert gnab == approx(gnref)

    def test_gradao_norm_co_ab(self, h1, S, aotwoint):

        Dc, Do = self.dco(h1, S)
        Da, Db = self.dab(h1, S)
        Fc, Fo = self.fco(h1, S, aotwoint)
        Fa, Fb = self.fab(h1, S, aotwoint)
        gnref = 2.75678
        gco = rohf.gradao(S, Dc, Do, Fc, Fo)
        gab = rohf.gradao(S, Da, Db, Fa, Fb)
        g2co = gco & (S.I@gco@S.I) * .5
        g2ab = gab & (S.I@gab@S.I) * .5
        assert math.sqrt(g2co) == approx(gnref)
        assert math.sqrt(g2ab) == approx(gnref)

    def test_feffs(self, h1, S, aotwoint):
        Dc, Do = self.dco(h1, S)
        Fc, Fo = self.fco(h1, S, aotwoint)
        Da, Db = self.dab(h1, S)
        Fa, Fb = self.fab(h1, S, aotwoint)
        Feff_co = rohf.jensen(S, Dc, Do, Fc, Fo)
        Feff_ab = S@rohf.Feff(Da@S, Db@S, S.I@Fa, S.I@Fb)
        npt.assert_allclose(
            np.array(Feff_co),
            np.array(Feff_ab),
        )
