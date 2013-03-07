
import os
import numpy as np
from daltools import one, dens
from dalmisc import rohf, fockab


def setup():
    global h1, S, EN, aooneint, aotwoint, na, nb
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_rohf.d')
    aooneint = os.path.join(suppdir, "AOONEINT")
    aotwoint = os.path.join(suppdir, "AOTWOINT")
    h1 = one.read("ONEHAMIL", aooneint).unpack().unblock()
    S = one.read("OVERLAP", aooneint).unpack().unblock()
    EN = one.readhead(aooneint)["potnuc"]
    na = 5
    nb = 5

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)

def test_potnuc():
    assert np.allclose(EN, 9.055004525638)

def test_h1diag_initial_mo(): 
    import ref_rohf as ref
    Cmo = dens.cmo(h1, S)
    Ca = Cmo[:, :na]
    Cb = Cmo[:, :nb]
    assert_(Ca, ref.Ca)

def test_h1diag_initial_dens(): 
    import ref_rohf as ref
    Cmo = dens.cmo(h1, S)
    Ca = Cmo[:, :na]
    Cb = Cmo[:, :nb]
    Da  = Ca*Ca.T
    Db  = Cb*Cb.T
    assert_(Da, ref.Da)

def test_h1diag_initial_fock(): 
    import ref_rohf as ref
    Cmo = dens.cmo(h1, S)
    Ca = Cmo[:, :na]
    Cb = Cmo[:, :nb]
    Da  = Ca*Ca.T
    Db  = Cb*Cb.T
    assert_(Da, ref.Da)


def test_h1diag_initial_energy(): 
    Cmo = dens.cmo(h1, S)
    Ca = Cmo[:, :na]
    Cb = Cmo[:, :nb]
    Da  = Ca*Ca.T
    Db  = Cb*Cb.T
    Fa, Fb = fockab.fockab((Da, Db), filename=aotwoint)
    E = EN + rohf.energy(Da, Db, h1, Fa, Fb)
    Eref = -73.240064311328
    assert_(E, Eref)

if __name__ == "__main__":
    setup()
