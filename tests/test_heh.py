import numpy.testing as npt
from pytest import approx, fixture

from dalmisc.scf_iter import RoothanIterator


@fixture
def rohf():
    roo = RoothanIterator(
        electrons=3,
        tmpdir='tests/test_heh.d',
        threshold=1e-5,
        max_iterations=20,
        ms=1/2,
        )
    return roo


def test_initial_energy(rohf):
    iter(rohf)
    energy, gradient = next(rohf)
    assert energy == approx(-3.2691909239)
    energy, gradient = next(rohf)
    assert energy == approx(-3.2697129169)
    # energy, gradient = next(rohf)
    # assert energy == approx(-3.2697227469)


def test_initial_gradient(rohf):
    iter(rohf)
    energy, gradient = next(rohf)
    assert gradient == approx(4.413441177e-02)
    energy, gradient = next(rohf)
    assert gradient == approx(6.055260228e-03)
    # energy, gradient = next(rohf)
    # assert gradient == approx(0.000572018305)


def test_initial_mo(rohf):
    assert rohf.C is None
    iter(rohf)
    npt.assert_allclose(
        rohf.C[0],
        [[1.00029353, -0.0810076],
         [-0.00355856,  1.00356202]]
    )
    assert rohf.C[0] is rohf.C[1]
    next(rohf)
    npt.assert_allclose(
        rohf.C[1],
        [[0.99769772, -0.10839113],
         [0.02394631,  1.00328259]]
    )
    next(rohf)
    npt.assert_allclose(
        rohf.C[1],
        [[0.9980969, -0.10465161],
         [0.02018645,  1.00336528]],
        atol=1e-8,
    )


def test_initial_dens(rohf):
    assert rohf.Da is None
    irohf = iter(rohf)
    assert rohf.Da is None
    next(irohf)
    npt.assert_allclose(
        irohf.Da,
        [[1.00714939, -0.08485576],
         [-0.08485576,  1.00714939]]
    )
    npt.assert_allclose(
        irohf.Db,
        [[1.00058715e+00, -3.55960477e-03],
         [-3.55960477e-03,  1.26633508e-05]]
    )
    next(irohf)
    npt.assert_allclose(
        irohf.Da,
        [[1.00714939, -0.08485576],
         [-0.08485576,  1.00714939]]
    )
    npt.assert_allclose(
        irohf.Db,
        [[9.95400748e-01, 2.38911792e-02],
         [2.38911792e-02, 5.73425774e-04]]
    )


def test_initial_fock(rohf):
    assert rohf.Fa is None
    irohf = iter(rohf)
    assert rohf.Fa is None
    next(irohf)
    npt.assert_allclose(
        irohf.h1 + irohf.Fa,
        [[-0.87871515, -0.10484785],
         [-0.10484785, -0.46979131]]
    )
    npt.assert_allclose(
        irohf.h1 + irohf.Fb,
        [[-0.87621434, -0.09160072],
         [-0.09160072,  0.3045648]]
    )
    next(irohf)
    npt.assert_allclose(
        irohf.h1 + irohf.Fa,
        [[-0.88107255, -0.10487083],
         [-0.10487083, -0.46884594]]
    )
    npt.assert_allclose(
        irohf.h1 + irohf.Fb,
        [[-0.87606003, -0.09915466],
         [-0.09915466,  0.3031376]]

    )


def test_initial_gab(rohf):
    assert rohf.ga is None
    irohf = iter(rohf)
    assert rohf.ga is None
    next(irohf)
    npt.assert_allclose(
        irohf.ga,
        [[0., 0], [0., 0.]],
        atol=1e-8,
    )
    npt.assert_allclose(
        irohf.gb,
        [[0., -0.02198874], [0.02198874, 0.]],
        atol=1e-8,
    )
    next(irohf)
    npt.assert_allclose(
        irohf.ga,
        [[0., 0], [0., 0.]],
        atol=1e-8,
    )
    npt.assert_allclose(
        irohf.gb,
        [[0., 0.00301686495], [-0.00301686495, 0.]],
        atol=1e-8,
    )


def test_rohf(rohf):
    final_energy, final_gradient = list(rohf)[-1]
    assert final_energy == approx(-3.269722925573)
    assert final_gradient < 1e-5
