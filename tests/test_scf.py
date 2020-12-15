import os
import pathlib

import pytest

from dalmisc.scf_iter import (
    SCFIterator,
    RoothanIterator,
    URoothanIterator,
    DaltonOutputIterator,
)


@pytest.fixture
def scf():
    return SCFIterator()


def test_fixture_scf(scf):
    assert isinstance(scf, SCFIterator)


@pytest.fixture
def roothan():
    _roothan = RoothanIterator(
        electrons=10,
        max_iterations=15,
        threshold=1e-5,
        tmpdir=os.path.splitext(__file__)[0] + ".d",
    )
    return _roothan


@pytest.fixture
def rohf_roothan():
    return RoothanIterator(
        electrons=9,
        max_iterations=15,
        threshold=1e-5,
        tmpdir=os.path.splitext(__file__)[0] + ".d",
        ms="1/2",
    )


@pytest.fixture
def uhf_roothan():
    return URoothanIterator(
        electrons=9,
        max_iterations=15,
        threshold=1e-5,
        tmpdir=os.path.splitext(__file__)[0] + ".d",
        ms="1/2",
    )


@pytest.fixture
def dalton_output():
    test_dir = pathlib.Path(__file__).with_suffix('.d')
    return DaltonOutputIterator(
        dalton_output=open(test_dir/"hf_h2o.out")
    )


def test_fixture_roothan(roothan):
    assert isinstance(roothan, RoothanIterator)
    assert isinstance(roothan, SCFIterator)


def test_rohf_fix(rohf_roothan):
    assert isinstance(rohf_roothan, RoothanIterator)
    assert isinstance(rohf_roothan, SCFIterator)


def test_defaults():
    roothan = RoothanIterator()
    assert roothan.max_iterations == 20
    assert roothan.threshold == 1e-5


def test_setup(roothan):
    assert roothan.max_iterations == 15
    assert roothan.threshold == 1e-5


def test_rohf_setup(rohf_roothan):
    assert rohf_roothan.max_iterations == 15
    assert rohf_roothan.threshold == 1e-5


def test_converged(roothan):
    roothan.it = 0
    assert not roothan.converged()
    roothan.it = 1
    roothan.gn = lambda: 0
    assert roothan.converged()


def test_rohf_converged(rohf_roothan):
    rohf_roothan.it = 0
    assert not rohf_roothan.converged()
    rohf_roothan.it = 1
    rohf_roothan.gn = lambda: 0
    assert rohf_roothan.converged()


def test_iter(roothan):
    assert iter(roothan) is roothan


def test_rohf_iter(rohf_roothan):
    assert iter(rohf_roothan) is rohf_roothan


def test_initial(roothan):
    assert roothan.C is None


def test_rohf_initial(rohf_roothan):
    assert rohf_roothan.C is None


def test_stop_threshold(roothan):
    with pytest.raises(StopIteration):
        roothan.it = 1
        roothan.gn = lambda: 0
        next(roothan)


def test_rohf_stop_iterations(roothan):
    with pytest.raises(StopIteration):
        roothan.gn = lambda: 1
        roothan.it = roothan.max_iterations
        next(roothan)


def test_initial_guess(roothan):
    initial_energy, _ = next(iter(roothan))
    assert initial_energy == pytest.approx(-73.2292918615)


def test_rohf_initial_guess(rohf_roothan):
    initial_energy, _ = next(iter(rohf_roothan))
    assert initial_energy == pytest.approx(-73.4538472272)


def test_one_fockit(roothan):
    scf = iter(roothan)
    next(scf)
    next(scf)
    assert scf.energy() == pytest.approx(-74.946960167351)


def test_rohf_one_fockit(rohf_roothan):
    scf = iter(rohf_roothan)
    next(scf)
    next(scf)
    assert scf.energy() == pytest.approx(-74.64791006861331)


def test_initial_electrons(roothan):
    initial = iter(roothan)
    assert initial.na == 5
    assert initial.nb == 5


def test_initial_rohf_electrons(rohf_roothan):
    initial = iter(rohf_roothan)
    assert initial.na == 5
    assert initial.nb == 4


def test_Z(roothan):
    assert roothan.Z is None
    assert iter(roothan).Z == pytest.approx(9.263515863231)


def test_overlap(roothan):
    assert roothan.S is None
    assert iter(roothan).S is not None


def test_densities(roothan):
    initial = iter(roothan)
    next(initial)
    assert initial.Da & initial.S == pytest.approx(5)
    assert initial.Db & initial.S == pytest.approx(5)


def test_rohf_densities(rohf_roothan):
    initial = iter(rohf_roothan)
    next(initial)
    assert initial.Da & initial.S == pytest.approx(5)
    assert initial.Db & initial.S == pytest.approx(4)


def test_h1(roothan):
    scf = iter(roothan)
    next(scf)
    assert scf.h1 & (scf.Da + scf.Db) == pytest.approx(-127.45439681043854)


def test_rhf_final(roothan):
    final_energy, _ = list(iter(roothan))[-1]
    assert final_energy == pytest.approx(-74.961598442034)


def test_rohf_final(rohf_roothan):
    final_energy, _ = list(iter(rohf_roothan))[-1]
    assert final_energy == pytest.approx(-74.651129646549)


def test_uhf_final(uhf_roothan):
    final_energy, _ = list(iter(uhf_roothan))[-1]
    assert final_energy == pytest.approx(-74.6531282076)


def test_daltin_initial(dalton_output):
    initial_energy, initial_gradient = next(iter(dalton_output))
    assert initial_energy == pytest.approx(-73.2292918615)
    assert initial_gradient == pytest.approx(2.98017)


def test_daltin_final(dalton_output):
    final_energy, final_gradient = list(dalton_output)[-1]
    assert final_energy == pytest.approx(-74.9615984420)
    assert final_gradient == pytest.approx(6.72407e-06)
