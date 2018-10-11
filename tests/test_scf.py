import pytest
from pytest import approx

from dalmisc.scf_iter import *

@pytest.fixture
def scf():
    return SCFIterator()

def test_fix(scf):
    assert isinstance(scf, SCFIterator)

@pytest.fixture
def roothan():
    return RoothanIterator(
        electrons=10,
        max_iterations=15,
        threshold=1e-5,
        tmpdir=os.path.splitext(__file__)[0] + ".d",
    )

def test_fix2(roothan):
    assert isinstance(roothan, RoothanIterator)

def test_defaults():
    roothan = RoothanIterator()
    assert roothan.max_iterations == 10
    assert roothan.threshold == 1e-3

def test_setup(roothan):
    assert roothan.max_iterations == 15
    assert roothan.threshold == 1e-5

def test_iter(roothan):
    assert iter(roothan) is roothan

def test_initial(roothan):
    assert roothan.C is None

def test_stop_threshold(roothan):
    with pytest.raises(StopIteration):
        roothan.gn = lambda: 0
        next(roothan)

def test_stop_iterations(roothan):
    with pytest.raises(StopIteration):
        roothan.gn = lambda: 1
        roothan.it = roothan.max_iterations
        next(roothan)

def test_initial_guess(roothan):
    initial = iter(roothan)
    assert initial.energy() == approx(-73.2292918615)

def test_one_fockit(roothan):
    initial = iter(roothan)
    assert next(initial)[0] == approx(-74.946960167351)

def test_initial_electrons(roothan):
    initial = iter(roothan)
    assert initial.na == 5
    assert initial.nb == 5

def test_Z(roothan):
    assert roothan.Z is None
    assert iter(roothan).Z == approx(9.263515863231)

def test_h1(roothan):
    assert roothan.h1 is None
    assert iter(roothan).h1 is not None

def test_overlap(roothan):
    assert roothan.S is None
    assert iter(roothan).S is not None

def test_densities(roothan):
    initial = iter(roothan)
    Da, Db = initial.densities()
    assert Da&initial.S == approx(5)
    assert Db&initial.S == approx(5)

def test_h1(roothan):
    initial = iter(roothan)
    Da, Db = initial.densities()
    assert initial.h1&(Da + Db) == approx(-127.45439681043854)
    




