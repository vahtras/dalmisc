from pytest import approx, mark
from dalmisc.scf_iter import RoothanIterator


def test_rhf():
    roo = RoothanIterator(
        electrons=10,
        tmpdir='tests/test_h2o.d',
        threshold=1e-5,
        max_iterations=20,
        )
    final_energy, final_gradient = list(roo)[-1]
    assert final_energy == approx(-74.9615984420)
    assert final_gradient < 1e-5


@mark.skip()
def test_rohf():
    roo = RoothanIterator(
        electrons=9,
        tmpdir='tests/test_h2o.d',
        threshold=1e-5,
        max_iterations=20,
        ms=1/2,
        )
    final_energy, final_gradient = list(roo)[-1]
    assert final_energy == approx(-74.651129646549)
    assert final_gradient < 1e-5
