from pytest import approx
from dalmisc.scf_iter import RoothanIterator


def test_rhf():
    roo = RoothanIterator(
        electrons=10,
        tmpdir='tests/test_h2o.d',
        threshold=1e-5,
        max_iterations=20,
        )
    for i, (e, gn) in enumerate(roo, start=1):
        print(f'{i:2d}: {e:14.10f} {gn:.3e}')
    assert e == approx(-74.9615984420)


def test_rohf():
    roo = RoothanIterator(
        electrons=9,
        tmpdir='tests/test_h2o.d',
        threshold=1e-5,
        max_iterations=20,
        ms=1/2,
        )
    for i, (e, gn) in enumerate(roo, start=1):
        print(f'{i:2d}: {e:14.10f} {gn:.3e}')
    assert e == approx(-74.651129646549)
