import os
import numpy as np
import pytest

from dalmisc.expect import value

@pytest.fixture
def suppdir():
    n, e = os.path.splitext(__file__)
    return n + ".d"

def test_z(suppdir):
    zdiplen, = value('ZDIPLEN', tmpdir=suppdir)
    assert np.allclose(zdiplen, -1.91204774)
