import unittest
import os
import numpy
from ..sirrst import SiriusRestart

class TestSirRst(unittest.TestCase):

    def setUp(self):
        self.suppdir = os.path.splitext(__file__)[0] + ".d"
        self.sirrst = SiriusRestart(os.path.join(self.suppdir, 'SIRIUS.RST'))
        self.ref_cmo = numpy.array([
           [ 0.71551428,  -0.72497592,    0.00000000,    0.00000000,    0.20543401],
           [ 0.00000000,   0.00000000,    0.00000000,    1.00000000,    0.00000000],
           [ 0.00000000,   0.00000000,    1.00000000,    0.00000000,    0.00000000],
           [-0.03156984,   0.09811871,    0.00000000,    0.00000000,    1.09913926],
           [ 0.55260971,   0.83525754,    0.00000000,    0.00000000,   -0.54355525]
            ])

    def test_cmo(self):
        numpy.testing.assert_almost_equal(self.sirrst.cmo[0], self.ref_cmo)

    def test_hf_S_symmetry(self):
        sirius_restart = SiriusRestart(os.path.join(self.suppdir, 'hf_S.SIRIUS.RST'))
        cmo = sirius_restart.cmo
        

