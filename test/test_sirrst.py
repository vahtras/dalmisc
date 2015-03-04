import unittest
import os
import numpy
from ..sirrst import SiriusRestart
from util.blocked import BlockDiagonalMatrix

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

    def test_dens(self):
        occupied = self.ref_cmo[:, :1]
        numpy.testing.assert_almost_equal(
            self.sirrst.get_rhf_density(), 
            2*occupied*occupied.T
        )

    def test_dens_symmetry(self):
        sir = SiriusRestart(os.path.join(self.suppdir, 'hf_S.SIRIUS.RST'))
        occupied = BlockDiagonalMatrix(sir.basinfo.nbas, sir.basinfo.nrhf)
        for nrhf, occ, cmo in zip(sir.basinfo.nrhf, occupied, sir.cmo):
            if nrhf > 0:
                occ[:, :] = cmo[:, :nrhf]
        full_occupied = occupied.unblock()
        full_density = 2*full_occupied*full_occupied.T
        numpy.testing.assert_almost_equal(
            sir.get_rhf_density(), full_density
        )


    def test_hf_S_symmetry(self):
        sirius_restart = SiriusRestart(os.path.join(self.suppdir, 'hf_S.SIRIUS.RST'))
        cmo = sirius_restart.cmo
        

