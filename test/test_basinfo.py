import unittest
import os
import numpy
from ..basinfo import BasInfo


class TestBasInfo(unittest.TestCase):
    
    def setUp(self):
        n, _ = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.bas_info = BasInfo(os.path.join(suppdir, 'SIRIUS.RST'))

    def test_nsym(self): 
        self.assertEqual(self.bas_info.nsym, 1)

    def test_nbas(self): 
        numpy.testing.assert_equal(self.bas_info.nbas, [5,0,0,0,0,0,0,0])

    def test_nbast(self): 
        self.assertEqual(self.bas_info.nbast, 5)

    def test_norb(self): 
        numpy.testing.assert_equal(self.bas_info.norb, [5,0,0,0,0,0,0,0])

    def test_norbt(self): 
        self.assertEqual(self.bas_info.norbt, 5)

    def test_ncmot(self):
        self.assertEqual(self.bas_info.ncmot, 25)

    def test_str(self):
        ref = """\
NSYM   :   1
NBAS   :   5  0  0  0  0  0  0  0
NORB   :   5  0  0  0  0  0  0  0
NRHF   :   1  0  0  0  0  0  0  0
IOPRHF :   0
"""
        self.assertEqual(str(self.bas_info), ref)

    def test_main(self):
        pass

if __name__ == "__main__":
    unittest.main()
