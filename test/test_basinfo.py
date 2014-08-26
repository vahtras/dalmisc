import unittest
import os
from dalmisc.basinfo import BasInfo

class TestBasInfo(unittest.TestCase):
    
    def setUp(self):
        n, _ = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.bas_info = BasInfo(os.path.join(suppdir, 'SIRIUS.RST'))

    def test_nsym(self): 
        self.assertEqual(self.bas_info.nsym, 1)

    def test_nbas(self): 
        self.assertEqual(self.bas_info.nbas[0], 5)

    def test_nbast(self): 
        self.assertEqual(self.bas_info.nbast, 5)

    def test_ncmot(self):
        self.assertEqual(self.bas_info.ncmot, 25)

if __name__ == "__main__":
    unittest.main()
