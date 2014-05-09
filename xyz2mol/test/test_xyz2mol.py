import unittest
from ..xyz2mol import xyz2mol, extract_element

xyz = """9
Ethanol
C1    9.136479   11.724970   51.500927
H1    9.294462   12.628635   50.910741
H2    8.160480   11.311733   51.245960
H3    9.126987   11.998458   52.556137
C2   10.230421   10.715737   51.222941
H4   11.209335   11.140253   51.472204
H5   10.240910   10.452121   50.159449
O     9.976946    9.554726   52.018512
H6   10.663980    8.903366   51.853128
"""

molref = """BASIS
cc-pVDZ
Ethanol
=======
Atomtypes=3
Charge=6 Atoms=2
C1    9.136479   11.724970   51.500927
C2   10.230421   10.715737   51.222941
Charge=1 Atoms=6
H1    9.294462   12.628635   50.910741
H2    8.160480   11.311733   51.245960
H3    9.126987   11.998458   52.556137
H4   11.209335   11.140253   51.472204
H5   10.240910   10.452121   50.159449
H6   10.663980    8.903366   51.853128
Charge=8 Atoms=1
O     9.976946    9.554726   52.018512
"""

class Test2Mol(unittest.TestCase):

    def test_extract_numbered_element(self):
        atom_line = "He1 0 0 0"
        element = extract_element(atom_line)
        self.assertEqual(element, "He")

    def test_extract_unnumbered_element(self):
        atom_line = "O 0 0 0"
        element = extract_element(atom_line)
        self.assertEqual(element, "O")

    def test_raise_unknown_element(self):
        self.assertRaises(Exception, extract_element, "Yo 1 2 3")

    def test_xyz2mol(self):
        mol = xyz2mol(xyz)
        print mol, molref
        assert mol == molref


if __name__ == "__main__":
    test_xyz2mol()
