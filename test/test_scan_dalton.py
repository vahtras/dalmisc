from unittest import TestCase
import numpy as np
from ..scan_dalton import *

class TestScan(TestCase):

    def setUp(self):
        self.filename = '/tmp/testscan.out'
        out_string = \
"""
...
  Atomic type no.    1
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Basis set file used for this atomic type with Z =   8 :
     "/home/olav/dev/dalton/git/build/master/basis/cc-pVDZ"

  Atomic type no.    2
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Basis set file used for this atomic type with Z =   1 :
...
  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:    9
  O       :     1  x   0.0000000000    2  y   0.0000000000    3  z  -0.2249058930
  H1      :     4  x   1.4523500000    5  y   0.0000000000    6  z   0.8996230000
  H2      :     7  x  -1.4523500000    8  y   0.0000000000    9  z   0.8996230000
...
     XDIPLEN  total        :-2.03844916D-15
     YDIPLEN  total        : 1.53769063D-15
     ZDIPLEN  total        :    -0.81457755
...
 Center-of-mass coordinates (a.u.):    0.000000    0.000000   -0.099054
...
@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.00000):     7.21103625278    
@ QRLRVE:  << YDIPLEN  ; XDIPLEN  >> (   0.00000):    3.052337643628E-15
@ QRLRVE:  << ZDIPLEN  ; XDIPLEN  >> (   0.00000):   -2.229805424566E-15
@ QRLRVE:  << XDIPLEN  ; YDIPLEN  >> (   0.00000):    3.072792478754E-15
@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.00000):     3.03446384360    
@ QRLRVE:  << ZDIPLEN  ; YDIPLEN  >> (   0.00000):    4.752225289844E-15
@ QRLRVE:  << XDIPLEN  ; ZDIPLEN  >> (   0.00000):   -6.118700306251E-15
@ QRLRVE:  << YDIPLEN  ; ZDIPLEN  >> (   0.00000):    4.387168195962E-15
@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.00000):     5.22710462524    
... 
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;X,X) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;X,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;X,X) =    -18.48299798
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,X) = beta(Y,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Y,X) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Y,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,X) = beta(Z,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,X) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,X) =    -18.48299918
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;X,Y) = beta(Y,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;X,Y) = beta(Y,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;X,Y) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,Y) = beta(Y,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,Y) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Y,Y) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Y,Y) =     -2.33649400
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Y) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Y) = beta(Z,Y,Y)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Y) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Y) =     -2.33649395
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,Y) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;X,Z) = beta(Z,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;X,Z) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;X,Z) = beta(Z,Z,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,Z) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Y,Z) = beta(Z,Y,Y)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Y,Z) = beta(Z,Z,Y)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Z) = beta(Z,Z,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Z) = beta(Z,Z,Y)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Z) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Z) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,Z) =    -11.17349291
"""

        with open(self.filename, 'w') as out_file:
            out_file.write(out_string)

    def test_coordinates(self):
        np.testing.assert_almost_equal(
            get_coordinates(self.filename), (
            [[0.0000000000, 0.0000000000, -0.2249058930],
             [1.4523500000, 0.0000000000, 0.8996230000],
             [-1.4523500000, 0.0000000000, 0.8996230000]
            ],
            )
        )

    def test_coordinate_sets(self):
        coor = get_coordinates(self.filename, self.filename)
        np.testing.assert_almost_equal(
            coor[0],
            [[0.0000000000, 0.0000000000, -0.2249058930],
             [1.4523500000, 0.0000000000, 0.8996230000],
             [-1.4523500000, 0.0000000000, 0.8996230000]]
        )
        np.testing.assert_almost_equal(
            coor[1],
            [[0.0000000000, 0.0000000000, -0.2249058930],
             [1.4523500000, 0.0000000000, 0.8996230000],
             [-1.4523500000, 0.0000000000, 0.8996230000]]
        )


    def test_nuclear_charges(self):
        np.testing.assert_almost_equal(
            get_nuclear_charges(self.filename),
            ([8, 1, 1],)
        )

    def test_total_dipole_moment(self):
        np.testing.assert_almost_equal(
            get_total_dipole_moment(self.filename), 
            ([ 2.03844916e-15,   -1.53769063e-15,  0.814576406],)
            )

    def test_nuclear_dipole_moment(self):
        PN = get_nuclear_dipole_moment(self.filename)
        np.testing.assert_almost_equal(PN, ([0, 0, -1.144e-6],))

    def test_electronic_diplen(self):
        np.testing.assert_almost_equal(
            get_electronic_dipole_moment(self.filename), 
            ([ 2.03844916e-15,   -1.53769063e-15,  8.14577550e-01],)
            )

    def test_pol(self):
        np.testing.assert_almost_equal(
            get_polarizability(self.filename),
            ([
                [  7.21103625e+00,  3.07279248e-15, -6.11870031e-15],
                [  3.07279248e-15,  3.03446384e+00,  4.38716820e-15],
                [ -6.11870031e-15,  4.38716820e-15,  5.22710463e+00]
            ],)
            )

    def test_hyp(self):
        np.testing.assert_almost_equal(
            get_hyperpolarizability(self.filename), (
            [[[  0.,         -0.,        -18.48299918],
              [ -0.,          0.,         -0.        ],
              [-18.48299918, -0.,         -0.        ]],

             [[ -0.,          0.,         -0.        ],
              [  0.,         -0.,         -2.33649395],
              [ -0.,         -2.33649395, -0.        ]],

             [[-18.48299918, -0.,         -0.        ],
              [ -0.,         -2.33649395, -0.        ],
              [ -0.,         -0.,        -11.17349291]]]
            ,)
            )


    def test_xyz_to_tuple(self):
        self.assertTupleEqual(xyz_to_tuple('xyz'), (0, 1, 2))

    def test_pot_output(self):
        pot = generate_pot_output(self.filename, fmt="%8.3f").strip()
        ref = "1   0.000   0.000  -0.099   0.000   0.000  -0.000   0.815   7.211   0.000  -0.000   3.034   0.000   5.227   0.000  -0.000 -18.483   0.000  -0.000  -0.000  -0.000  -2.336  -0.000 -11.173"
        self.assertEqual(pot, ref)

    def test_get_cm(self):
        ref = ((.000000,   0.000000,  -0.099054),)
        cm = get_center_of_mass(self.filename)
        np.testing.assert_allclose(cm ,ref)

