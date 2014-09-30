import unittest
import numpy as np
from ..scan_dalton import *

class TestScan(unittest.TestCase):

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
 Sym       Hartree-Fock orbital energies

1 A     -20.55302081    -1.32934124    -0.69144670    -0.56387102    -0.49212361
          0.18304191     0.25403508     0.77887068     0.84069361     1.16451072

...

@ Excitation energy :  0.20348808     au
@ Excitation energy :  0.27095312     au
@ Excitation energy :  0.27500777     au
@ Excitation energy :  0.31446588     au
@ Excitation energy :  0.32292011     au

...
@    Final HF energy:             -76.025681483940
...
     XDIPLEN  total        :-2.03844916D-15
     YDIPLEN  total        : 1.53769063D-15
     ZDIPLEN  total        :    -0.81457755
     XXSECMOM total        :     7.26004973
     XYSECMOM total        : 2.11454908D-16
     XZSECMOM total        : 2.11454908D-16
     YYSECMOM total        :     5.25404533
     YZSECMOM total        :-1.84889923D-18
     ZZSECMOM total        :     6.40304140
...
     KINENERG inactive part:     0.00000000
     KINENERG active part  :    10.19430190
     KINENERG total        :    10.19430190
...
     D1-SO XX inactive part:     0.00000000
     D1-SO XX active part  :    11.25069842
     D1-SO XX total        :    11.25069842
...
@G GC2    -0.000671 -0.000807 -0.000807 -0.000056 -0.000055  0.000056  0.000056
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

@ Transition operator type:    XDIPLEN
@ STATE NO:    1 *TRANSITION MOMENT: -7.59491062E-16 *ENERGY(eV):   9.0419086
@ STATE NO:    2 *TRANSITION MOMENT:  4.32992333E-16 *ENERGY(eV):   9.0543615
@ STATE NO:    3 *TRANSITION MOMENT:  2.62890021E-17 *ENERGY(eV):   15.594527
@ STATE NO:    4 *TRANSITION MOMENT:  1.61762067E-16 *ENERGY(eV):   15.763878

@ Transition operator type:    YDIPLEN
@ STATE NO:    1 *TRANSITION MOMENT:  0.23546966     *ENERGY(eV):   9.0419086
@ STATE NO:    2 *TRANSITION MOMENT:  0.43727594     *ENERGY(eV):   9.0543615
@ STATE NO:    3 *TRANSITION MOMENT:  2.99603676E-03 *ENERGY(eV):   15.594527
@ STATE NO:    4 *TRANSITION MOMENT:  2.47044333E-03 *ENERGY(eV):   15.763878

@ Transition operator type:    ZDIPLEN
@ STATE NO:    1 *TRANSITION MOMENT:  1.16578641E-15 *ENERGY(eV):   9.0419086
@ STATE NO:    2 *TRANSITION MOMENT:  5.73781327E-16 *ENERGY(eV):   9.0543615
@ STATE NO:    3 *TRANSITION MOMENT:  1.41349350E-17 *ENERGY(eV):   15.594527
@ STATE NO:    4 *TRANSITION MOMENT:  6.63118919E-15 *ENERGY(eV):   15.763878

...
@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< XANGMOM  ; XANGMOM  >> =  7.647374881166D+01
@ -<< XANGMOM  ; X1SPNORB >> =  2.065248721972D-01
@ -<< X1SPNORB ; X1SPNORB >> =  6.963789419286D-04

@ -<< XANGMOM  ; X2SPNORB >> = -8.426341638486D-02
...
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

    def test_nuclear_quadrupole_moment(self):
        PN = get_nuclear_quadrupole_moment(self.filename)
        np.testing.assert_almost_equal(PN, (
            [4.21864105, 0.00000,  0.00000,  0.00000,  0.00000,  1.66922574],
            ))

    def test_electronic_diplen(self):
        np.testing.assert_almost_equal(
            get_electronic_dipole_moment(self.filename), 
            ([ 2.03844916e-15,   -1.53769063e-15,  8.14577550e-01],)
            )

    def test_electronic_no_diplen_raises_exception(self):
        self.assertRaises(
            NotFoundError, get_electronic_dipole_moment, '/dev/null'
            )

    def test_electronic_quadrupole(self):
        np.testing.assert_almost_equal(
            get_electronic_quadrupole_moment(self.filename), 
            ([ -7.26004973, 0, 0, -5.25404533, 0, -6.40304140],)
            )

    def test_electronic_no_quadrupole_raises_exception(self):
        self.assertRaises(
            NotFoundError, get_electronic_quadrupole_moment, '/dev/null'
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

    def test_nopol(self):
        self.assertRaises(NotFoundError, get_polarizability, "/dev/null")

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

    def notest_nuclear_quadrupole(self):
        ref = (4.21864,  0.00000,  0.00000,  0.00000,  0.00000,  2.02330)
        q = get_nuclear_quadrupole_moment(self.filename)
        np.testing.assert_allclose(q, ref)

    def test_get_final_hf_energy(self):
        ref = (-76.02568148394,)
        e, = get_final_energy(self.filename)
        self.assertAlmostEqual(e, ref)

    def test_not_found(self):
        def func():
            e = get_g_rmc('/dev/null')
        self.assertRaises(NotFoundError, func)


    def test_get_g_rmc(self):
        ref = -0.001087
        e, = get_g_rmc(self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_get_g_gc1(self):
        ref = ((0.000599, 0, 0), (0, 0, 0), (0, 0, 0))
        e, = get_g_gc1(self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_get_g_gc2(self):
        ref = (
            (-0.000671, -0.000056, 0.000056), 
            (-0.000055, -0.000807, 0), 
            (0.0000560, 0, -0.000807)
        )

        e, = get_g_gc2(self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_get_g_oz1(self):
        ref = ((0.41304974, 0, 0), (0, 0, 0), (0, 0, 0))
        e, = get_g_oz1(self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_get_g_oz2(self):
        ref = ((-0.16852683, 0, 0), (0, 0, 0), (0, 0, 0))
        e, = get_g_oz2(self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_get_orbital_energies(self):
        ref = (
            -20.55302081, -1.32934124, -0.69144670, -0.56387102, -0.49212361,
              0.18304191,  0.25403508,  0.77887068,  0.84069361,  1.16451072
            )
        e, = get_orbital_energies(self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_get_excitation_energies(self):
        ref = (0.20348808, 0.27095312, 0.27500777, 0.31446588, 0.32292011)
        e, = get_excitation_energies(self.filename)
        np.testing.assert_almost_equal(e, ref)


    def notest_get_transition_moments(self):
        ref = (0.23546966, 0.43727594, 2.99603676E-03, 2.47044333E-03)
        e, = get_transition_moments("YDIPLEN", self.filename)
        np.testing.assert_almost_equal(e, ref)

    def test_transition_operator_pattern(self):
        self.assertIsNotNone(
            re.match(
                transition_operator_pattern("XDIPLEN"),
                "@ Transition operator type:    XDIPLEN"
                )
            )
            


if __name__ == "__main__":
    unittest.main()
