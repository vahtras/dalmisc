import  unittest
import os
import numpy
from util.full import matrix
from .. import two

class TwoTest(unittest.TestCase):

    def setUp(self):
        n, e = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.aotwoint = os.path.join(suppdir, "AOTWOINT")

    def test_basinfo(self):
        info = two.info(self.aotwoint)
        self.assertEqual(info["nsym"], 1)
        numpy.testing.assert_equal(info["nbas"], [7,0,0,0,0,0,0,0])
        self.assertEqual(info["lbuf"], 600)
        self.assertEqual(info["nibuf"], 1)
        self.assertEqual(info["nbits"], 8)


    def test_first_integral(self):
        for ig, g in two.list_integrals(self.aotwoint):
            break
        self.assertAlmostEqual(g, 4.78506540471)
        numpy.testing.assert_equal(ig, [1,1,1,1])


class TestBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        root, ext = n, e = os.path.splitext(__file__)
        cls.base_dir = root + ".d"
        
class TestH2O(TestBase):

    @classmethod
    def setUpClass(cls):
        super(TestH2O, cls).setUpClass()
        cls.subdir = "H2O"
        cls.dal = "hf"
        cls.mol = "H2O_ccpVDZ"
        cls.filename = "%s_%s.AOTWOINT" % (cls.dal, cls.mol)
        cls.tmpdir = os.path.join(cls.base_dir, cls.subdir)
        cls.aotwoint = os.path.join(cls.tmpdir, cls.filename)
        if not os.path.exists(cls.aotwoint):
            os.chdir(cls.tmpdir)
            os.system('dalton -get AOTWOINT hf H2O_ccpVDZ')

    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 11412)

    def setUp(self):
        self.d = matrix((24, 24))
        self.f = matrix((24, 24))

        self.d[:11, :5] = numpy.array([
             [ 2.00160706,   -0.01036708,    -0.00181429,     0.00677187,    -0.00241386],
             [-0.01036708,    1.69869523,    -0.09633899,    -0.27114604,    -0.10920586],
             [-0.00181429,   -0.09633899,     0.08874347,    -0.30570693,     0.01352778],
             [ 0.00677187,   -0.27114604,    -0.30570693,     1.28122865,    -0.01084007],
             [-0.00241386,   -0.10920586,     0.01352778,    -0.01084007,     0.00767067],
             [ 0.00006681,   -0.00823265,    -0.00653547,     0.02831209,    -0.00008716],
             [-0.00018734,   -0.00235146,     0.00228316,    -0.00791291,     0.00034070],
             [-0.00389444,    0.35145795,    -0.19737922,     0.62805954,    -0.03821221],
             [ 0.00283175,   -0.18331099,     0.08114304,    -0.24351160,     0.01800991],
             [ 0.00146798,   -0.04495921,     0.01951429,    -0.05823360,     0.00438197],
             [ 0.00123945,   -0.03923112,     0.00210759,     0.00671245,     0.00251033] 
            ])

        self.d[:11, 5:10] = numpy.array([
            [ 0.00006681,   -0.00018734,    -0.00389444,     0.00283175,     0.00146798],
            [-0.00823265,   -0.00235146,     0.35145795,    -0.18331099,    -0.04495921],
            [-0.00653547,    0.00228316,    -0.19737922,     0.08114304,     0.01951429],
            [ 0.02831209,   -0.00791291,     0.62805954,    -0.24351160,    -0.05823360],
            [-0.00008716,    0.00034070,    -0.03821221,     0.01800991,     0.00438197],
            [ 0.00062869,   -0.00016935,     0.01321749,    -0.00506051,    -0.00120868],
            [-0.00016935,    0.00005876,    -0.00506668,     0.00207976,     0.00050004],
            [ 0.01321749,   -0.00506668,     0.45083132,    -0.18868220,    -0.04545457],
            [-0.00506051,    0.00207976,    -0.18868220,     0.07988865,     0.01926631],
            [-0.00120868,    0.00050004,    -0.04545457,     0.01926631,     0.00464710],
            [ 0.00019991,    0.00005121,    -0.00787025,     0.00413572,     0.00101526] 
            ])

        self.d[:15, 10:15] = numpy.array([
            [ 0.00123945,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [-0.03923112,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00210759,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00671245,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00251033,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00019991,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00005121,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [-0.00787025,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00413572,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00101526,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00090670,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00000000,     1.69432975,     0.13331289,     0.03312803,     0.05610831],
            [ 0.00000000,     0.13331289,     0.01048930,     0.00260657,     0.00441470],
            [ 0.00000000,     0.03312803,     0.00260657,     0.00064773,     0.00109705],
            [ 0.00000000,     0.05610831,     0.00441470,     0.00109705,     0.00185805] 
            ])

        self.d[15:22, 15:20] = numpy.array([
            [ 1.01124471,    -0.14392643,     0.03731214,     0.78643508,    -0.26241828],
            [-0.14392643,     0.02048448,    -0.00531049,    -0.11193017,     0.03734895],
            [ 0.03731214,    -0.00531049,     0.00137671,     0.02901728,    -0.00968251],
            [ 0.78643508,    -0.11193017,     0.02901728,     0.61160284,    -0.20408012],
            [-0.26241828,     0.03734895,    -0.00968251,    -0.20408012,     0.06809762],
            [-0.03329508,     0.00473876,    -0.00122850,    -0.02589325,     0.00864008],
            [-0.04612493,     0.00656478,    -0.00170188,    -0.03587091,     0.01196943]
            ])

        self.d[15:22, 20:24] = numpy.array([
            [-0.03329508,    -0.04612493,     0.00000000,     0.00000000],
            [ 0.00473876,     0.00656478,     0.00000000,     0.00000000],
            [-0.00122850,    -0.00170188,     0.00000000,     0.00000000],
            [-0.02589325,    -0.03587091,     0.00000000,     0.00000000],
            [ 0.00864008,     0.01196943,     0.00000000,     0.00000000],
            [ 0.00109624,     0.00151866,     0.00000000,     0.00000000],
            [ 0.00151866,     0.00210385,     0.00000000,     0.00000000]
            ])

        self.f[:11, :5] = numpy.array([
            [12.45937161,    -0.63704828,     1.95558879,     0.00360350,    -0.00052098],
            [-0.63704828,     6.82369804,     5.69718410,     0.10818033,     0.06591144],
            [ 1.95558879,     5.69718410,     5.79457260,     0.08665063,     0.06521453],
            [ 0.00360350,     0.10818033,     0.08665063,     6.90560486,     4.26716496],
            [-0.00052098,     0.06591144,     0.06521453,     4.26716496,     4.62120096],
            [ 0.00038575,     0.00182652,     0.00087351,     0.05815816,     0.04421790],
            [-0.00469206,     0.00404940,     0.00442371,     0.03779009,     0.03337247],
            [ 1.26213056,     5.56856988,     6.25911109,     2.47075981,     3.05828625],
            [ 1.34757129,     5.72503437,     6.48855702,     1.44115327,     2.00553367],
            [-1.67621948,    -4.67125327,    -4.25488310,    -4.05369485,    -3.01795390],
            [-1.29869229,    -3.54818988,    -3.22870004,     0.15474749,     1.36984177] 
            ])

        self.f[:11, 5:10] = numpy.array([
            [ 0.00038575,    -0.00469206,     1.26213056,     1.34757129,    -1.67621948],
            [ 0.00182652,     0.00404940,     5.56856988,     5.72503437,    -4.67125327],
            [ 0.00087351,     0.00442371,     6.25911109,     6.48855702,    -4.25488310],
            [ 0.05815816,     0.03779009,     2.47075981,     1.44115327,    -4.05369485],
            [ 0.04421790,     0.03337247,     3.05828625,     2.00553367,    -3.01795390],
            [ 7.24153911,     0.01827008,     0.12245273,     0.02494436,    -2.15163129],
            [ 0.01827008,     7.23051852,    -0.91671029,    -0.12495988,    -0.00063511],
            [ 0.12245273,    -0.91671029,    11.80281335,    11.26687675,    -3.48548743],
            [ 0.02494436,    -0.12495988,    11.26687675,    12.29454973,    -3.31072622],
            [-2.15163129,    -0.00063511,    -3.48548743,    -3.31072622,    12.99380806],
            [ 2.43758337,     2.36581597,    -1.39184409,    -1.31458787,     0.73448332] 
            ])

        self.f[:15, 10:15] = numpy.array([
            [-1.29869229,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [-3.54818988,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [-3.22870004,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.15474749,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 1.36984177,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 2.43758337,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 2.36581597,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [-1.39184409,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [-1.31458787,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.73448332,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [10.28246104,     0.00000000,     0.00000000,     0.00000000,     0.00000000],
            [ 0.00000000,     6.83817059,     4.24466177,     0.05123594,     3.21192054],
            [ 0.00000000,     4.24466177,     4.61821686,     0.04020765,     3.63810087],
            [ 0.00000000,     0.05123594,     0.04020765,     7.22649702,     2.36576057],
            [ 0.00000000,     3.21192054,     3.63810087,     2.36576057,     9.74473836] 
            ])

        self.f[15:22, 15:20] = numpy.array([
            [ 6.96210169,     4.29189102,    -0.02339239,     3.11382994,     1.78195656],
            [ 4.29189102,     4.62997234,    -0.02821633,     3.86436194,     2.49255893],
            [-0.02339239,    -0.02821633,     7.23092994,     1.42794681,     0.20727286],
            [ 3.11382994,     3.86436194,     1.42794681,     5.18731515,     3.04609968],
            [ 1.78195656,     2.49255893,     0.20727286,     3.04609968,     2.51926982],
            [-1.87285341,    -0.16875709,    -2.40405120,     0.62958017,     0.96777694],
            [-3.98148145,    -2.94544371,    -0.61408162,    -0.68318845,    -0.38697081]
            ])

        self.f[15:24, 20:24] = numpy.array([
            [-1.87285341,   -3.98148145,     0.00000000,     0.00000000],
            [-0.16875709,   -2.94544371,     0.00000000,     0.00000000],
            [-2.40405120,   -0.61408162,     0.00000000,     0.00000000],
            [ 0.62958017,   -0.68318845,     0.00000000,     0.00000000],
            [ 0.96777694,   -0.38697081,     0.00000000,     0.00000000],
            [ 7.07648330,    0.26707594,     0.00000000,     0.00000000],
            [ 0.26707594,    9.16223276,     0.00000000,     0.00000000],
            [ 0.00000000,    0.00000000,     7.23712021,     3.03286880],
            [ 0.00000000,    0.00000000,     3.03286880,     8.88253901]
            ])

    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(two.fock(self.d, filename=self.aotwoint, f2py=False), self.f)

    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(two.fock(self.d, filename=self.aotwoint, f2py=True), self.f)

if __name__ == "__main__":
    unittest.main()

    
