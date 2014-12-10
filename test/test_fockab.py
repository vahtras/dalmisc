import os
import numpy as np
from daltools.util.full import matrix, init
from daltools import prop
from ..two import fockab, fock

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"

faref = init([
   [2.02818057, 0.26542036,  0.00000000,  0.06037429, 0.00000000, 0.00000000],
   [0.26542036, 0.74551226,  0.00000000,  0.28646605, 0.00000000, 0.00000000],
   [0.00000000, 0.00000000,  1.01061061, -0.47343699, 0.00000000, 0.00000000],
   [0.06037429, 0.28646605, -0.47343699,  0.86039561, 0.00000000, 0.00000000],
   [0.00000000, 0.00000000,  0.00000000,  0.00000000, 1.01061061, 0.00000000],
   [0.00000000, 0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.01061061],
   ])

fbref = init([
   [2.07812203, 0.35828051,  0.00000000,  0.09571548, 0.00000000, 0.00000000],
   [0.35828051, 1.03607456,  0.00000000,  0.40151695, 0.00000000, 0.00000000],
   [0.00000000, 0.00000000,  1.07479494, -0.50700370, 0.00000000, 0.00000000],
   [0.09571548, 0.40151695, -0.50700370,  0.93374109, 0.00000000, 0.00000000],
   [0.00000000, 0.00000000,  0.00000000,  0.00000000, 1.07479494, 0.00000000],
   [0.00000000, 0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.07479494],
   ])

def test_fab_p(): 
    da = matrix((6, 6))
    db = matrix((6, 6))
    da[0, 0] = 1 
    da[1, 1] = 1 
    db[0, 0] = 1 
    fa, fb = fockab((da, db), filename = os.path.join(suppdir, "AOTWOINT"), f2py=False)



    assert_(fa, faref)
    assert_(fb, fbref)

def test_fab_f(): 
    da = matrix((6, 6))
    db = matrix((6, 6))
    da[0, 0] = 1 
    da[1, 1] = 1 
    db[0, 0] = 1 
    fa, fb = fockab((da, db), filename = os.path.join(suppdir, "AOTWOINT"), f2py=True)


    assert_(fa, faref)
    assert_(fb, fbref)

def test_f_p(): 
    da = matrix((6, 6))
    db = matrix((6, 6))
    da[0, 0] = 1 
    da[1, 1] = 1 
    db[0, 0] = 1 
    d = da + db
    f = fock(d, filename = os.path.join(suppdir, "AOTWOINT"), f2py=False)

    fref = 0.5*(faref+fbref)
    assert_(f, fref)

def test_f_f(): 
    da = matrix((6, 6))
    db = matrix((6, 6))
    da[0, 0] = 1 
    da[1, 1] = 1 
    db[0, 0] = 1 
    d = da + db
    f = fock(d, filename = os.path.join(suppdir, "AOTWOINT"), f2py=True)

    fref = 0.5*(faref+fbref)
    assert_(f, fref)
    

if __name__ == "__main__":
    setup()
    test_fab_f()

#   faref = init([
#      [2.02818057, 0.26542036,  0.00000000,  0.06037429, 0.00000000, 0.00000000],
#      [0.26542036, 0.74551226,  0.00000000,  0.28646605, 0.00000000, 0.00000000],
#      [0.00000000, 0.00000000,  1.01061061, -0.47343699, 0.00000000, 0.00000000],
#      [0.06037429, 0.28646605, -0.47343699,  0.86039561, 0.00000000, 0.00000000],
#      [0.00000000, 0.00000000,  0.00000000,  0.00000000, 1.01061061, 0.00000000],
#      [0.00000000, 0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.01061061],
#      ])

#   fbref = init([
#      [2.07812203, 0.35828051,  0.00000000,  0.09571548, 0.00000000, 0.00000000],
#      [0.35828051, 1.03607456,  0.00000000,  0.40151695, 0.00000000, 0.00000000],
#      [0.00000000, 0.00000000,  1.07479494, -0.50700370, 0.00000000, 0.00000000],
#      [0.09571548, 0.40151695, -0.50700370,  0.93374109, 0.00000000, 0.00000000],
#      [0.00000000, 0.00000000,  0.00000000,  0.00000000, 1.07479494, 0.00000000],
#      [0.00000000, 0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.07479494],
#      ])

