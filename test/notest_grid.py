import os
import math
import numpy as np
from dalmisc import grid

def assert_(this,ref):
    print this
    print ref
    assert np.allclose(this, ref)

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"
    os.chdir(suppdir)

def teardown():
    pass

def test_readfirst():
    wRref = [
        -8.575690271798617, 
        -4.223255559902539, 
        -2.8653033453872707, 
        3.3506109330112563,
        ]
    wR = grid.readfirst('DALTON.QUAD')
    assert_(wR, wRref)

def test_TI_readfirst():
    wRref = [
        -8.575690271798617, 
        -4.223255559902539, 
        -2.8653033453872707, 
        3.3506109330112563,
        ]
    gr = grid.GridIntegrator('DALTON.QUAD')
    wR = gr.points[0,:]
    assert_(wR, wRref)

def test_generate_primitive_s():
    e = [1.0, 0.5]
    x = (1, 1, 1)
    X = (0, 0, 0)
    prim =  grid.generate_primitive(0, e, X)
    this,_ = prim(*x)
    ref =  [(2*a/math.pi)**0.75*math.exp(-a*3) for a in e]
    assert_(this, ref)

def test_integrate_primitive_s():
    e = [1.0 ]
    prim = grid.generate_primitive(0, e, (0,0,0))
    s = 0
    for x, y, z, w in grid.quaditer('DALTON.QUAD'):
        p, _ = prim(x, y, z)
        s += w*p[0]**2
    assert_(s, 1.0)


def test_generate_primitive_p():
    e = [1.0, 0.5]
    x = (1, 1, 1)
    X = (0, 0, 0)
    prim =  grid.generate_primitive(1, e, X)
    this, _  = prim(*x)
    ref =  [(2*a/math.pi)**0.75*math.sqrt(4*a)*c*math.exp(-a*3) for a in e for c in x]
    assert_(this, ref)

def notest_generate_contracted_s():
    c = [.1, .2]
    e = [1.0, 0.5]
    N = [(2*e_/math.pi)**0.75 for e_ in e]
    x = (1, 1, 1)
    X = (0, 0, 0)
    cont = grid.generate_contracted(0, e, c, X)
    ref =  0.03724419688727937
    assert_(cont(*x), ref)

def test_integrate_contracted():
    c = [.1, .2]
    e = [1.0, 0.5]
    cont = grid.generate_contracted(0, e, c, (0,0,0))
    s = 0
    for x, y, z, w in grid.quaditer('DALTON.QUAD'):
        f, _ = cont(x, y, z)
        s += w*f[0]**2
    assert_(s, 1.0)

def notest_generate_p():
    c = [.1, .2]
    e = [1.0, 0.5]
    N = [(2*e_/math.pi)**0.75 * (4*e_)**0.5 for e_ in e]
    x = (1, 1, 1)
    X = (0, 0, 0)
    p = grid.generate_contracted(1, e, c, X)
    px, py, pz = p(*x)
    refx = 0.453197088787
    assert_(px, refx)

def notest_generate_d():
    c = [.1, .2]
    e = [1.0, 0.5]
    N = [(2*e_/math.pi)**0.75 * ((8*e_)**2/12)**0.5 for e_ in e]
    x = (1, 1, 1)
    X = (0, 0, 0)
    d = grid.generate_contracted(2, e, c, X)
    dxx = d(*x)[0]
    refxx = 0.378386464369
    assert_(dxx, refxx)


def test_set_cfunc():
    cfunc = grid.set_cfunc('DALTON.BAS')
    assert_(len(cfunc),  16)

def test_eval_cfunc():
    cfunc = grid.set_cfunc('DALTON.BAS')
    res = grid.eval_cfunc(cfunc, (0,0,0))
    #print res
    assert_(len(res), 13)
    
    
def test_gendens():
    #rho = gendens(x, y, z)
    pass

def notest_electrons():
    wRref = 10.0
    wR = grid.count_electrons('DALTON.BAS', 'DALTON.QUAD')
    assert_(wR, wRref)
    pass

def test_integrator_electrons():
    wRref = (5.0, 5.0)
    def rhoa(ra, rb):
        return ra
    def rhob(ra, rb):
        return rb
    args = (rhoa, rhob)
    kwargs = {'bf':'DALTON.BAS', 'qf':'DALTON.QUAD'}
    wR = grid.integrator(*args, **kwargs)
    assert_(wR, wRref)

def test_energy_slater():
    PREF = -0.75*(6.0/math.pi)**(1.0/3.0)
    F43 = 4.0/3.0
    def slater(rhoa, rhob):
        return PREF*(rhoa**F43 + rhob**F43)

    args = (slater,)
    kwargs = {'bf':'DALTON.BAS', 'qf':'DALTON.QUAD'}
    ESlater = grid.integrator(*args, **kwargs)
    Eref = -8.080707409149
    assert_(ESlater, Eref)
    
    

if __name__ == "__main__":
    from profile import run
    setup()
    #run('test_integrator_electrons()')
    #import pdb; pdb.set_trace()
    run('test_GI_integrate_exp()')
    teardown
    
