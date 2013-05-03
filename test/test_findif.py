import os
import shutil
import numpy as np
from dalmisc.findif import *

def assert_(this,ref):
    print this
    print ref
    assert np.allclose(this, ref)
        

def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"
    if os.path.isdir(suppdir):
        shutil.rmtree(suppdir)
    os.mkdir(suppdir)
    os.chdir(suppdir)

def teardown():
    #shutil.rmtree(suppdir)
    pass

def test_EV_HF():
    ev = ExpVal('XDIPLEN')
    res = ev.exe()
    assert_(res, 3.13487390e-3)

def test_LR_HF(): 
    lr = LinResp('XANGMOM', 'X1SPNORB')
    res = lr.exe()
    assert_(res, -1.310518028294e-04)

def test_lr_HFx():
    lr = LinResp('XANGMOM', 'X1SPNORB', field='XDIPLEN', delta=0.0001)
    fd = FinDif(lr)
    res = fd.first()
    assert_(res, 0.00057501)

def test_QR_HF(): 
    qr = QuadResp('XANGMOM', 'X1SPNORB', 'XDIPLEN')
    res = qr.exe()
    assert_(res, 0.00057501)

def test_LRx_QR():
        wf="HF"
        lr = FinDif(
            LinResp('XANGMOM', 'X1SPNORB', wf=wf, field='XDIPLEN', delta=0.0001)
            ).first()
        qr = QuadResp('XANGMOM', 'X1SPNORB', 'XDIPLEN', wf=wf).exe()
        assert_(lr, qr)
    

if __name__ == "__main__":
    pass
