import os
import shutil
import numpy as np
from mol import inp
from dalmisc.findif import *
import dalinp

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
    ev = ExpVal('XDIPLEN', mol=inp["ch2+"])
    res = ev.exe()
    assert_(res, 3.1348738989E-03)

def test_LR_HF(): 
    lr = LinResp('XANGMOM', 'X1SPNORB', mol=inp["ch2+"])
    res = lr.exe()
    assert_(res, -0.000131051802829)

def test_lr_HFx():
    lr = LinResp('XANGMOM', 'X1SPNORB', mol=inp["ch2+"], field='XDIPLEN', delta=0.0001)
    fd = FinDif(lr)
    res = fd.first()
    assert_(res, 0.00057501)

def test_QR_HF(): 
    qr = QuadResp('XANGMOM', 'X1SPNORB', 'XDIPLEN', wf="DFT\nGGAKEY HF=1", mol=inp["ch2+"])
    res = qr.exe()
    assert_(res, 0.00057501)

def test_LRx_QR():
        wf="HF"
        lr = FinDif(
            LinResp('XANGMOM', 'X1SPNORB', wf=wf, mol=inp["ch2+"], field='XDIPLEN', delta=0.0001)
            ).first()
        qr = QuadResp('XANGMOM', 'X1SPNORB', 'XDIPLEN', mol=inp["ch2+"], wf=wf).exe()
        assert_(lr, qr)
    

if __name__ == "__main__":
    pas
