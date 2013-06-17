import os
import shutil
import numpy as np
from mol import inp
from dalmisc.findif import *
#import dalinp

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

def test_E_HF():
    ecal = RspCalc(mol=inp["ch2+"])
    res = ecal.exe()
    assert_(res, -39.072083833208)

def test_EV_HF():
    ev = RspCalc('XDIPLEN', mol=inp["ch2+"])
    res = ev.exe()
    assert_(res, 3.1348738989E-03)

def test_EV_DFT():
    ev = RspCalc('XDIPLEN', wf='DFT\nB3LYP', dal='b3lyp', mol=inp["ch2+"])
    res = ev.exe()
    assert_(res, 3.18476890e-03)

def test_LR_HF(): 
    lr = RspCalc('XANGMOM', 'X1SPNORB', mol=inp["ch2+"])
    res = lr.exe()
    assert_(res, -0.000131051802829)

def test_LR_DFT(): 
    lr = RspCalc('XANGMOM', 'X1SPNORB', wf='DFT\nLDA', dal='lda', mol=inp["ch2+"])
    res = lr.exe()
    assert_(res, -1.770146614464e-04)

def test_lr_HFx():
    lr = RspCalc('XANGMOM', 'X1SPNORB', mol=inp["ch2+"], field='XDIPLEN', delta=0.0001)
    fd = FinDif(lr)
    res = fd.first()
    assert_(res, 0.00057501)

def test_QR_HF(): 
    qr = RspCalc('XANGMOM', 'X1SPNORB', 'XDIPLEN', wf="DFT\nGGAKEY HF=1", mol=inp["ch2+"], parallel=False)
    res = qr.exe()
    assert_(res, 0.00057501)

def test_LRx_QR():
        wf="HF"
        lr = FinDif(
            RspCalc('XANGMOM', 'X1SPNORB', wf=wf, mol=inp["ch2+"], field='XDIPLEN', delta=0.0001)
            ).first()
        qr = RspCalc('XANGMOM', 'X1SPNORB', 'XDIPLEN', mol=inp["ch2+"], wf=wf, parallel=False).exe()
        assert_(lr, qr)
    

if __name__ == "__main__":
    setup()
    test_EV_HF()
    teardown()
