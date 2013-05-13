import os
import shutil
import numpy as np
from dalmisc.findif import *
from molinp import *

mol=h2op

def assert_(this,ref):
    print this, ref
    assert np.allclose(this, ref, rtol=1e-4)
        

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



def test_EVx_LR_HF():
    wf='HF'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)

def test_EVx_LR_Becke():
    wf='DFT\nBecke'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mBecke():
    wf='DFT\nmBecke'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B86x():
    wf='DFT\nB86x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B86mx():
    wf='DFT\nB86mx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B97():
    wf='DFT\nB97'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B97_1():
    wf='DFT\nB97-1'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B97_2():
    wf='DFT\nB97-2'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B97_3():
    wf='DFT\nB97-3'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B97_K():
    wf='DFT\nB97-K'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_Example2():
    wf='DFT\nExample2'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_Example():
    wf='DFT\nExample'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_DK87x():
    wf='DFT\nDK87x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_G96x():
    wf='DFT\nG96x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_KTx():
    wf='DFT\nKTx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LB94():
    wf='DFT\nLB94'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LG93x():
    wf='DFT\nLG93x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LRC95x():
    wf='DFT\nLRC95x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LYP():
    wf='DFT\nLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LYPr():
    wf='DFT\nLYPr'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH():
    wf='DFT\nHCTH'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH93():
    wf='DFT\nHCTH93'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH93m():
    wf='DFT\nHCTH93m'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH120():
    wf='DFT\nHCTH120'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH147():
    wf='DFT\nHCTH147'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH407():
    wf='DFT\nHCTH407'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_HCTH407p():
    wf='DFT\nHCTH407p'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_OPTX():
    wf='DFT\nOPTX'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPWx():
    wf='DFT\nmPWx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_P86c():
    wf='DFT\nP86c'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW86x():
    wf='DFT\nPW86x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91c():
    wf='DFT\nPW91c'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91nc():
    wf='DFT\nPW91nc'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91x():
    wf='DFT\nPW91x'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91x2():
    wf='DFT\nPW91x2'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW92c():
    wf='DFT\nPW92c'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW92ac():
    wf='DFT\nPW92ac'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PZ81():
    wf='DFT\nPZ81'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBEC():
    wf='DFT\nPBEC'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBEx():
    wf='DFT\nPBEx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_RPBEx():
    wf='DFT\nRPBEx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPBEx():
    wf='DFT\nmPBEx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_revPBEx():
    wf='DFT\nrevPBEx'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_Slater():
    wf='DFT\nSlater'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_VWN3():
    wf='DFT\nVWN3'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_VWN5():
    wf='DFT\nVWN5'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_VWNI():
    wf='DFT\nVWNI'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_VWN3I():
    wf='DFT\nVWN3I'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_VWN():
    wf='DFT\nVWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_XAlpha():
    wf='DFT\nXAlpha'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_Wigner():
    wf='DFT\nWigner'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_WL90c():
    wf='DFT\nWL90c'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B2PLYP():
    wf='DFT\nB2PLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B2TPLYP():
    wf='DFT\nB2TPLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW2PLYP():
    wf='DFT\nmPW2PLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW2KPLYP():
    wf='DFT\nmPW2KPLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B2GPPLYP():
    wf='DFT\nB2GPPLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B2PIPLYP():
    wf='DFT\nB2PIPLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBE0DH():
    wf='DFT\nPBE0DH'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B3LYP():
    wf='DFT\nB3LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B3LYPg():
    wf='DFT\nB3LYPg'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B3LYPGauss():
    wf='DFT\nB3LYPGauss'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B3P86():
    wf='DFT\nB3P86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B3P86g():
    wf='DFT\nB3P86g'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B3PW91():
    wf='DFT\nB3PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B1LYP():
    wf='DFT\nB1LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B1PW91():
    wf='DFT\nB1PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BHandH():
    wf='DFT\nBHandH'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BHandHLYP():
    wf='DFT\nBHandHLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B86VWN():
    wf='DFT\nB86VWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B86LYP():
    wf='DFT\nB86LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B86P86():
    wf='DFT\nB86P86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B86PW91():
    wf='DFT\nB86PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_B97_D():
    wf='DFT\nB97-D'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BVWN():
    wf='DFT\nBVWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BLYP():
    wf='DFT\nBLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BP86():
    wf='DFT\nBP86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BPW91():
    wf='DFT\nBPW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BW():
    wf='DFT\nBW'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_BFW():
    wf='DFT\nBFW'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_Camb3lyp():
    wf='DFT\nCamb3lyp'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_Combine():
    wf='DFT\nCombine'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_DBLYP():
    wf='DFT\nDBLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_DBP86():
    wf='DFT\nDBP86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_DBPW91():
    wf='DFT\nDBPW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_EDF1():
    wf='DFT\nEDF1'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_EDF2():
    wf='DFT\nEDF2'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_1_4():
    wf='DFT\n1/4'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)




def test_EVx_LR_G96VWN():
    wf='DFT\nG96VWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_G96LYP():
    wf='DFT\nG96LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_G96P86():
    wf='DFT\nG96P86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_G96PW91():
    wf='DFT\nG96PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_G961LYP():
    wf='DFT\nG961LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_KMLYP():
    wf='DFT\nKMLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_KT1():
    wf='DFT\nKT1'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_KT2():
    wf='DFT\nKT2'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_KT3():
    wf='DFT\nKT3'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LDA():
    wf='DFT\nLDA'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_LG1LYP():
    wf='DFT\nLG1LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPWVWN():
    wf='DFT\nmPWVWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPWLYP():
    wf='DFT\nmPWLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPWP86():
    wf='DFT\nmPWP86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPWPW91():
    wf='DFT\nmPWPW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW91():
    wf='DFT\nmPW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW1PW91():
    wf='DFT\nmPW1PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW3PW91():
    wf='DFT\nmPW3PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW1K():
    wf='DFT\nmPW1K'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW1N():
    wf='DFT\nmPW1N'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPW1S():
    wf='DFT\nmPW1S'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_OVWN():
    wf='DFT\nOVWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_OLYP():
    wf='DFT\nOLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_OP86():
    wf='DFT\nOP86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_OPW91():
    wf='DFT\nOPW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBE0():
    wf='DFT\nPBE0'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBE0PBE():
    wf='DFT\nPBE0PBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBE1PBE():
    wf='DFT\nPBE1PBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBE():
    wf='DFT\nPBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PBEPBE():
    wf='DFT\nPBEPBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_rCamb3lyp():
    wf='DFT\nrCamb3lyp'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_revPBE():
    wf='DFT\nrevPBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_RPBE():
    wf='DFT\nRPBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_mPBE():
    wf='DFT\nmPBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91():
    wf='DFT\nPW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_revPBE():
    wf='DFT\nrevPBE'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91VWN():
    wf='DFT\nPW91VWN'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91LYP():
    wf='DFT\nPW91LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91P86():
    wf='DFT\nPW91P86'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_PW91PW91():
    wf='DFT\nPW91PW91'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_SVWN3():
    wf='DFT\nSVWN3'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_SVWN5():
    wf='DFT\nSVWN5'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_XLYP():
    wf='DFT\nXLYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


def test_EVx_LR_X3LYP():
    wf='DFT\nX3LYP'
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)

if __name__ == "__main__":
    import sys
    setup()
    eval("test_EVx_LR_%s()"%sys.argv[1])
    teardown()
