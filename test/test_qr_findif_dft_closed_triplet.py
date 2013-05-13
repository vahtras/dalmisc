import os
import shutil
import numpy as np
from dalmisc.findif import *
from molinp import *

mol=h2o

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



def test_LRx_QR_HF():
    wf='HF'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HF():
    wf='HF'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Becke():
    wf='DFT\nBecke'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mBecke():
    wf='DFT\nmBecke'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B86x():
    wf='DFT\nB86x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B86mx():
    wf='DFT\nB86mx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B97():
    wf='DFT\nB97'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B97_1():
    wf='DFT\nB97-1'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B97_2():
    wf='DFT\nB97-2'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B97_3():
    wf='DFT\nB97-3'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B97_K():
    wf='DFT\nB97-K'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Example2():
    wf='DFT\nExample2'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Example():
    wf='DFT\nExample'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_DK87x():
    wf='DFT\nDK87x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_G96x():
    wf='DFT\nG96x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_KTx():
    wf='DFT\nKTx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LB94():
    wf='DFT\nLB94'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LG93x():
    wf='DFT\nLG93x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LRC95x():
    wf='DFT\nLRC95x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LYP():
    wf='DFT\nLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LYPr():
    wf='DFT\nLYPr'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH():
    wf='DFT\nHCTH'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH93():
    wf='DFT\nHCTH93'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH93m():
    wf='DFT\nHCTH93m'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH120():
    wf='DFT\nHCTH120'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH147():
    wf='DFT\nHCTH147'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH407():
    wf='DFT\nHCTH407'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_HCTH407p():
    wf='DFT\nHCTH407p'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_OPTX():
    wf='DFT\nOPTX'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPWx():
    wf='DFT\nmPWx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_P86c():
    wf='DFT\nP86c'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW86x():
    wf='DFT\nPW86x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91c():
    wf='DFT\nPW91c'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91nc():
    wf='DFT\nPW91nc'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91x():
    wf='DFT\nPW91x'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91x2():
    wf='DFT\nPW91x2'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW92c():
    wf='DFT\nPW92c'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW92ac():
    wf='DFT\nPW92ac'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PZ81():
    wf='DFT\nPZ81'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBEC():
    wf='DFT\nPBEC'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBEx():
    wf='DFT\nPBEx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_RPBEx():
    wf='DFT\nRPBEx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPBEx():
    wf='DFT\nmPBEx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_revPBEx():
    wf='DFT\nrevPBEx'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Slater():
    wf='DFT\nSlater'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_VWN3():
    wf='DFT\nVWN3'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_VWN5():
    wf='DFT\nVWN5'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_VWNI():
    wf='DFT\nVWNI'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_VWN3I():
    wf='DFT\nVWN3I'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_VWN():
    wf='DFT\nVWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_XAlpha():
    wf='DFT\nXAlpha'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Wigner():
    wf='DFT\nWigner'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_WL90c():
    wf='DFT\nWL90c'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B2PLYP():
    wf='DFT\nB2PLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B2TPLYP():
    wf='DFT\nB2TPLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW2PLYP():
    wf='DFT\nmPW2PLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW2KPLYP():
    wf='DFT\nmPW2KPLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B2GPPLYP():
    wf='DFT\nB2GPPLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B2PIPLYP():
    wf='DFT\nB2PIPLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBE0DH():
    wf='DFT\nPBE0DH'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B3LYP():
    wf='DFT\nB3LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B3LYPg():
    wf='DFT\nB3LYPg'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B3LYPGauss():
    wf='DFT\nB3LYPGauss'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B3P86():
    wf='DFT\nB3P86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B3P86g():
    wf='DFT\nB3P86g'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B3PW91():
    wf='DFT\nB3PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B1LYP():
    wf='DFT\nB1LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B1PW91():
    wf='DFT\nB1PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BHandH():
    wf='DFT\nBHandH'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BHandHLYP():
    wf='DFT\nBHandHLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B86VWN():
    wf='DFT\nB86VWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B86LYP():
    wf='DFT\nB86LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B86P86():
    wf='DFT\nB86P86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B86PW91():
    wf='DFT\nB86PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_B97_D():
    wf='DFT\nB97-D'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BVWN():
    wf='DFT\nBVWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BLYP():
    wf='DFT\nBLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BP86():
    wf='DFT\nBP86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BPW91():
    wf='DFT\nBPW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BW():
    wf='DFT\nBW'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_BFW():
    wf='DFT\nBFW'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Camb3lyp():
    wf='DFT\nCamb3lyp'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_Combine():
    wf='DFT\nCombine'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_DBLYP():
    wf='DFT\nDBLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_DBP86():
    wf='DFT\nDBP86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_DBPW91():
    wf='DFT\nDBPW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_EDF1():
    wf='DFT\nEDF1'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_EDF2():
    wf='DFT\nEDF2'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_1_4():
    wf='DFT\n1/4'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_G96VWN():
    wf='DFT\nG96VWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_G96LYP():
    wf='DFT\nG96LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_G96P86():
    wf='DFT\nG96P86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_G96PW91():
    wf='DFT\nG96PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_G961LYP():
    wf='DFT\nG961LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_KMLYP():
    wf='DFT\nKMLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_KT1():
    wf='DFT\nKT1'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_KT2():
    wf='DFT\nKT2'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_KT3():
    wf='DFT\nKT3'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LDA():
    wf='DFT\nLDA'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_LG1LYP():
    wf='DFT\nLG1LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPWVWN():
    wf='DFT\nmPWVWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPWLYP():
    wf='DFT\nmPWLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPWP86():
    wf='DFT\nmPWP86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPWPW91():
    wf='DFT\nmPWPW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW91():
    wf='DFT\nmPW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW1PW91():
    wf='DFT\nmPW1PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW3PW91():
    wf='DFT\nmPW3PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW1K():
    wf='DFT\nmPW1K'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW1N():
    wf='DFT\nmPW1N'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPW1S():
    wf='DFT\nmPW1S'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_OVWN():
    wf='DFT\nOVWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_OLYP():
    wf='DFT\nOLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_OP86():
    wf='DFT\nOP86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_OPW91():
    wf='DFT\nOPW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBE0():
    wf='DFT\nPBE0'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBE0PBE():
    wf='DFT\nPBE0PBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBE1PBE():
    wf='DFT\nPBE1PBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBE():
    wf='DFT\nPBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PBEPBE():
    wf='DFT\nPBEPBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_rCamb3lyp():
    wf='DFT\nrCamb3lyp'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_RPBE():
    wf='DFT\nRPBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_mPBE():
    wf='DFT\nmPBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91():
    wf='DFT\nPW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_revPBE():
    wf='DFT\nrevPBE'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91VWN():
    wf='DFT\nPW91VWN'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91LYP():
    wf='DFT\nPW91LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91P86():
    wf='DFT\nPW91P86'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_PW91PW91():
    wf='DFT\nPW91PW91'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_SVWN3():
    wf='DFT\nSVWN3'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_SVWN5():
    wf='DFT\nSVWN5'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_XLYP():
    wf='DFT\nXLYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


def test_LRx_QR_X3LYP():
    wf='DFT\nX3LYP'
    lr = FinDif(
        LinResp('ZANGMOM 1', 'ZANGMOM 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001, triplet=True)
        ).first()
    qr = QuadResp('ZANGMOM 1', 'ZANGMOM 1', 'XDIPLEN', wf=wf, mol=mol, triplet=True, aux=".ISPABC\n1 1 0").exe()
    assert_(lr, qr)


if __name__ == "__main__":
    import sys
    setup()
    wf="DFT\n%s"%sys.argv[1]
    ev = FinDif(
        ExpVal('XDIPLEN', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPLEN', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)

