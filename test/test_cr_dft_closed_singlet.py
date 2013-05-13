
import os
import shutil
import numpy as np
from dalmisc.findif import *
from mol import inp

mol=inp["h2o"]
#mol=inp["h2o+"]

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


def test_QRx_CR_HF():
    wf='HF'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Becke():
    wf='DFT\nBecke'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mBecke():
    wf='DFT\nmBecke'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B86x():
    wf='DFT\nB86x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B86mx():
    wf='DFT\nB86mx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B97():
    wf='DFT\nB97'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B97_1():
    wf='DFT\nB97-1'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B97_2():
    wf='DFT\nB97-2'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B97_3():
    wf='DFT\nB97-3'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B97_K():
    wf='DFT\nB97-K'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Example2():
    wf='DFT\nExample2'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Example():
    wf='DFT\nExample'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_DK87x():
    wf='DFT\nDK87x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_G96x():
    wf='DFT\nG96x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_KTx():
    wf='DFT\nKTx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LB94():
    wf='DFT\nLB94'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LG93x():
    wf='DFT\nLG93x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LRC95x():
    wf='DFT\nLRC95x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LYP():
    wf='DFT\nLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LYPr():
    wf='DFT\nLYPr'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH():
    wf='DFT\nHCTH'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH93():
    wf='DFT\nHCTH93'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH93m():
    wf='DFT\nHCTH93m'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH120():
    wf='DFT\nHCTH120'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH147():
    wf='DFT\nHCTH147'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH407():
    wf='DFT\nHCTH407'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_HCTH407p():
    wf='DFT\nHCTH407p'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_OPTX():
    wf='DFT\nOPTX'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPWx():
    wf='DFT\nmPWx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_P86c():
    wf='DFT\nP86c'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW86x():
    wf='DFT\nPW86x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91c():
    wf='DFT\nPW91c'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91nc():
    wf='DFT\nPW91nc'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91x():
    wf='DFT\nPW91x'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91x2():
    wf='DFT\nPW91x2'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW92c():
    wf='DFT\nPW92c'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW92ac():
    wf='DFT\nPW92ac'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PZ81():
    wf='DFT\nPZ81'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBEC():
    wf='DFT\nPBEC'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBEx():
    wf='DFT\nPBEx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_RPBEx():
    wf='DFT\nRPBEx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPBEx():
    wf='DFT\nmPBEx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_revPBEx():
    wf='DFT\nrevPBEx'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Slater():
    wf='DFT\nSlater'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_VWN3():
    wf='DFT\nVWN3'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_VWN5():
    wf='DFT\nVWN5'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_VWNI():
    wf='DFT\nVWNI'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_VWN3I():
    wf='DFT\nVWN3I'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_VWN():
    wf='DFT\nVWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_XAlpha():
    wf='DFT\nXAlpha'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Wigner():
    wf='DFT\nWigner'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_WL90c():
    wf='DFT\nWL90c'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B2PLYP():
    wf='DFT\nB2PLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B2TPLYP():
    wf='DFT\nB2TPLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW2PLYP():
    wf='DFT\nmPW2PLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW2KPLYP():
    wf='DFT\nmPW2KPLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B2GPPLYP():
    wf='DFT\nB2GPPLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B2PIPLYP():
    wf='DFT\nB2PIPLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBE0DH():
    wf='DFT\nPBE0DH'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B3LYP():
    wf='DFT\nB3LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B3LYPg():
    wf='DFT\nB3LYPg'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B3LYPGauss():
    wf='DFT\nB3LYPGauss'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B3P86():
    wf='DFT\nB3P86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B3P86g():
    wf='DFT\nB3P86g'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B3PW91():
    wf='DFT\nB3PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B1LYP():
    wf='DFT\nB1LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B1PW91():
    wf='DFT\nB1PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BHandH():
    wf='DFT\nBHandH'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BHandHLYP():
    wf='DFT\nBHandHLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B86VWN():
    wf='DFT\nB86VWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B86LYP():
    wf='DFT\nB86LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B86P86():
    wf='DFT\nB86P86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B86PW91():
    wf='DFT\nB86PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_B97_D():
    wf='DFT\nB97-D'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BVWN():
    wf='DFT\nBVWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BLYP():
    wf='DFT\nBLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BP86():
    wf='DFT\nBP86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BPW91():
    wf='DFT\nBPW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BW():
    wf='DFT\nBW'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_BFW():
    wf='DFT\nBFW'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Camb3lyp():
    wf='DFT\nCamb3lyp'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_Combine():
    wf='DFT\nCombine'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_DBLYP():
    wf='DFT\nDBLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_DBP86():
    wf='DFT\nDBP86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_DBPW91():
    wf='DFT\nDBPW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_EDF1():
    wf='DFT\nEDF1'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_EDF2():
    wf='DFT\nEDF2'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_1_4():
    wf='DFT\n1/4'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_G96VWN():
    wf='DFT\nG96VWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_G96LYP():
    wf='DFT\nG96LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_G96P86():
    wf='DFT\nG96P86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_G96PW91():
    wf='DFT\nG96PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_G961LYP():
    wf='DFT\nG961LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_KMLYP():
    wf='DFT\nKMLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_KT1():
    wf='DFT\nKT1'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_KT2():
    wf='DFT\nKT2'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_KT3():
    wf='DFT\nKT3'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LDA():
    wf='DFT\nLDA'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_LG1LYP():
    wf='DFT\nLG1LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPWVWN():
    wf='DFT\nmPWVWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPWLYP():
    wf='DFT\nmPWLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPWP86():
    wf='DFT\nmPWP86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPWPW91():
    wf='DFT\nmPWPW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW91():
    wf='DFT\nmPW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW1PW91():
    wf='DFT\nmPW1PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW3PW91():
    wf='DFT\nmPW3PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW1K():
    wf='DFT\nmPW1K'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW1N():
    wf='DFT\nmPW1N'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPW1S():
    wf='DFT\nmPW1S'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_OVWN():
    wf='DFT\nOVWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_OLYP():
    wf='DFT\nOLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_OP86():
    wf='DFT\nOP86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_OPW91():
    wf='DFT\nOPW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBE0():
    wf='DFT\nPBE0'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBE0PBE():
    wf='DFT\nPBE0PBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBE1PBE():
    wf='DFT\nPBE1PBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBE():
    wf='DFT\nPBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PBEPBE():
    wf='DFT\nPBEPBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_rCamb3lyp():
    wf='DFT\nrCamb3lyp'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_RPBE():
    wf='DFT\nRPBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_mPBE():
    wf='DFT\nmPBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91():
    wf='DFT\nPW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_revPBE():
    wf='DFT\nrevPBE'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91VWN():
    wf='DFT\nPW91VWN'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91LYP():
    wf='DFT\nPW91LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91P86():
    wf='DFT\nPW91P86'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_PW91PW91():
    wf='DFT\nPW91PW91'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_SVWN3():
    wf='DFT\nSVWN3'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_SVWN5():
    wf='DFT\nSVWN5'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_XLYP():
    wf='DFT\nXLYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


def test_QRx_CR_X3LYP():
    wf='DFT\nX3LYP'
    qr = FinDif(
        QuadResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', wf=wf, mol=mol, field='YDIPLEN', delta=0.0001, triplet=True)
        ).first()
    cr = CubResp('ZANGMOM', 'ZANGMOM', 'YDIPLEN', 'YDIPLEN', wf=wf, mol=mol).exe()
    assert_(qr, cr)


if __name__ == "__main__":
    import sys
    setup()
    #eval("test_EVx_LR_%s()"%sys.argv[1])
    #eval("test_LRx_QR_%s()"%sys.argv[1])
    eval("test_QRx_CR_%s()"%sys.argv[1])
    teardown()

