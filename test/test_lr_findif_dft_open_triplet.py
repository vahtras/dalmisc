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
        ExpVal('XDIPVEL 1', wf=wf, mol=mol, field='XDIPLEN', delta=0.0001)
        ).first()
    lr = LinResp('XDIPVEL 1', 'XDIPLEN', wf=wf, mol=mol).exe()
    assert_(ev, lr)


if __name__ == "__main__":
    import sys
    setup()
    eval("test_EVx_LR_%s()"%sys.argv[1])
    teardown()
