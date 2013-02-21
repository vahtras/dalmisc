import os
from dalmisc.basinfo import BasInfo

def setup():
    global suppdir
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_basinfo.d')

def test_me(): 
    bi = BasInfo(os.path.join(suppdir, 'SIRIUS.RST'))
    assert bi.nsym == 1
    assert bi.nbas[0] == 5
    assert len(bi.nbas) == 8

if __name__ == "__main__":
    setup()
    test_me()
