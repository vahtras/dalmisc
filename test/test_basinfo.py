import os
from dalmisc.basinfo import BasInfo

def setup():
    global bas_info
    thisdir  = os.path.dirname(__file__)
    suppdir = os.path.join(thisdir, 'test_basinfo.d')
    bas_info = BasInfo(os.path.join(suppdir, 'SIRIUS.RST'))

def test_nsym(): 
    assert bas_info.nsym == 1
def test_nbas(): 
    assert bas_info.nbas[0] == 5
def test_nbast(): 
    assert bas_info.nbast == 5

if __name__ == "__main__":
    setup()
    test_nsym()
