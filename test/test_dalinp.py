from .. import rspinp
def setup():
    pass

def teardown():
    pass

def test_ave():
    ave = "<XDIPLEN>"
    inp = rspinp.dalinp(ave)
    assert "PROPAV\nXDIPLEN" in inp

def test_lr():
    ave = "<<XDIPLEN; YDIPLEN>>"
    inp = rspinp.dalinp(ave)
    assert "PROPRT\nXDIPLEN\n.PROPRT\nYDIPLEN" in inp

def test_qr():
    ave = "<<XDIPLEN; YDIPLEN, ZDIPLEN>>"
    inp = rspinp.dalinp(ave)
    assert "APROP\nXDIPLEN\n.BPROP\nYDIPLEN\n.CPROP\nZDIPLEN" in inp

def test_wf():
    wf = rspinp.wfinp()
    assert "WAVE FUNCTION\n.HF" in wf
