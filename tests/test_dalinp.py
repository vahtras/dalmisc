from dalmisc import dalinp

X = "XDIPLEN"
Y = "YDIPLEN"
Z = "ZDIPLEN"


def setup():
    pass


def teardown():
    pass


def test_ave():
    inp = dalinp.dalinp(X)
    assert "PROPAV\nXDIPLEN" in inp


def test_lr():
    inp = dalinp.dalinp(X, Y)
    assert "PROPRT\nXDIPLEN\n.PROPRT\nYDIPLEN" in inp


def test_qr():
    inp = dalinp.dalinp(X, Y, Z)
    assert "APROP\nXDIPLEN\n.BPROP\nYDIPLEN\n.CPROP\nZDIPLEN" in inp


def test_wf():
    wf = dalinp.wfinp()
    assert "WAVE FUNCTION\n.HF" in wf
