import pytest
from dalmisc import dalinp

X = "XDIPLEN"
Y = "YDIPLEN"
Z = "ZDIPLEN"


@pytest.mark.parametrize(
    'operator, expected',
    [
        (X, f"PROPAV\n{X}"),
        (Y, f"PROPAV\n{Y}"),
        (Z, f"PROPAV\n{Z}"),
    ]
)
def test_ave(operator, expected):
    inp = dalinp.dalinp(operator)
    assert expected in inp


@pytest.mark.parametrize(
    'operators, expected',
    [
        ((X, Y), "PROPRT\nXDIPLEN\n.PROPRT\nYDIPLEN"),
        ((X, Z), "PROPRT\nXDIPLEN\n.PROPRT\nZDIPLEN"),
    ]
)
def test_lr(operators, expected):
    inp = dalinp.dalinp(*operators)
    assert expected in inp


@pytest.mark.parametrize(
    'operators, expected',
    [
        ((X, X, X), "APROP\nXDIPLEN\n.BPROP\nXDIPLEN\n.CPROP\nXDIPLEN"),
        ((X, Y, Y), "APROP\nXDIPLEN\n.BPROP\nYDIPLEN\n.CPROP\nYDIPLEN"),
        ((X, Y, Z), "APROP\nXDIPLEN\n.BPROP\nYDIPLEN\n.CPROP\nZDIPLEN"),
    ]
)
def test_qr(operators, expected):
    inp = dalinp.dalinp(*operators)
    assert expected in inp


def test_wf():
    wf = dalinp.wfinp()
    assert "WAVE FUNCTION\n.HF" in wf
