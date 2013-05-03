import os

mol = """ATOMBASIS
STO-3G
------
    2   1 Z
        6.    1 Basis=STO-3G
C1      0.0006122714    0.0000000000    0.0000000000   
        1.    2 Basis=STO-3G
H1      1.5162556382   -1.3708721537    0.0000000000
H2     -0.7584339548    0.6854360769    1.7695110698
"""

class FinDif: 
    """Take objects with an exe method taking a 
    float argument returning a numerical object supporting subtraction
    and scalar multiplication allowing for finite differetiation"""

    def __init__(self, obj):
        self.obj = obj

    def first(self):
        d = self.obj.delta
        exe = self.obj.exe
        ret1 = exe(0.5*d)
        ret2 = exe(-0.5*d)
        return  (1/d) * (ret1 - ret2)

    def second(self):
        d = self.obj.delta
        exe = self.obj.exe
        ret1 = exe(0.5*d)
        ret2 = exe(-0.5*d)
        ret0 = exe(0.0)
        return  (4/d**2) * (ret1 - 2*ret0  + ret2)

class ExpVal:
    """Execute dalton LR"""
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.A, = args
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)

    def exe(self, delta=None):
        if self.field and delta:
            self.ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            self.ff = "###"
        dal = """**DALTON INPUT
.RUN RESPONSE
**WAVE FUNCTIONS
.%s
*SCF INPUT
.THRESHOLD
1e-12
%s
**RESPONSE
.PROPAV
%s
**END OF DALTON
"""%(self.wf, self.ff, self.A)

        dalfile = open('DALTON.INP', 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open('MOLECULE.INP', 'w')
        molfile.write(mol)
        molfile.close()

        os.system('rm -f RSPVEC RESULTS.RSP; dalton.x > log 2> err')

        result = None
        for line in open('DALTON.OUT'):
            if "total" in line and self.A in line:
                data = line.split(':')[1].replace('D', 'E')
                result = float(data)
                break
        return result

class LinResp:
    """Execute dalton LR"""
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.A, self.B = args
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)

    def exe(self, delta=None):
        if self.field and delta:
            self.ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            self.ff = "###"
        dal = """**DALTON INPUT
.RUN RESPONSE
**WAVE FUNCTIONS
.%s
*SCF INPUT
.THRESHOLD
1e-12
%s
**RESPONSE
*LINEAR
.THCLR
1e-9
.PROPRT
%s
.PROPRT
%s
**END OF DALTON
"""%(self.wf, self.ff, self.A, self.B)

        dalfile = open('DALTON.INP', 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open('MOLECULE.INP', 'w')
        molfile.write(mol)
        molfile.close()

        os.system('rm -f RSPVEC RESULTS.RSP; dalton.x > log 2> err')

        result = None
        for line in open('DALTON.OUT'):
            if "@" in line and self.A in line and self.B in line:
                data = line.split('=')[1].replace('D', 'E')
                result = -float(data)
                break
        return result

class QuadResp:
    """Execute dalton LR"""
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.A, self.B, self.C = args
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)

    def exe(self, delta=None):
        if self.field and delta:
            self.ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            self.ff = "###"
        dal = """**DALTON INPUT
.RUN RESPONSE
**WAVE FUNCTIONS
.%s
*SCF INPUT
.THRESHOLD
1e-12
%s
**RESPONSE
*QUADRATIC
.THCLR
1e-9
.APROP 
%s
.BPROP
%s
.CPROP
%s
**END OF DALTON
"""%(self.wf, self.ff, self.A, self.B, self.C)

        dalfile = open('DALTON.INP', 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open('MOLECULE.INP', 'w')
        molfile.write(mol)
        molfile.close()

        os.system('rm -f RSPVEC RESULTS.RSP; dalton.x > log 2> err')

        result = None
        for line in open('DALTON.OUT'):
            if "@omega" in line:
                data = line.split()[-1]
                result = float(data)
                break
        return result
        


