import os
import dalinp


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
        #print ret1, ret2
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
        self.mol = kwargs.get('mol', None)
        self.triplet = kwargs.get('triplet', False)

    def exe(self, delta=None):

        if self.field and delta:
            ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            ff = "###"

        if self.triplet:
            trpflg = ".TRPFLG"
        else:
            trpflg = "#"

        dal = """**DALTON INPUT
.RUN RESPONSE
**WAVE FUNCTIONS
.%s
*SCF INPUT
.THRESHOLD
1e-12
%s
**RESPONSE
%s
.PROPAV
%s
**END OF DALTON
"""%(self.wf, ff, trpflg, self.A)

        wf = self.wf.split('\n')[-1].replace(' ','_').replace('/','_')
        dalfile = open(wf + ".dal", 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open(wf + ".mol", 'w')
        molfile.write(self.mol)
        molfile.close()

        os.system("dalton -N 8 -d -t /tmp/ExpVal_%s %s > log 2>&1 " % (wf, wf))

        result = None
        print "###",self.A.split()
        A  = self.A.split()[0]
        for line in open(wf + ".out"):
            if "total" in line and A in line:
                data = line.split(':')[1].replace('D', 'E')
                result = float(data)
                break
        assert result is not None
        return result

class LinResp:
    """Execute dalton LR"""
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.A, self.B = args
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)
        self.mol = kwargs.get('mol')
        self.triplet = kwargs.get('triplet', False)
        self.aux = kwargs.get('aux', '#')

    def exe(self, delta=None):
        if self.field and delta:
            ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            ff = "###"

        if self.triplet:
            trpflg = ".TRPFLG"
        else:
            trpflg = "#"

        dal = """**DALTON INPUT
.RUN RESPONSE
**WAVE FUNCTIONS
.%s
*SCF INPUT
.THRESHOLD
1e-12
%s
**RESPONSE
%s
*LINEAR
.THCLR
1e-9
.PROPRT
%s
.PROPRT
%s
%s
**END OF DALTON
"""%(self.wf, ff, trpflg, self.A, self.B, self.aux)

        wf = self.wf.split('\n')[-1].replace(' ','_').replace('/','_')
        dalfile = open(wf + ".dal", 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open(wf + ".mol", 'w')
        molfile.write(self.mol)
        molfile.close()

        os.system("dalton -N 8 -d -t /tmp/LinResp_%s  %s > log 2>&1 " % (wf, wf))

        result = None
        A = self.A.split()[0]
        B = self.B.split()[0]
        for line in open(wf + ".out"):
            if "@" in line and A in line and B in line:
                data = line.split('=')[1].replace('D', 'E')
                result = -float(data)
                break
        assert result is not None
        return result

class QuadResp:
    """Execute dalton QR"""
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.A, self.B, self.C = args
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)
        self.mol = kwargs.get('mol')
        self.triplet = kwargs.get('triplet', False)
        self.aux = kwargs.get('aux', '#')

    def exe(self, delta=None):
        if self.field and delta:
            ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            ff = "###"

        if self.triplet:
            trpflg = ".TRPFLG"
        else:
            trpflg = "#"

        dal = """**DALTON INPUT
.RUN RESPONSE
**WAVE FUNCTIONS
.%s
*SCF INPUT
.THRESHOLD
1e-12
%s
**RESPONSE
%s
*QUADRATIC
.THCLR
1e-9
.APROP 
%s
.BPROP
%s
.CPROP
%s
%s
**END OF DALTON
"""%(self.wf, ff, trpflg, self.A, self.B, self.C, self.aux)

        wf = self.wf.split('\n')[-1].replace(' ','_').replace('/','_')
        dalfile = open(wf + ".dal", 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open(wf + ".mol", 'w')
        molfile.write(self.mol)
        molfile.close()

        os.system("dalton -N 8 -d -t /tmp/QuadResp_%s  %s > log 2>&1 " % (wf, wf))

        result = None
        for line in open(wf + ".out"):
            if "@omega" in line:
                data = line.split()[-1]
                result = float(data)
                break
        assert result is not None
        return result

class CubResp:
    """Execute dalton CR"""
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.A, self.B, self.C, self.D = args
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)
        self.mol = kwargs.get('mol')
        self.trpflg = kwargs.get('trpflg', '#')
        self.aux = kwargs.get('aux', '#')

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
%s
*CUBIC
.THCLR
1e-9
.APROP 
%s
.BPROP
%s
.CPROP
%s
.DPROP
%s
%s
**END OF DALTON
"""%(self.wf, self.ff, self.trpflg, self.A, self.B, self.C, self.D, self.aux)

        wf = self.wf.split('\n')[-1].replace(' ','_').replace('/','_')
        dalfile = open(wf + ".dal", 'w')
        dalfile.write(dal)
        dalfile.close()
    
        molfile = open(wf + ".mol", 'w')
        molfile.write(self.mol)
        molfile.close()

        os.system("dalton -N 8 -d -t /tmp/QuadResp_%s  %s > log 2>&1 " % (wf, wf))

        result = None
        for line in open(wf + ".out"):
            if "@ << A; B, C, D >>" in line:
                data = line.split()[-1]
                result = float(data)
                break
        assert result is not None
        return result
        


