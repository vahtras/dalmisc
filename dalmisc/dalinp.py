
def dalinp(*args, **kwargs):
    """
    >>> print(dalinp("XDIPLEN"))
    **DALTON
    .RUN RESPONSE
    **WAVE FUNCTION
    .HF
    **DALTON
    **RESPONSE
    .PROPAV
    XDIPLEN
    **END OF INPUT
    """
    return wfinp(*args, **kwargs) + rspinp(*args, **kwargs)

def wfinp(*args, **kwargs):
    wf = kwargs.get('wf', 'HF')
    wfout  = """**DALTON
.RUN RESPONSE
**WAVE FUNCTION
.%s
""" % wf
    scfinp = kwargs.get('scfinp')
    if scfinp:
        wfout += "*SCF INPUT\n" 
        for k in scfinp:
            wfout += ".%s\n"%k
            if scfinp[k]:
                wfout += "%s\n" % scfinp[k]
    return wfout

def rspinp(*args, **kwargs):

    if len(args) == 1:
        A, = args
        inpfile = """**DALTON
**RESPONSE
.PROPAV
%s
**END OF INPUT
""" % A
        return inpfile

    elif len(args) == 2:
        A, B = args
        inpfile = """**DALTON
**RESPONSE
*LINEAR
.PROPRT
%s
.PROPRT
%s
**END OF INPUT
""" % (A, B)
        return inpfile

    elif len(args) == 3:
        A, B, C =  args
        inpfile = """**DALTON
**RESPONSE
*QUADRATIC
.APROP
%s
.BPROP
%s
.CPROP
%s
**END OF INPUT
""" % (A, B, C)
        return inpfile
    else:
        raise Exception("Input error")

    return None

class DalExe:
    """General Dalton executer"""    
    def __init__(self, *args, **kwargs):
        self.wf = kwargs.get('wf', 'HF')
        self.dal = dalinp(*args, **kwargs)
        self.field = kwargs.get('field', None)
        self.delta = kwargs.get('delta', 0)
        self.mol = kwargs.get('mol', None)

    def exe(self, delta=None):
        if self.field and delta:
            self.ff = "*HAMILTON\n.FIELD\n%f\n%s"%(delta, self.field)
        else:
            self.ff = "###"

        wf = self.wf.split('\n')[-1].replace(' ','_')
        dalfile = open(wf + ".dal", 'w')
        dalfile.write(self.dal)
        dalfile.close()
    
        molfile = open(wf + ".mol", 'w')
        molfile.write(self.mol)
        molfile.close()

        os.system("dalton -d -t /tmp/ExpVal_%s -N 4 %s > log 2>&1 " % (wf, wf))

        result = None
        for line in open(wf + ".out"):
            if "total" in line and self.A in line:
                data = line.split(':')[1].replace('D', 'E')
                result = float(data)
                break
        assert result is not None
        return result
