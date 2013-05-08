import re

def dalinp(*args, **kwargs):
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
    rsp, = args
    print rsp
    evmatch = re.match("< *(\w+) *>", rsp)
    lrmatch = re.match("<< *(\w+) *; *(\w+) *>>", rsp)
    qrmatch = re.match("<< *(\w+) *; *(\w+) *, *(\w+) *>>", rsp)

    if evmatch:
        A, =  evmatch.groups()
        print A
        inpfile = """**DALTON
**RESPONSE
.PROPAV
%s
**END OF INPUT
""" % A
        return inpfile

    if lrmatch:
        A, B =  lrmatch.groups()
        print A, B
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

    if qrmatch:
        A, B, C =  qrmatch.groups()
        print A, B, C
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

    return None

    
