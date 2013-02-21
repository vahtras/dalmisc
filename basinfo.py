#!/usr/bin/env python
"""Extract info from DALTON/SIRIUS restart file stored under label BASINFO"""
import numpy
from util.unformatted import FortranBinary

class BasInfo():
    """Simple class for BASINFO data"""
    label = "BASINFO"
    def __init__(self, name="SIRIUS.RST"):
        self.name = name
        sirrst = FortranBinary(name)
        sirrst.find(BasInfo.label)
        sirrst.next()
        self.nsym, = sirrst.readbuf(1,'i')
        self.nbas = numpy.array(sirrst.readbuf(8,'i'))
        self.norb = numpy.array(sirrst.readbuf(8,'i'))
        self.nrhf = numpy.array(sirrst.readbuf(8,'i'))
        self.ioprhf, = sirrst.readbuf(1,'i')
        sirrst.close()

    def __repr__(self):
        """Print method for BasInfo objects"""
        pr = lambda v:  len(v)*"%3d" % tuple(v)
        retstr = ""
        retstr += "NSYM   : %3d\n" % self.nsym
        retstr += "NBAS   : %s\n" % pr(self.nbas)
        retstr += "NORB   : %s\n" % pr(self.norb)
        retstr += "NRHF   : %s\n" % pr(self.nrhf)
        retstr += "IOPRHF : %3d\n" % self.ioprhf
        return retstr

if __name__ == "__main__":
    import sys
    try:
        basinfo = BasInfo(sys.argv[1])
    except IndexError:
        print "Usage: %s [path]/SIRIFC" % sys.argv[0]
        sys.exit(1)

    print basinfo
