#!/usr/bin/env python
import numpy
from util import unformatted,blocked,full

class sirrst(unformatted.FortranBinary):
    def __init__(self,name="SIRIUS.RST"):
        unformatted.FortranBinary.__init__(self,name)
        self.getdata()
        self.getcmo()
        self.close()

    def __str__(self):
        retstr=""
        retstr+="\nNSYM : " + str(self.nsym)
        retstr+="\nNBAS : " + str(self.nbas)
        retstr+="\nNBAST: " + str(self.nbast)
        retstr+="\nNORB : " + str(self.norb)
        retstr+="\nNORBT: " + str(self.norbt)
        retstr+="\nNRHF : " + str(self.nrhf)
        retstr+="\nIOPRHF:" + str(self.ioprhf)
        retstr+="\nCMO:   " + str(self.cmo)
        return retstr

    def getdata(self):
        self.find("BASINFO")
        self.readrec()
        self.nsym, = self.readbuf(1,'i')
        self.nbas = numpy.array(self.readbuf(8,'i'))#[:self.nsym]
        self.nbast = self.nbas.sum()
        self.norb = numpy.array(self.readbuf(8,'i'))#[:self.nsym]
        self.norbt = self.norb.sum()
        self.nrhf = numpy.array(self.readbuf(8,'i'))#[:self.nsym]
        self.ncmo = self.nbas*self.norb
        self.ioprhf, = self.readbuf(1,'i')

    def getcmo(self):
        self.find("NEWORB")
        self.ncmot=numpy.dot(self.nbas, self.norb)
        ncmot4=max(self.ncmot,4)
        dbl=self.readrec()
        n=0
        self.cmo=blocked.matrix(self.nbas,self.norb)
        for isym in range(self.nsym):
            cmoi = numpy.array(self.readbuf(self.ncmo[isym],'d')
                   ).reshape((self.nbas[isym], self.norb[isym]), order='F')
            self.cmo.subblock[isym] = cmoi.view(full.matrix)
            n += self.ncmo[isym]
        assert(n == self.ncmot)

if __name__ == "__main__":
    import os, sys
    try:
        rst=sirrst(sys.argv[1])
    except IndexError:
        print "Usage: %s [<path>/]SIRIUS.RST"
        sys.exit(1)
    print rst

