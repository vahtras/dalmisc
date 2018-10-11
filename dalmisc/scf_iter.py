#!/usr/bin/env python
import os
import math
from daltools import one, dens
from two.core import fockab

class SCFIterator():
    pass

class RoothanIterator(SCFIterator):

    def __init__(self, *args, **kwargs):

        self.max_iterations = kwargs.get('max_iterations', 10)
        self.threshold = kwargs.get('threshold', 1e-3)
        self.C = kwargs.get('C0', None)
        self.tmpdir = kwargs.get('tmpdir', '/tmp')
        self.Z = None
        self.h1 = None
        self.S = None
        self.nel = kwargs.get('electrons', 0)
        self.ms = kwargs.get('ms', 0)
        #self.unrest = kwargs.get('unrest', False)

        #assert os.path.isdir(self.tmpdir)
        #AOONEINT = os.path.join(self.tmpdir, 'AOONEINT')
        #assert os.path.isfile(AOONEINT)

        #self.h = one.read(label='ONEHAMIL', filename=AOONEINT).unpack().unblock()
        #self.S = one.read(label='OVERLAP', filename=AOONEINT).unpack().unblock()
        #self.Z = one.readhead(AOONEINT)['potnuc']
#
        # default initialization H1DIAG
        #if self.C is None:
        #    Ca = dens.cmo(self.h, self.S)[:na]
        #    if self.unrest:
        #        Cb = dens.cmo(self.h, self.S)[:nb]
        #    else:
        #        Cb = Ca[:nb]
        #    self.C = (Ca, Cb)

    def __iter__(self):
        """This is a self-iterator"""
        self.it = 0
#
# n = na + nb            na = (n + 2ms)/2
# ms = (na - nb)/2       nb = (n - 2ms)/2
#
        self.na = (self.nel + 2*self.ms)//2
        self.nb = (self.nel - 2*self.ms)//2

        AOONEINT = os.path.join(self.tmpdir, 'AOONEINT')
        self.Z = one.readhead(AOONEINT)['potnuc']
        self.h1 = one.read(label='ONEHAMIL', filename=AOONEINT).unpack().unblock()
        self.S = one.read(label='OVERLAP', filename=AOONEINT).unpack().unblock()
#
        #default initialization H1DIAG
        if self.C is None:
            Ca = dens.cmo(self.h1, self.S)[:, :self.na]
            Cb = Ca[:, :self.nb]
            self.C = (Ca, Cb)
        return self

    def __next__(self):
        if self.gn() > self.threshold and  self.it < self.max_iterations:
            self.it += 1
            Fa, Fb = self.focks()
            F = self.h1 + (Fa + Fb)/2
            Ca = dens.cmo(F, self.S)[:, :self.na]
            Cb = Ca[:, :self.nb]
            self.C = Ca, Cb
            e = self.energy()
            gn = self.gn()
            return e, gn
        else:
            raise StopIteration

    def densities(self):
        Ca, Cb = self.C
        Da = Ca*Ca.T
        Db = Ca*Ca.T
        return Da, Db

    def focks(self):
        Da, Db = self.densities()
        AOTWOINT = os.path.join(self.tmpdir, 'AOTWOINT')
        Fa, Fb = fockab((Da, Db), filename=AOTWOINT)
        return Fa, Fb


    def energy(self):
        Da, Db = self.densities()
        e1 = self.h1&(Da + Db)
        Fa, Fb = self.focks()
        e2 = 0.5*((Da&Fa) + (Db&Fb))
        return e1 + e2 + self.Z

    def gn(self):
        AOTWOINT = os.path.join(self.tmpdir, 'AOTWOINT')
        Ca, Cb = self.C
        Da = Ca*Ca.T
        Db = Cb*Cb.T
        Fa, Fb = fockab((Da, Db), filename=AOTWOINT)
        Fa += self.h1
        Fb += self.h1
        ga = Da*Fa - self.S.I*Fa*Da*self.S
        gb = Db*Fb - self.S.I*Fb*Db*self.S
        gn = ((ga + gb)**2).tr()
        return  gn
        



if __name__ == "__main__":
    import sys
    import optparse 
    op = optparse.OptionParser()
    op.add_option(
        '-t', '--tmpdir', dest='tmpdir', default='/tmp',
        help='Scratch directory'
        )
    o, a = op.parse_args(sys.argv[1:])
    na = int(a[0])
    nb = int(a[1])

    for E, gn in Roothan(5, 5):
        print(E, gn)
    
