#!/usr/bin/env python
import os
import math
import numpy
from fractions import Fraction
from daltools import one, dens
from two.core import fockab

class SCFIterator():
    pass

class RoothanIterator(SCFIterator):

    def __init__(self, *args, **kwargs):

        self.it = 0
        self.max_iterations = kwargs.get('max_iterations', 10)
        self.threshold = kwargs.get('threshold', 1e-3)
        self.C = kwargs.get('C0', None)
        self.tmpdir = kwargs.get('tmpdir', '/tmp')
        self.nel = kwargs.get('electrons', 0)
        self.ms = Fraction(kwargs.get('ms', 0))
        self.Z = None
        self.h1 = None
        self.S = None
        self.Fa = None
        self.Fb = None
        self.Da = None
        self.Db = None

    def __iter__(self):
        """
        Initial setup for SCF iterations
        """
        self.na = (self.nel + 2*self.ms)//2
        self.nb = (self.nel - 2*self.ms)//2

        AOONEINT = os.path.join(self.tmpdir, 'AOONEINT')
        self.Z = one.readhead(AOONEINT)['potnuc']
        self.h1 = one.read(label='ONEHAMIL', filename=AOONEINT).unpack().unblock()
        self.S = one.read(label='OVERLAP', filename=AOONEINT).unpack().unblock()
        self.Fa = self.h1*0
        self.Fb = self.Fa
        return self

    def __next__(self):
        if not self.converged() and self.it < self.max_iterations:
            self.it += 1
            self.update_mo()
            self.set_densities()
            self.set_focks()
            return self.energy(), self.gn()
        else:
            raise StopIteration
            
    def update_mo(self):
        F = self.h1 + (self.Fa + self.Fb)/2
        if self.ms != 0 and self.Da is not None:
            Fs = (self.Fa - self.Fb)/2
            D = self.Da + self.Db
            Ds = self.Da - self.Db
            ID = numpy.eye(D.shape[0]) - D
            F += (Ds*Fs*ID + ID*Fs*Ds)/2
        Ca = dens.cmo(F, self.S)[:, :self.na]
        Cb = Ca[:, :self.nb]
        self.C = Ca, Cb

    def converged(self):
        if self.it == 0:
            return False
        else:
            return self.gn() < self.threshold

    def set_densities(self):
        Ca, Cb = self.C
        self.Da = Ca*Ca.T
        self.Db = Cb*Cb.T

    def set_focks(self):
        AOTWOINT = os.path.join(self.tmpdir, 'AOTWOINT')
        self.Fa, self.Fb = fockab((self.Da, self.Db), filename=AOTWOINT)


    def energy(self):
        e1 = self.h1&(self.Da + self.Db)
        e2 = 0.5*((self.Da&self.Fa) + (self.Db&self.Fb))
        return e1 + e2 + self.Z

    def gn(self):
        Fa = self.h1 + self.Fa
        Fb = self.h1 + self.Fb
        ga = self.Da*Fa - self.S.I*Fa*self.Da*self.S
        gb = self.Db*Fb - self.S.I*Fb*self.Db*self.S
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
    
