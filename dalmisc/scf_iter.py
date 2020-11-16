#!/usr/bin/env python
import os
import math
from itertools import combinations
from fractions import Fraction

import numpy

from daltools import one, dens
from two.core import fockab


class SCFIterator():
    pass


class RoothanIterator(SCFIterator):

    def __init__(self, *args, **kwargs):

        self.it = 0
        self.max_iterations = kwargs.get('max_iterations', 20)
        self.threshold = kwargs.get('threshold', 1e-5)
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
        self.ga = None
        self.gb = None

    def __iter__(self):
        """
        Initial setup for SCF iterations
        """
        self.na = (self.nel + 2*self.ms)//2
        self.nb = (self.nel - 2*self.ms)//2

        AOONEINT = os.path.join(self.tmpdir, 'AOONEINT')
        self.Z = one.readhead(AOONEINT)['potnuc']
        self.h1 = one.read(
            label='ONEHAMIL', filename=AOONEINT
        ).unpack().unblock()
        self.S = one.read(
            label='OVERLAP', filename=AOONEINT
        ).unpack().unblock()
        self.Ca = dens.cmo(self.h1, self.S)
        self.Cb = self.Ca
        self.C = (self.Ca, self.Cb)
        return self

    def __next__(self):
        """
        Updates for in a SCF iteration
        """
        if not self.converged() and self.it < self.max_iterations:
            self.it += 1
            self.set_densities()
            self.set_focks()
            self.update_mo()
            return self.energy(), self.gn()
        else:
            raise StopIteration
    
    def set_densities(self):
        Ca, Cb = self.C
        self.Da = dens.C1D(Ca, self.na)
        self.Db = dens.C1D(Cb, self.nb)

    def set_focks(self):
        AOTWOINT = os.path.join(self.tmpdir, 'AOTWOINT')
        (self.Fa, self.Fb), = fockab((self.Da, self.Db), filename=AOTWOINT)

    def energy(self):
        e1 = self.h1 & (self.Da + self.Db)
        e2 = 0.5*((self.Da & self.Fa) + (self.Db & self.Fb))
        return e1 + e2 + self.Z

    def gn(self):
        Fa = self.h1 + self.Fa
        Fb = self.h1 + self.Fb
        Da = self.Da
        Db = self.Db
        S = self.S

        self.ga = ga = S@Da@Fa - Fa@Da@S
        self.gb = gb = S@Db@Fb - Fb@Db@S

        gn = ((ga + gb) & (S.I@(ga + gb)@S.I))

        return math.sqrt(gn)

    def update_mo(self):
        F = self.h1 + (self.Fa + self.Fb)/2
        S = self.S

        Da = self.Da
        Db = self.Db
        Fa = self.h1 + self.Fa
        Fb = self.h1 + self.Fb

        ga = S @ Da @ Fa - Fa @ Da @ S
        gb = S @ Db @ Fb - Fb @ Db @ S
        g = ga + gb

        if self.ms != 0 and self.Da is not None:
            V = sum(
                S@P@(g - F)@Q@S
                for P, Q in combinations((Db, Da - Db, S.I - Da), 2)
            )
            F += V + V.T

        Ca = dens.cmo(F, self.S)
        Cb = Ca
        self.C = Ca, Cb

    def converged(self):
        if self.it == 0:
            return False
        else:
            return self.gn() < self.threshold



class URoothanIterator(RoothanIterator):

    def update_mo(self):
        Ca = dens.cmo(self.h1 + self.Fa, self.S)[:, :self.na]
        Cb = dens.cmo(self.h1 + self.Fb, self.S)[:, :self.nb]
        self.C = Ca, Cb

    def gn(self):
        Fa = self.h1 + self.Fa
        Fb = self.h1 + self.Fb
        ga = self.Da@Fa - self.S.I@Fa@self.Da@self.S
        gb = self.Db@Fb - self.S.I@Fb@self.Db@self.S
        Da = self.Da
        Db = self.Db
        Δa = self.S.I - Da
        Δb = self.S.I - Db
        gn = (ga@Δa@ga.T@Da + gb@Δb@gb.T@Db).tr()
        return math.sqrt(gn)


if __name__ == "__main__":

    kwargs = dict(
        electrons=3,
        tmpdir='tests/test_heh.d',
        threshold=1e-5,
        max_iterations=20,
        ms=1/2,
        )

    roo = RoothanIterator(**kwargs)
    uroo = URoothanIterator(**kwargs)

    for i, (e, gn) in enumerate(roo, start=1):
        print(f'{i:2d}: {e:14.10f} {gn:.5e}')
