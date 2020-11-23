#!/usr/bin/env python
from fractions import Fraction
import math
import os
import re
from collections import deque

import numpy


from daltools import one, dens
from two.core import fockab
from util.full import Matrix

from dalmisc import rohf


class SCFIterator():
    def __init__(self):
        self.energies = []
        self.gradient_norms = []


class DaltonOutputIterator(SCFIterator):

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.output = kwargs.get('dalton_output', 'DALTON.OUT')
        self.f = None

    def __iter__(self):
        regex = r'@\s+\d+\s+(-?[\d.]+)\s+([\d.D+-]+)'
        with open(self.output) as f:
            lines = f.read()
            self.lines_iterator = iter(re.findall(regex, lines))
        return self

    def __next__(self):
        e, n = next(self.lines_iterator)
        e = float(e)
        n = float(n.replace('D', 'e'))
        self.energies.append(e)
        self.gradient_norms.append(n)
        return e, n


class RoothanIterator(SCFIterator):

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.it = 0
        self.max_iterations = kwargs.get('max_iterations', 20)
        self.threshold = kwargs.get('threshold', 1e-5)
        self.C = kwargs.get('C0', None)
        if isinstance(self.C, numpy.ndarray):
            self.C = (self.C, self.C)
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
        self.rohf_factors = kwargs.get('factors', (1.0, 1.0))

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
        if self.C is None:
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
            e = self.energy()
            gn = self.gn()
            self.update_mo()
            self.energies.append(e)
            self.gradient_norms.append(gn)
            return (e, gn)
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

        gn = 2*((ga + gb) & (S.I@(ga + gb)@S.I))

        return math.sqrt(gn)

    def update_mo(self):

        F = self.Feff()
        Ca = dens.cmo(F, self.S)
        Cb = Ca
        self.C = Ca, Cb

    def Feff(self):

        F = self.h1 + (self.Fa + self.Fb)/2
        S = self.S

        Da = self.Da
        Db = self.Db
        Fa = self.h1 + self.Fa
        Fb = self.h1 + self.Fb

        ga = S @ Da @ Fa - Fa @ Da @ S
        gb = S @ Db @ Fb - Fb @ Db @ S
        g = ga + gb

        inactive = Db
        active = Da - Db
        virtual = S.I - Da

        # factor_pq = (1.0, 1.0)
        projectors = [(inactive, active), (active, virtual)]
        factors = self.rohf_factors

        if self.ms != 0 and self.Da is not None:
            V = sum(
                S@P@(f*g - F)@Q@S - S@Q@(f*g + F)@P@S
                for f, (P, Q) in zip(factors, projectors)
                # for P, Q in ((Db, Da-Db), (Da-Db, S.I-Da))
            )
            F += V
            # F = S@rohf.Feff(Da@S, Db@S, S.I@Fa, S.I@Fb)
        return F

    def Feff2(self):
        Fa = self.S.I@(self.h1 + self.Fa)
        Fb = self.S.I@(self.h1 + self.Fb)
        Da = self.Da@self.S
        Db = self.Db@self.S
        return self.S@rohf.Feff(Da, Db, Fa, Fb)

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


class DiisIterator(RoothanIterator):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.evecs = deque()
        self.vecs = deque()
        self.max_vecs = kwargs.get('max_vecs', 2)

    def __next__(self):
        """
        Updates for in a SCF iteration
        """
        if not self.converged() and self.it < self.max_iterations:
            self.it += 1
            self.set_densities()
            self.set_focks()
            e = self.energy()
            gn = self.gn()

            self.vecs.append(self.Feff())
            self.evecs.append(self.ga + self.gb)
            if len(self.vecs) > self.max_vecs:
                self.vecs.popleft()
                self.evecs.popleft()

            self.update_mo()
            self.energies.append(e)
            self.gradient_norms.append(gn)
            return (e, gn)
        else:
            raise StopIteration

    def B(self):
        dim = len(self.evecs) + 1
        Bmat = numpy.ones((dim, dim))
        for i, vi in enumerate(self.evecs):
            for j, vj in enumerate(self.evecs):
                Bmat[i, j] = 4*(vi & (self.S.I@vj@self.S.I))

        Bmat[-1, -1] = 0
        return Bmat

    def c(self):
        rhs = numpy.zeros(len(self.evecs) + 1)
        rhs[-1] = 1.0
        return numpy.linalg.solve(self.B(), rhs)[:-1]

    def Fopt(self):
        return sum(
            c*e
            for c, e in zip(self.c(), self.vecs)
        )

    def update_mo(self):

        F = self.Fopt()
        Ca = dens.cmo(F, self.S)
        Cb = Ca
        self.C = Ca, Cb



if __name__ == "__main__":

    if 1:
        print("HeH\n---\n")
        kwargs = dict(
            electrons=3,
            tmpdir='tests/test_heh.d',
            threshold=1e-5,
            max_iterations=20,
            ms=1/2,
            )

        roo = RoothanIterator(**kwargs)
        uroo = URoothanIterator(**kwargs)
        diis = DiisIterator(**kwargs)

        for i, (e, gn) in enumerate(roo, start=1):
            print(f'{i:2d}: {e:14.10f} {gn:.5e}')
        print()
        for i, (e, gn) in enumerate(diis, start=1):
            print(f'{i:2d}: {e:14.10f} {gn:.5e}')
    if 0:
        print("H2O\n---\n")
        kwargs = dict(
            electrons=9,
            tmpdir='tests/test_h2o.d',
            threshold=1e-5,
            max_iterations=20,
            ms=1/2,
            )

        roo = RoothanIterator(**kwargs)
        uroo = URoothanIterator(**kwargs)

        for i, (e, gn) in enumerate(roo, start=1):
            print(f'{i:2d}: {e:14.10f} {gn:.5e}')
