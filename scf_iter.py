#!/usr/bin/env python
import os
import math
from daltools import one, dens
from dalmisc import two

class Roothan(object):

    def __init__(self, *args, **kwargs):

        na, nb = args

        self.C = kwargs.get('C0', None)
        self.threshold = kwargs.get('threshold', 1e-3)
        self.max_iterations = kwargs.get('max_iterations', 10)
        self.tmpdir = kwargs.get('tmpdir', '/tmp')
        self.unrest = kwargs.get('unrest', False)

        assert os.path.isdir(self.tmpdir)
        AOONEINT = os.path.join(self.tmpdir, 'AOONEINT')
        assert os.path.isfile(AOONEINT)

        self.h = one.read(label='ONEHAMIL', filename=AOONEINT).unpack().unblock()
        self.S = one.read(label='OVERLAP', filename=AOONEINT).unpack().unblock()
        self.Z = one.readhead(AOONEINT)['potnuc']

        # default initialization H1DIAG
        if self.C is None:
            Ca = dens.cmo(self.h, self.S)[:na]
            if self.unrest:
                Cb = dens.cmo(self.h, self.S)[:nb]
            else:
                Cb = Ca[:nb]
            self.C = (Ca, Cb)

    def __iter__(self):
        """This is a self-iterator"""
        self.it = 0
        return self

    def next(self):
        if self.gn > self.threshold and  self.it < self.max_iterations:
            self.it += 1
            Ca, Cb = self.C
            Da = Ca*Ca.T
            Db = Cb*Cb.T
            Fa, Fb = two.fockab(Da, Db)
            Fa += self.h
            Fb += self.h
            E = 0.5*((Da&(h+Fa) + (Db&(H+Fb)))) + self.Z
            ga = Da*Fa - self.S.I*Fa*Da*self.S
            gb = Db*Fb - self.S.I*Fb*Db*self.S
            if self.unrest:
                gn = math.sqrt(-(ga*ga + gb*gb).tr())
            else:
                gn = math.sqrt(-(ga+gb)**2)/2
            return E, gn
        raise StopIteration


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
        print E, gn
    
