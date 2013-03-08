#!/usr/bin/env python
"""Module for getting stuff from Dalton two-electron integral file AOTWOINT"""

import numpy as np
from util.unformatted import FortranBinary as FB
from util.full import matrix

def info(filename="AOTWOINT"):
    """Extract info block BASINFO on AOTWOINT"""
    _aotwoint = FB(filename)
    _aotwoint.find("BASINFO")
    rec = next(_aotwoint)
    fileinfo = rec.read(12,'i')
    retinfo = {
        "nsym":fileinfo[0],
        "nbas":fileinfo[1:9],
        "lbuf":fileinfo[9],
        "nibuf":fileinfo[10],
        "nbits":fileinfo[11]
        }
    return retinfo


def list_buffers(filename="AOTWOINT"):
    """ Return integral buffers in AOTWOINT"""
    _aotwoint = FB(filename)
    _aotwoint.find("BASTWOEL")

    for rec in _aotwoint:
        lbuf = (_aotwoint.reclen-4)/12

        buf = np.array(_aotwoint.readbuf(lbuf,'d'))
        ibuf = np.array(_aotwoint.readbuf(4*lbuf,'b')).reshape(lbuf, 4)
        length = _aotwoint.readbuf(1,'i')[0]

        if length < 0: raise StopIteration
        yield buf, ibuf, length

def list_integrals(filename="AOTWOINT"):
    """ List two-electron spin-orbit integrals in file """

    for buf, ibuf, length in list_buffers(filename):
        for g, ig in zip(buf[:length], ibuf[:length]):
            yield ig, g


def yolist_integrals(filename="AOTWOINT"):
    """ List two-electron spin-orbit integrals in file """
    _aotwoint = FB(filename)
    _aotwoint.find("BASTWOEL")

    for rec in _aotwoint:
        lbuf = (_aotwoint.reclen-4)/12

        buf = np.array(_aotwoint.readbuf(lbuf,'d'))
        ibuf = np.array(_aotwoint.readbuf(4*lbuf,'b')).reshape(lbuf, 4)
        length = _aotwoint.readbuf(1,'i')[0]
        #
        # Negative length marks end of file
        #
        if length < 0: break
            #
        for g, ig in zip(buf[:length], ibuf[:length]):
            yield ig, g

    _aotwoint.close()


def fockab(Dab, filename="AOTWOINT", hfc=1, hfx=1, f2py=True):
    """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

    Da, Db = Dab
    D = Da + Db

    if f2py is True:
        try:
            import sirfck
        except(ImportError):
            f2py = False
            print "Warning: non-existent sirfck.so wanted - reverting to python"


    J = matrix(D.shape)
    Ka = matrix(D.shape)
    Kb = matrix(D.shape)

    if f2py:
        for buf, ibuf, length in list_buffers(filename):
            J, Ka, Kb = sirfck.fckab(
                J, Ka, Kb, Da, Db, buf[:length], ibuf.T[:, :length]
                )
    else:
        for ig, g in list_integrals(filename):
            p, q, r, s = ig
            s, r, q, p = (p-1, q-1, r-1, s-1)
            if ( p == q ): g *= .5
            if ( r == s ): g *= .5
            if ( p == r and q == s ): g *= .5

            Jadd = g * (D[r, s] + D[s, r])
            J[p, q] += Jadd
            J[q, p] += Jadd
            Jadd = g * (D[p, q] + D[q, p])
            J[r, s] += Jadd
            J[s, r] += Jadd
            Ka[p, s] += g*Da[r, q]
            Ka[p, r] += g*Da[s, q]
            Ka[q, s] += g*Da[r, p]
            Ka[q, r] += g*Da[s, p]
            Ka[r, q] += g*Da[p, s]
            Ka[s, q] += g*Da[p, r]
            Ka[r, p] += g*Da[q, s]
            Ka[s, p] += g*Da[q, r]
            Kb[p, s] += g*Db[r, q]
            Kb[p, r] += g*Db[s, q]
            Kb[q, s] += g*Db[r, p]
            Kb[q, r] += g*Db[s, p]
            Kb[r, q] += g*Db[p, s]
            Kb[s, q] += g*Db[p, r]
            Kb[r, p] += g*Db[q, s]
            Kb[s, p] += g*Db[q, r]

    Fa = hfc*J-hfx*Ka
    Fb = hfc*J-hfx*Kb
    return (Fa, Fb)

def fock(D, filename="AOTWOINT", hfc=1, hfx=1, f2py=True):
    """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

    if f2py is True:
        try:
            import sirfck
        except(ImportError):
            f2py = False
            print "Warning: non-existent sirfck.so wanted - reverting to python"

    J = matrix(D.shape)
    K = matrix(D.shape)

    if f2py:
        for buf, ibuf, length in list_buffers(filename):
            J, K = sirfck.fck(
                J, K,  D, D, buf[:length], ibuf.T[:, :length]
                )
    else:
        for ig, g in list_integrals(filename):
            p, q, r, s = ig
            s, r, q, p = (p-1, q-1, r-1, s-1)
            if ( p == q ): g *= .5
            if ( r == s ): g *= .5
            if ( p == r and q == s ): g *= .5

            Jadd = g * (D[r, s] + D[s, r])
            J[p, q] += Jadd
            J[q, p] += Jadd
            Jadd = g * (D[p, q] + D[q, p])
            J[r, s] += Jadd
            J[s, r] += Jadd
            K[p, s] += g*D[r, q]
            K[p, r] += g*D[s, q]
            K[q, s] += g*D[r, p]
            K[q, r] += g*D[s, p]
            K[r, q] += g*D[p, s]
            K[s, q] += g*D[p, r]
            K[r, p] += g*D[q, s]
            K[s, p] += g*D[q, r]

    return hfc*J - 0.5*hfx*K

if __name__ == "__main__":
    D = matrix((6, 6))
    D[0, 0] = 1.0
    Fa, Fb = fockab((D, D), "test/test_fockab.d/AOTWOINT", f2py=False)
    print Fa
    Fa, Fb = fockab((D, D), "test/test_fockab.d/AOTWOINT", f2py=True)
    print Fa

if __name__ == "__main__":
    import sys, os, optparse 
    parser = optparse.OptionParser()
    parser.add_option(
       '-d','--directory', dest='dir', default='/tmp', 
       help='Directory containing Dalton job files')
    parser.add_option(
       '-l','--list', dest='list', default=False, action='store_true', 
       help='List integrals on file')

    opt, arg = parser.parse_args(sys.argv[1:])

    #
    # Get ao basis dimension from one-integral file
    #


    aotwoint = os.path.sep.join([opt.dir, "AOTWOINT"])

    if opt.list:
        print "List integrals"
        for ig, g in list_integrals(aotwoint):
            print ig, g


    sys.exit(0)
