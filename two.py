#!/usr/bin/env python
"""Module for getting stuff from Dalton two-electron integral file AOTWOINT"""

import numpy as np
from util.unformatted import FortranBinary as fb
from fockab import fockab

def info(filename="AOTWOINT"):
    """Extract info block BASINFO on AOTWOINT"""
    _aotwoint = fb(filename)
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



def list_integrals(filename="AOTWOINT"):
    """ List two-electron spin-orbit integrals in file """
    _aotwoint = fb(filename)
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

def fock(*args, **kwargs):
    """Placeholder"""
    return fockab(*args, **kwargs)

if __name__ == "__main__":
    import sys, os, optparse 
    parser = optparse.OptionParser()
    parser.add_option(
       '-d','--directory', dest='dir', default='/tmp', 
       help='Directory containing Dalton job files')
    parser.add_option(
       '-f','--fock', dest='fock', default=False, action='store_true', 
       help='Calculate the [inactive] Fock matrix')
    parser.add_option(
       '-g','--grad', dest='grad', default=0, type=int, 
       help='Calculate the [inactive] Fock matrix')
    parser.add_option(
       '-l','--list', dest='list', default=False, action='store_true', 
       help='List integrals on file')
    parser.add_option(
       '-o','--occ', dest='occ', default=None, 
       help='List integrals on file')
    parser.add_option(
       '-a','--fock-ab', dest='fock_ab', default=False, action='store_true', 
       help='Calculate the alpha,beta Fock matrices')
    parser.add_option(
       '-J','--hfc', dest='hfc', default=1.0, type=float, 
       help='Hartree-Fock coulomb factor')
    parser.add_option(
       '-K','--hfx', dest='hfx', default=1.0, type=float, 
       help='Hartree-Fock exchange factor')
    parser.add_option(
       '-S','--spin-density', dest='S', default=1, type=int, 
       help='Choose input density (0,1)')
    parser.add_option(
       '-O','--optimized', dest='O', default=False, action='store_true', 
       help='Calculate Fock matrix with compiled Fortran module (sirfck.so)')

    opt, arg = parser.parse_args(sys.argv[1:])

    #
    # Get ao basis dimension from one-integral file
    #



    if opt.fock:
        if opt.grad not in  [0, 1]:
            print "Valid opt: 0,1"
            sys.exit(1)
        if opt.grad:
            print "Get triplet Fock matrix"
        else:
            print "Get singlet Fock matrix"

    if opt.occ:
        occ = [int(i) for i in opt.occ.split(",")]
        print occ

    #
    # Check file
    #

    aotwoint = os.path.sep.join([opt.dir, "AOTWOINT"])

    if opt.list:
        print "List integrals"
        for ig, g in list_integrals(aotwoint):
            print ig, g

    if opt.fock or opt.fock_ab:
##
## M from from ifc
##
        SIRIUS_RST = os.path.sep.join([opt.dir, "SIRIUS.RST"])
        from sirrst import sirrst
        rst = sirrst(SIRIUS_RST)
        cmo = rst.cmo.unblock()

        if opt.occ:
            na, nb = occ
        else:
            na = nb = cmo.shape[1]

        cmoa = cmo[:, :na]
        cmob = cmo[:, :nb]

        Da = cmoa*cmoa.T
        Db = cmob*cmob.T

    if opt.fock:

        D = Da + opt.S*Db
        F = fock(D, 
            file=aotwoint,  hfc=opt.hfc, hfx=opt.hfx, f2py=opt.O
            )
        print D, F

    if opt.fock_ab:
        Fa, Fb = fockab(
            (Da, Db), filename=aotwoint, hfc=opt.hfc, hfx=opt.hfx, 
            f2py=opt.O
            )
        print  "Fa", Fa, "Fb", Fb

    sys.exit(0)
