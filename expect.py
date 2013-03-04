#!/usr/bin/env python
"""Experimental module for expectation value"""
import os
import sys
from daltools import dens, prop

def value(label, tmpdir='/tmp'):
    """Calculate expectation of operator with label"""   
    aoproper = os.path.join(tmpdir, 'AOPROPER')
    sirifc = os.path.join(tmpdir, 'SIRIFC')
    A = prop.read(label, propfile=aoproper).unpack()
    dc, do = dens.ifc(filename=sirifc)
    return A&(dc + do)


if __name__ == "__main__":
    from optparse import OptionParser

    op = OptionParser()
    op.add_option(
        '-t', '--tmpdir', dest='tmpdir', default='/tmp', 
        help='Dalton scratch directory'
        )
    opt, arg = op.parse_args(sys.argv[1:])
    try:
        A, = arg
    except:
        print "Usage:%s [-t <tmpdir>] A" % sys.argv[0]
        raise SystemExit
    
    print value(A, tmpdir = opt.tmpdir)

   

