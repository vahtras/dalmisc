#!/usr/bin/env python
"""
Module for evaluating quadratic response functions at SCF level from LR vectors
"""
import os, sys
from daltools import dens, prop, lr, sirifc, rspvec, one
from util.full import matrix
#local
import two


def E3(pB, pC, ifc, **kwargs):
    """ Emulate the so called E3 contribution to a quadratic response function
        <<A; B, C>> = NA  E3 (NB NC +  NC NB) + A2 (NB NC + NC NB) + NA (B2 NC + C2 NB)

        Emulation of current Dalton implementation in terms of high spin fock matrices
        Closed and open shell matrices
        Dc = inactive 
        Do = -active
        Fc = Fa+Q
        Fo = ? CHECK

        Formulas
        1/2*[qa, [kb, [kc, H]]] + P(b,c)

        [kc, H] = (p~q|rs)H(Sc, 0) + (pq|r~s) H(0, Sc)
        [kb, [kc, H]] 
                = (p~~q|rs)H(SbSc, 0) + (p~q|r~s) H(Sb, Sc)
                + (p~q|r~s)H(Sc, Sb) + (pq|r~~s) H(0, SbSc)

        and for 
        H(S1, S2) generates Fock from (D(S1) g D(S2) - Da g Da - Db g Db
        F(S1, S2) = E(S1) g D(S2) - D(S1) g E(S2) - Ea g Da - Da g Ea - Eb g Db - Db g Eb
                  = Ea [ g D(S2) - D(S1) g - g Da - Da g ]
                  + Eb [ S1 g D(S2) - D(S1) g S2  - g Db - Db g ]
    """


    tmpdir = kwargs.get('tmpdir', '/tmp')
    AOONEINT = os.path.join(tmpdir, "AOONEINT")
    h = one.read(label='ONEHAMIL', filename=AOONEINT).unpack().unblock()
    S = one.read(label='OVERLAP', filename=AOONEINT).unblock().unpack()

    AOTWOINT = os.path.join(tmpdir, "AOTWOINT")
    kwargs['filename'] = AOTWOINT


    cmo = ifc.cmo.unblock()
    kB = cmo@pB["kappa"]@cmo.T
    kC = cmo@pC["kappa"]@cmo.T
    kB_ = kB@S
    _kB = S@kB
    kC_ = kC@S
    _kC = S@kC

    sB = pB.get("spin", 1)
    sC = pC.get("spin", 1)
   
    #
    # Fock matrices
    #
    da, db = dens.Dab(ifc_=ifc)
    (fa, fb), = two.fockab((da, db), **kwargs)
    fa += h
    fb += h
    Bfa, Bfb = [_kB@f - f@kB_ for f in (fa, sB*fb)]
    Cfa, Cfb = [_kC@f - f@kC_ for f in (fa, sC*fb)]
    
    BCfa, BCfb = [_kB@Cf - Cf@kB_ for Cf in (Cfa, sB*Cfb)]
    CBfa, CBfb = [_kC@Bf - Bf@kC_ for Bf in (Bfa, sC*Bfb)]

    daB, dbB = [_kB.T@d - d@kB_.T  for d in (da, sB*db)]
    (faB, fbB), = two.fockab((daB, dbB), **kwargs)
    CfaB, CfbB = [_kC@fB - fB@kC_ for fB in (faB, sC*fbB)]

    daC, dbC = [_kC.T@d - d@kC_.T  for d in (da, sC*db)]
    (faC, fbC), = two.fockab((daC, dbC), **kwargs)
    BfaC, BfbC = [_kB@fC - fC@kB_ for fC in (faC, sB*fbC)]
    
    daBC, dbBC  = (_kB.T@dC - dC@kB_.T for dC in (daC, sB*dbC))
    daCB, dbCB  = (_kC.T@dB - dB@kC_.T for dB in (daB, sC*dbB))
    daBC = 0.5*(daBC + daCB)
    dbBC = 0.5*(dbBC + dbCB)
    (faBC, fbBC), = two.fockab((daBC, dbBC), **kwargs)

#
# Add all focks
#
    fa = faBC + BfaC + CfaB + .5*(BCfa + CBfa)
    fb = fbBC + BfbC + CfbB + .5*(BCfb + CBfb)

    G = cmo.T@(S@(da@fa.T + db@fb.T) - (fa.T@da + fb.T@db)@S)@cmo
    Gv = rspvec.tovec(G, ifc)

    return Gv


def B2C(*args, **kwargs):

    pB, pC, ifc = args
    tmpdir = kwargs.get("tmpdir", ".")

    AOONEINT = os.path.join(tmpdir, "AOONEINT")
    S = one.read(label="OVERLAP", filename=AOONEINT).unblock().unpack()

    cmo = ifc.cmo.unblock()
    mB = pB["matrix"]
    mC = pC["matrix"]
    kB = cmo@pB["kappa"]@cmo.T
    kC = cmo@pC["kappa"]@cmo.T

    kBmC = S@kB@mC - mC@kB@S
    kCmB = S@kC@mB - mB@kC@S
    BC = kBmC + kCmB

    da, db = dens.Dab(ifc_=ifc)
    G = cmo.T@(S@(da@BC.T + db@BC.T) - (BC.T@da + BC.T@db)@S)@cmo
    Gv = rspvec.tovec(G, ifc)
    return Gv


def A2B(*args, **kwargs):
   
    pA, pB, ifc = args
    tmpdir = kwargs.get("tmpdir", ".")

    AOONEINT = os.path.join(tmpdir, "AOONEINT")
    S = one.read(label = "OVERLAP", filename = AOONEINT).unblock().unpack()

    cmo = ifc.cmo.unblock()
    mA = pA["matrix"]
    kB = cmo@pB["kappa"]@cmo.T

    BA = S@kB@mA - mA@kB@S

    da, db = dens.Dab(ifc_=ifc)
    G = cmo.T@(S@(da@BA.T + db@BA.T) - (BA.T@da + BA.T@db)@S)@cmo
    Gv = rspvec.tovec(G, ifc)
    return Gv

    
    

def main(*args, **kwargs):

    labs = args
    tmpdir = kwargs.get("tmpdir", ".")
    ranks = kwargs.get('rank', (0, 0, 0))
    pars = [ (-1)**r for r in ranks]

    ifc = sirifc.sirifc(os.path.join(tmpdir, "SIRIFC"))
    cmo = ifc.cmo.unblock()

    AOONEINT = os.path.join(tmpdir, "AOONEINT")
    AOPROPER = os.path.join(tmpdir, "AOPROPER")
    RSPVEC   = os.path.join(tmpdir, "RSPVEC")

    vecs = rspvec.read(*labs, propfile=RSPVEC)[0]
    kappa = [rspvec.tomat(vec, ifc).T for vec in vecs]
    kappa[0] = kappa[0].T
    a, b, c = [cmo.T@x@cmo for x in prop.read(*labs, filename=AOPROPER)]

    NA = vecs[0]
    kA, kB, kC = kappa
    pA, pB, pC = [{ "lab":lab, "rank":rank, 'kappa':k} for lab, rank, k in zip(labs, ranks, kappa)]

    da, db = dens.Dab(ifc_=ifc)
    d = da + db
    
    S = one.read(label = "OVERLAP", filename = AOONEINT).unblock().unpack()
    D = cmo.T@S@d@S@cmo

    E3BC = E3(pB, pC, ifc, tmpdir=tmpdir)
    AE3BC = -NA&E3BC
    B2C = (-(kA^(kC^b))&D)
    C2B = (-(kA^(kB^c))&D)
    A2B = (.5*(kC^(kB^a))&D)
    A2C = (.5*(kB^(kC^a))&D)
    #print "E3BC",E3BC
    val = AE3BC
    print("E3  %14.8f %14.8f" % (AE3BC, val))
    val += B2C
    print("B2C %14.8f %14.8f" % (B2C, val))
    val += C2B
    print("C2B %14.8f %14.8f" % (C2B, val))
    val += A2B
    print("A2B %14.8f %14.8f" % (A2B, val))
    val += A2C
    print("A2C %14.8f %14.8f" % (A2C, val))
    
    return val

def a2bc(A, B, C):

    AOONEINT = os.path.join(tmp, "AOONEINT")
    AOPROPER = os.path.join(tmp, "AOPROPER")
    RSPVEC   = os.path.join(tmp, "RSPVEC")
    SIRIFC   = os.path.join(tmp, "SIRIFC")

    NB = rspvec.read(B, RSPVEC)
    NC = rspvec.read(C, RSPVEC)

    ifc = sirifc.sirifc(SIRIFC)
    cmo = ifc.cmo.unblock()
    dc,do = dens.ifc(ifc=ifc)
    d = dc+do
    a = cmo.T@prop.read(A, AOPROPER).unpack()@cmo
    kB = rspvec.tomat(NB, ifc, tmpdir = tmp).T
    kC = rspvec.tomat(NC, ifc, tmpdir = tmp).T

    S = one.read(label = "OVERLAP", filename = AOONEINT).unblock().unpack()
    D = cmo.T@S@d@S@cmo
    A2B = (.5*(kC^(kB^a))&D)
    A2C = (.5*(kB^(kC^a))&D)

    return A2B + A2C

if __name__ == "__main__":
    import argparse
    global tmpdir

    parser = argparse.ArgumentParser()
    parser.add_argument('A', help="Property A")
    parser.add_argument('B', help="Perturbation B")
    parser.add_argument('C', help="Perturbation C")
    parser.add_argument('--tmpdir', default='.', help="Dalton tmp directory")

    args = parser.parse_args()
    tmpdir = args.tmpdir

    print("QR %14.8f" % main(args.A, args.B, args.C, tmpdir=args.tmpdir))
