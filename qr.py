#!/usr/bin/env python
"""
Module for evaluating quadratic response functions at SCF level from LR vectors
"""
import os, sys
from daltools import dens, prop, lr, sirifc, rspvec, one
from util.full import matrix
#local
import two


def E3(pB, pC, ifc, tmpdir='/tmp'):
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


    AOONEINT = os.path.join(tmpdir, "AOONEINT")
    h = one.read(label='ONEHAMIL', filename=AOONEINT).unpack().unblock()
    S = one.read(label='OVERLAP', filename=AOONEINT).unblock().unpack()

    AOTWOINT = os.path.join(tmpdir, "AOTWOINT")


    kB = pB["kappa"]
    kC = pC["kappa"]

    sB = pB.get("spin", 1)
    sC = pC.get("spin", 1)
   
    cmo = ifc.cmo.unblock()
  
    kb = cmo*kB*cmo.T
    kc = cmo*kC*cmo.T
 

    _kb = S*kb
    kb_ = kb*S
    _kc = S*kc
    kc_ = kc*S
    #
    # Fock matrices
    #
    da, db = dens.Dab(ifc_=ifc)

    daB = da*S*kb - kb*S*da 
    daC = da*S*kc - kc*S*da 
    daBC = (daB*S*kc - kc*S*daB + daC*S*kb - kb*S*daC)/2

    dbB = (db*S*kb - kb*S*db)*sB
    dbC = (db*S*kc - kc*S*db)*sC
    dbBC = ((dbB*S*kc - kc*S*dbB)*sC + (dbC*S*kb - kb*S*dbC)*sB)/2

    fa, fb = two.fockab((da, db), filename=AOTWOINT)
    fa += h
    fb += h

    faB, fbB = two.fockab((daB, dbB), filename = AOTWOINT)
    faC, fbC = two.fockab((daC, dbC), filename = AOTWOINT)
    faBC, fbBC = two.fockab((daBC, dbBC), filename = AOTWOINT)

    Bfa = S*kb*fa-fa*kb*S
    Cfa = S*kc*fa-fa*kc*S
    BCfa = (S*kb*Cfa-Cfa*kb*S + S*kc*Bfa-Bfa*kc*S)/2
    BfaC =  S*kb*faC-faC*kb*S + S*kc*faB-faB*kc*S
   
    Bfb = (S*kb*fb-fb*kb*S)*sB
    Cfb = (S*kc*fb-fb*kc*S)*sC
    BCfb = (S*kb*Cfb-Cfb*kb*S*sB + S*kc*Bfb-Bfb*kc*S*sC)/2
    BfbC = S*kb*fbC-fbC*kb*S*sB + S*kc*fbB-fbB*kc*S*sC
 #
 # Add all focks
 #
 
    fa = faBC+BfaC+BCfa
    fb = fbBC+BfbC+BCfb

    G = (cmo.T*((fa*da+fb*db)*S - S*(da*fa+db*fb))*cmo).T

    Gv = rspvec.tovec(G, ifc)

    return Gv


def main(*args, **kwargs):

    labs = args
    tmpdir = kwargs.get("tmpdir", "/tmp")
    ranks = kwargs.get('rank', (0, 0, 0))
    pars = [ (-1)**r for r in ranks]

    global tmp

    ifc = sirifc.sirifc(os.path.join(tmpdir, "SIRIFC"))
    cmo = ifc.cmo.unblock()

    AOONEINT = os.path.join(tmpdir, "AOONEINT")
    AOPROPER = os.path.join(tmpdir, "AOPROPER")
    RSPVEC   = os.path.join(tmpdir, "RSPVEC")

    vecs = [rspvec.read(lab, RSPVEC) for lab in labs]
    kappa = [rspvec.tomat(vec, ifc) for vec in vecs]
    ops = [cmo.T*prop.read(lab, AOPROPER) for lab in labs]

    kA, kB, kC = kappa
    a, b, c = ops
    pA, pB, pC = [{ "lab":lab, "rank":rank, 'kappa':k} for lab, rank, k in zip(labs, ranks, kappa)]

    dc, do = dens.ifc(ifc_=ifc)
    d = dc+do
    
    S = one.read(label = "OVERLAP", filename = AOONEINT).unblock().unpack()
    D = cmo.T*S*d*S*cmo

    E3BC = E3(pB, pC, ifc, tmpdir)
    AE3BC = -NA&E3BC
    B2C = (-(kA^(kC^b))&D)
    C2B = (-(kA^(kB^c))&D)
    A2B = (.5*(kC^(kB^a))&D)
    A2C = (.5*(kB^(kC^a))&D)
    #print "E3BC",E3BC
    val = AE3BC
    print "E3  %14.8f %14.8f" % (AE3BC, val)
    val += B2C
    print "B2C %14.8f %14.8f" % (B2C, val)
    val += C2B
    print "C2B %14.8f %14.8f" % (C2B, val)
    val += A2B
    print "A2B %14.8f %14.8f" % (A2B, val)
    val += A2C
    print "A2C %14.8f %14.8f" % (A2C, val)
    return (val, AE3BC, B2C, C2B, A2B, A2C)

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
    a = cmo.T*prop.read(A, AOPROPER).unpack()*cmo
    kB = rspvec.tomat(NB, ifc, tmpdir = tmp).T
    kC = rspvec.tomat(NC, ifc, tmpdir = tmp).T

    S = one.read(label = "OVERLAP", filename = AOONEINT).unblock().unpack()
    D = cmo.T*S*d*S*cmo
    A2B = (.5*(kC^(kB^a))&D)
    A2C = (.5*(kB^(kC^a))&D)

    return A2B + A2C

if __name__ == "__main__":
    global tmpdir
    try:
        A = sys.argv[1]
        B = sys.argv[2]
        C = sys.argv[3]
    except(IndexError):
        print "Usage: %s A B C" % sys.argv[0]
        sys.exit(1)
    print "QR %14.8f" % main(*sys.argv[1:])[0]
