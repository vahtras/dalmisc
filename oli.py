"""Orbital linear transformation module E2*N, S2*N"""
import os
from daltools import sirifc, one, dens, rspvec
from dalmisc import two

def e2n(N, tmpdir='/tmp', hfx=1, Sg=1, Sv=1):
    """
    E[2]*N linear transformation:

     [k,F] = [k(pq)E(Sv)(pq), F(rs)E(rs)]
           = k(pq)F(rs) (E(ps)d(rq) - E(rq)d(ps))
           = k(pq)F(ps)E(ps) - k(pq)F(rp)E(rq)
           = [k,F](pq)E(pq)
           = kF
    <[q,kF]> = <[E_{pq}, kF(rs)E(rs)]>
             = kF(rs) <E(ps)d(rq) - E(rq)d(ps)>
             = kF(qs)D(ps) - kF(rp)D(rq)
             = [D, kF.T](pq)
    kD = <[k, E(pq)]>
       = <[k(rs) E(rs), E(pq)]>
       = k(rs) (E(rq)d(ps) - E(ps)d(rq))
       = k(rp)D(rq) - k(qs)D(ps)
       = [k.T, D](p,q)
    Fk = F[kD]
    """

    SIRIFC = os.path.join(tmpdir, 'SIRIFC')
    AOONEINT = os.path.join(tmpdir, 'AOONEINT')
    AOTWOINT = os.path.join(tmpdir, 'AOTWOINT')
    LUINDF = os.path.join(tmpdir, 'LUINDF')

    ifc = sirifc.sirifc(SIRIFC)
    cmo = ifc.cmo.unblock()

    h = one.read('ONEHAMIL', filename=AOONEINT).unblock().unpack()
    S = one.read('OVERLAP',  filename=AOONEINT).unblock().unpack()


    dct, dot = dens.ifc(SIRIFC)
    da = dot + dct/2
    db = dct/2

    kN = rspvec.tomat(N, ifc, tmpdir=tmpdir).T #transpose  = (q, q+) to (q+/q)
    kn = cmo*kN*cmo.T

    dak = (kn.T*S*da - da*S*kn.T)
    dbk = (kn.T*S*db - db*S*kn.T)
    #print "dak",dak,"dbk",dbk


    fa, fb = two.fockab((da, db),  filename=AOTWOINT, hfx=hfx)
    fa += h; fb += h
    #print "electronic energy", .5*(((h+fa)&da) + ((h+fb)&db))
    fka, fkb = two.fockab((dak, dbk), filename=AOTWOINT, hfx=hfx)

    kfa = S*kn*fa - fa*kn*S
    kfb = S*kn*fb - fb*kn*S

    fat = fka + kfa
    fbt = fkb + kfb

    gao = S*(da*fat.T + (Sg*Sv)*db*fbt.T) - (fat.T*da + (Sg*Sv)*fbt.T*db)*S
    #print "gao",gao
    gm = cmo.T*gao*cmo
    #print "gm",gm

    # sign convention <[q,[k,F]]> = -E[2]*N
    gv = - rspvec.tovec(gm, ifc)

    return gv

def s2n(N, tmpdir='/tmp', Sg=1, Sv=1):
    """
    S[2]*N linear transformation:
    <[q,k]> = <[E(Sg)_{ij}, k_{kl}, E(Sv)_{kl}]> = 
            = k_{kl} [E_{il}d(kj) - E(kj)d(il)]
            = k(jl)E(SgSv)(il) - k(ki)E(SgSv)(kj)
            = D(SgSv)k.T(ij) - k.TD(SgSv)(ij)
            = [D(SgSv), k.T](ij)
    """


    SIRIFC = os.path.join(tmpdir, 'SIRIFC')
    AOONEINT = os.path.join(tmpdir, 'AOONEINT')
    AOTWOINT = os.path.join(tmpdir, 'AOTWOINT')
    LUINDF = os.path.join(tmpdir, 'LUINDF')

    ifc = sirifc.sirifc(SIRIFC)
    cmo = ifc.cmo.unblock()

    S = one.read('OVERLAP',  filename=AOONEINT).unblock().unpack()


    dct, dot = dens.ifc(SIRIFC)
    da = dot + dct/2
    db = dct/2

    kN = rspvec.tomat(N, ifc, tmpdir=tmpdir).T
    kn = cmo*kN*cmo.T

    dak = (kn.T*S*da - da*S*kn.T)
    dbk = (kn.T*S*db - db*S*kn.T)

    gv = -rspvec.tovec(cmo.T*S*(dak+(Sg*Sv)*dbk)*S*cmo, ifc)

    return gv



if __name__ == "__main__":
    pass
