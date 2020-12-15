#!/usr/bin/env python

import math
import os
import pathlib
from itertools import combinations

import numpy as np

from util.full import Matrix
from daltools import one
import two


"""A lot of cleanup to do"""


def grad(S, C, Dc, Do, Fc, Fo):
    g = 2 * C.T @ S @ (Dc @ Fc + Do @ Fo) @ C
    g = g - g.T
    return g


def gradao(S, Dc, Do, Fc, Fo):
    g = 2 * S * (Dc * Fc + Do * Fo)
    g = g - g.T
    return g


def gradvec(S, C, Dc, Do, Fc, Fo, nisht=0, nasht=0, nocct=0):
    g = grad(S, C, Dc, Do, Fc, Fo)
    nisht = int((S & Dc) / 2 + 0.5)
    nasht = int((S & Do) + 0.5)
    norbt, _ = C.shape
    nocct = nisht + nasht
    # print nisht,nasht,nocct,norbt
    dim = nisht * (norbt - nisht) + nasht * (norbt - nocct)
    gv = full.matrix((dim, 2))
    ij = 0
    for i in range(nisht):
        for j in range(nisht, nocct):
            # print ij+1,i+1,j+1
            gv[ij, 0] = -g[i, j]
            gv[ij, 1] = -g[j, i]
            ij += 1
    for i in range(nisht):
        for j in range(nocct, norbt):
            # print ij+1,i+1,j+1
            gv[ij, 0] = -g[i, j]
            gv[ij, 1] = -g[j, i]
            ij += 1
    for i in range(nisht, nocct):
        for j in range(nocct, norbt):
            # print ij+1,i+1,j+1
            gv[ij, 0] = -g[i, j]
            gv[ij, 1] = -g[j, i]
            ij += 1
    return gv


def gradmatrix(p, nisht, nasht, norbt):
    new = full.matrix(norbt, norbt)
    nocct = nisht + nasht
    ij = 0
    for i in range(nisht):
        for j in range(nisht, norbt):
            new[i, j] = -p[ij, 0]
            new[j, i] = p[ij, 0]
            ij += 1
    for i in range(nisht, nocct):
        for j in range(nocct, norbt):
            new[i, j] = -p[ij, 0]
            new[j, i] = p[ij, 0]
            ij += 1
    return new


def gradnorm(S, C, Dc, Do, Fc, Fo, nisht, nasht):
    g = grad(S, C, Dc, Do, Fc, Fo)
    norbt = C.shape[1]
    gsum = 0.0
    for i in range(nisht):
        for j in range(nisht, norbt):
            gsum += g[i, j] ** 2
    for i in range(nisht, nisht + nasht):
        for j in range(nisht + nasht, norbt):
            gsum += g[i, j] ** 2
    return math.sqrt(gsum)


def jensen(S, Dc, Do, Fc, Fo, C=None, h=None):
    Id = full.unit(len(S))
    corr = S * Do * (Fc - Fo) * ((Dc + Do) * S - Id)
    F = Fc + corr + corr.T
    return F


# def jensen_scaled(S,Dc,Do,Fc,Fo):
#   I=full.unit(S.rdim)
#   corr=S*Dc*(Fc-Fo)*((Dc+Do)*S-I)
#   F=Fc
#   for i in range(nisht):
#      for j in range(nisht+nasht,F.rdim):
#         F[i,j]*=2
#         F[j,i]*=2
#   return F+corr+corr.T
# def filatov(S,Dc,Do,Fc,Fo):
#   I=full.unit(S.rdim)
#   f=0.5
#   beta=1./(1-f)
#   PS = Do*S
#   SG = S*Dc/2 - (1/beta)*I + (1/(2*beta))*S*Do
#   F=Fc
#   corr=beta*SG*(Fc-Fo)*PS
#   F+=corr+corr.T
#   return F
# def okazaki(S,Dc,Do,Fc,Fo):
#   Dc=0.5*Dc
#   Fo=0.5*Fo
#   I=full.unit(S.rdim)
#   F=(I-S*Do)*Fc*(I-Do*S) + (I-S*Dc)*Fo*(I-Dc*S)\
#    +S*(Dc*(Fc-Fo)*Do + Do*(Fc-Fo)*Dc)*S
#   return F
#


class Increase(Exception):
    pass


class Stop(Exception):
    pass


class Converged(Exception):
    """To be removed"""

    def __init__(self, value):
        super().__init__()
        self.value = value

    def __str__(self):
        return repr(self.value)


# class BackstepFail(Exception):
#   def __init__(self,value):
#      self.value=value
#   def __str__(self):
#      return repr(self.value)

from math import sqrt
from util import full
from daltools import one, dens
from two import fockab


def energy(Da, Db, h1, Fa, Fb):
    e1 = h1 & (Da + Db)
    e2 = 0.5 * ((Da & Fa) + (Db & Fb))
    return e1 + e2


def uroothan(
    Ca, Cb, na, nb, hfx=1, iters=10, threshold=1e-6, unrest=False, wrkdir="/tmp"
):
    """
    Open-shell Roothan iterations, restricted or unrestricted
    """
    if unrest:
        print("Unrestricted HF Na=%d Nb=%d\n" % (na, nb))
    else:
        print("Restricted RHF Na=%d Nb=%d\n" % (na, nb))
    E0 = 0.0
    aooneint = os.path.join(wrkdir, "AOONEINT")
    aotwoint = os.path.join(wrkdir, "AOTWOINT")
    h = one.read("ONEHAMIL", aooneint).unpack().unblock()
    S = one.read("OVERLAP", aooneint).unpack().unblock()
    potnuc = one.readhead(aooneint)["potnuc"]
    print(potnuc)
    iterinf = []
    try:
        for i in range(iters):
            Da = dens.C1D(Ca, na)
            Db = dens.C1D(Cb, nb)
            (Fa, Fb), = fockab((Da, Db), hfx=hfx, filename=aotwoint)
            Fa += h
            Fb += h
            E = 0.5 * ((Da & (h + Fa)) + (Db & (h + Fb))) + potnuc
            ga = S @ Da @ Fa - Fa @ Da @ S
            gb = S @ Db @ Fb - Fb @ Db @ S
            g = ga + gb
            if unrest:
                g2 = -(ga @ ga + gb @ gb)
            else:
                g2 = (ga + gb) @ (S.I@(ga + gb)@S.I).T
            gn = sqrt(g2.tr())
            gn = sqrt(2*(ga + gb) & (S.I@(ga + gb)@S.I))
            iterinf.append((E, gn))
            print("%2d:E=%16.10f %16.5e %16.2e" % (i + 1, E, gn, E - E0))
            if gn < threshold:
                raise Converged(gn)
            if unrest:
                Ca = dens.cmo(Fa, S)
                Cb = dens.cmo(Fb, S)
            else:
                D = Da + Db
                Ds = Da - Db
                Fs = Fa - Fb
                ID = S.I - D
                F0 = (Fa + Fb)/2
                V = sum(
                    S@P@(g - F0)@Q@S
                    for P, Q in combinations((Db, Da - Db, S.I - Da), 2)
                )
                V += V.T
                V = jensen(S, 2*Db, Da - Db, F0, Fa)
                # F = ((Fa + Fb) + Ds @ Fs @ ID + ID @ Fs @ Ds) / 2
                F = F0 + V
                F = S@Feff(Da@S, Db@S, S.I@Fa, S.I@Fb)
                Ca = dens.cmo(F, S)
                Cb = Ca
    except Converged:
        print("-Converged-")
    if unrest:
        return (Ca, Cb)
    else:
        return Ca


def mkB(vecs):
    B = full.Matrix((len(vecs) + 1, len(vecs) + 1))
    for i in range(len(vecs)):
        for j in range(len(vecs)):
            B[i, j] = vecs[i] & vecs[j]
        B[i, len(vecs)] = 1
        B[len(vecs), i] = 1
    return B


def mkB2(vecs):
    B = full.matrix(len(vecs) + 1, len(vecs) + 1)
    for i in range(len(vecs)):
        for j in range(len(vecs)):
            B[i, j] = -(vecs[i] * vecs[j]).tr()
        B[i, len(vecs)] = -1
        B[len(vecs), i] = -1
    return B


def mkB3(avecs, bvecs, unrest):
    dim = len(avecs)
    B = full.matrix((dim + 1, dim + 1))
    for i in range(dim):
        for j in range(dim):
            if unrest:
                B[i, j] = -(avecs[i] @ avecs[j]).tr() - (bvecs[i] @ bvecs[j]).tr()
            else:
                vecsi = avecs[i] + bvecs[i]
                vecsj = avecs[j] + bvecs[j]
                B[i, j] = (vecsi @ vecsj.T).tr()
        B[i, dim] = -1
        B[dim, i] = -1
    return B


#
# def Eg(C,na,nb,hfx=1):
#   import one
#   potnuc=one.info()[0]["potnuc"]
#   S=one.read('OVERLAP','AOONEINT').unpack().unblock()
#   Si=S.inv()
#   h=Si*one.read('ONEHAMIL','AOONEINT').unpack().unblock()
#   Da=C[:,:na]*C[:,:na].T*S
#   Db=C[:,:nb]*C[:,:nb].T*S
#   Fa,Fb=fab(Da,Db,Si=Si)
#   Fa+=h
#   Fb+=h
#   E=0.5*((Da*(h+Fa)) + (Db*(h+Fb))).tr() + potnuc
#   g=Da*Fa-Fa*Da + Db*Fb-Fb*Db
#   g2=-g*g/2
#   gn=math.sqrt(g2.tr())
#   return E,g,gn
#
# def qnr(C,na,nb,iter=10,hfx=1,threshold=1e-6):
#   C0=C
##
## Initial energy and gradient, unit inverse Hessian
##
#   E0,g0,gn0=Eg(C0,na,nb); H0=full.unit(g0.rdim)
##
#   print "%2d:E=%16.12f %16.5e"%(0,E0,gn0)
##
## Initial direction
##
#   norbt=C.rdim
#   p0=-g0
#   pm0=gradmatrix(p0,nb,na-nb,norbt)
#   G0=E0
#   Gp0=(g0[:,0]&p0)
##
## Main loop
##
#   try:
#      for i in range(iter):
#         print "--- Iteration %d"%(i+1)
#         ## line search
#         x=pm0
#         v=p0
#         U=((-x)).func(cmath.exp)
#         C=C0*U
#         E,g,gn=Eg(C,na,nb)
#         if E < E0:
#            print "%2d:E=%16.12f %16.5e %16.2e"%(i+1,E,gn,E-E0)
#         else:
#            print "%2d:E=%16.12f %16.5e %16.2e"%(i+1,E,gn,E-E0)
#            print "=backstep="
#            G1=E
#            l=-(Gp0)/(2*(G1-G0-Gp0))
#            print "= G0 G1 Gp0 l",G0,G1,Gp0,l
#            x=l*pm0
#            v=l*p0
#            U=(-x).func(cmath.exp)
#            C=C0*U
#            E,g,gn=Eg(C,na,nb)
#            print "%2d:E=%16.12f %16.5e %16.2e"%(i+1,E,gn,E-E0)
#            if E < E0:
#               print "=backstep accept="
#            else:
#               print "=backstep 2="
#               G1=E
#               l=-(Gp0)/(2*(G1-G0-Gp0))
#               print "= G0 G1 Gp0 l",G0,G1,Gp0,l
#               x=-l*pm0
#               v=-l*p0
#               U=(-x).func(cmath.exp)
#               C=C0*U
#               E,g,gn=Eg(C,na,nb)
#               print "%2d:E=%16.12f %16.5e %16.2e"%(i+1,E,gn,E-E0)
#               if E < E0:
#                  print "=backstep 2 accept="
#               else:
#                  print "=backstep 2 fail="
#                  raise BackstepFail(None)
#         if  (gn  < threshold): raise Converged(gn)
##
## Hessian update
##
#         if 1:
#            #
#            # DFP
#            #
#            dg=g-g0
#            Hg=H0*dg
#            gHg=dg&Hg
#            H=H0
#            H += v*v.T/(v&dg)
#            H -= Hg*Hg.T/gHg
#            #
#            # BFGS
#            #
#            u=v/(v&dg) - Hg/gHg
#            H+=gHg*u*u.T
#         elif nasht != 0:
#            raise "open shell not implemented"
#         else:
#            print "no update"
#            H=H0
##
## Next direction
##
#         p=-H*g[:,0]
#         pm=gradmatrix(p,nisht,nasht,norbt)
##
## Save iteration information
##
#         E0=E
#         g0=g
#         H0=H
#         C0=C
#         p0=p
#         pm0=pm
#         Gp0=(g[:,0]&p)
#         G0=E
#
#   except Converged:
#      print "=Converged="
#   except "stop":
#      print "=STOP="
#
# def hinit(n,type="unit"):
#   print "hinit n type",n,type
#   if type == "unit":
#      return full.unit(n)
##      ij=0
##     for i in range(nisht):
##        for j in range(nocct,norbt):
##           print "i j de",i,j,ev[i]-ev[j]
##           de=ev[j]-ev[i]
##           gv[ij,0] /= 2*de
##           gv[ij,1] /= 2*de
#


def diis(
    C,
    nisht,
    nasht,
    iters=10,
    fock=jensen,
    hfx=1,
    threshold=1e-6,
    maxerr=2,
    wrkdir="/tmp",
):
    C1 = C
    E = 0

    aooneint = pathlib.Path(wrkdir) / "AOONEINT"
    aotwoint = pathlib.Path(wrkdir) / "AOTWOINT"
    potnuc = one.readhead(aooneint)["potnuc"]
    vecs = []
    evecs = []
    h = one.read("ONEHAMIL", aooneint).unpack().unblock()
    S = one.read("OVERLAP", aooneint).unpack().unblock()
    try:
        for i in range(iters):
            Di, Da = dens.C2D(C, nisht, nasht)
            D = Di + Da
            Fc = h + two.fock(D, filename=aotwoint, hfx=hfx)
            Fo = two.fock(Da, hfc=0, filename=aotwoint) + Fc
            F = fock(S, Di, Da, Fc, Fo, C, h)
            E0 = E
            E = ((h + Fc) & D) / 2 + ((Fo - Fc) & Da) / 2 + potnuc
            g = grad(S, C1, Di, Da, Fc, Fo)
            # gao = gradao(S, C, Di, Da, Fc, Fo)
            gn = gradvec(S, C, Di, Da, Fc, Fo, nisht, nasht).norm2()
            gn /= math.sqrt(2)
            print("%2d:E = %16.12f %16.5e %16.2e" % (i + 1, E, gn, E - E0))
            if gn < threshold:
                raise Converged(gn)
            vecs.append(F)
            evecs.append(g)
            edim = min(len(evecs), maxerr)
            ev = evecs[-edim:]
            fv = vecs[-edim:]
            B = mkB(ev)
            rhs = full.matrix((edim + 1, 1))
            rhs[-1, 0] = 1
            c = rhs / B
            subevecs = full.matrix(g.shape)
            subvecs = full.matrix(F.shape)
            for i in range(edim):
                subevecs += c[i, 0] * ev[i]
                subvecs += c[i, 0] * fv[i]
            update = -subevecs
            upd = update.lower()
            upd.anti = 0
            update = upd.unpack()
            F = subvecs  # +update
            C = dens.cmo(F, S)
    except Converged:
        print("-Converged-")
    except Stop:
        print("-STOP-")


#
# def fab(Da,Db,Si=None,hfx=1):
#   if Si: # mixed representation in/out
#      Fs=Si*two.fock((Da+Db)*Si,hfx=hfx)
#      Ft=Si*two.fock((Da-Db)*Si,hfx=hfx,hfc=0)
#   else:
#      Fs=two.fock(Da+Db,hfx=hfx)
#      Ft=two.fock(Da-Db,hfx=hfx,hfc=0)
#   Fa=Fs+Ft
#   Fb=Fs-Ft
#   # print ((Da*Fa).tr() + (Db*Fb).tr())/2
#   return Fa,Fb
#


def Feff(Da, Db, Fa, Fb):
    I_n = full.unit(len(Da))
    D = Da + Db
    Ds = Da - Db
    ID = I_n - D
    Fs = Fa - Fb
    F = ((Fa + Fb) + Ds * Fs * ID + ID * Fs * Ds) / 2
    return F


def udiis(
    Ca,
    Cb,
    na,
    nb,
    iters=10,
    fock=jensen,
    hfx=1,
    threshold=1e-6,
    maxerr=2,
    unrest=False,
    wrkdir="/tmp",
):
    print(Ca, Cb)
    saveD = 1
    saveC = 0
    E = 0
    aooneint = os.path.join(wrkdir, "AOONEINT")
    aotwoint = os.path.join(wrkdir, "AOTWOINT")
    potnuc = one.readhead(aooneint)["potnuc"]
    vecs = []
    vecsa = []
    vecsb = []
    evecsa = []
    evecsb = []
    Eit = []
    S = one.read("OVERLAP", aooneint).unpack().unblock()
    h = S.I @ one.read("ONEHAMIL", aooneint).unpack().unblock()
    Da = dens.C1D(Ca, na) @ S
    Db = dens.C1D(Cb, nb) @ S

    try:
        for i in range(iters):
            Da = dens.C1D(Ca, na) @ S
            Db = dens.C1D(Cb, nb) @ S
            print("D", (Da + Db) @ S.I)
            (Fa, Fb), = two.fockab((Da, Db), hfx=hfx, filename=aotwoint)
            Fa = h + S.I @ Fa
            Fb = h + S.I @ Fb
            E0 = E
            E = ((Da @ (h + Fa)) + (Db @ (h + Fb))).tr() / 2 + potnuc
            D = Da + Db
            print("hd", (h @ D).tr(), h & D, h & D.T)
            print("FD", Fa & Da)
            Eit.append(E)
            ga = Da @ Fa - Fa @ Da
            gb = Db @ Fb - Fb @ Db
            if unrest:
                g2 = -(ga @ ga + gb @ gb)
            else:
                g2 = -(ga + gb) @ (ga + gb)
            gn = math.sqrt(g2.tr())
            print("%2d:E = %16.12f %16.5e %16.2e" % (i + 1, E, gn, E - E0))
            if gn < threshold:
                raise Converged(gn)
            # if E > E0:
            #     raise Exception("Energy increase")
            if unrest:
                Ca = dens.cmo(Fa)
                Cb = dens.cmo(Fb)
                # Ca = Ca*Ua
                # Cb = Cb*Ub
            else:
                Ca = dens.cmo(Feff(Da, Db, Fa, Fb), S)
                Cb = Ca[:, :]
            Da = dens.C1D(Ca, na) @ S
            Db = dens.C1D(Cb, nb) @ S
            if saveD:
                vecsa.append(Da)
                vecsb.append(Db)
                evecsa.append(ga @ Da - Da @ ga)
                evecsb.append(gb @ Db - Db @ gb)
            elif saveC:
                vecsa.append(Ca)
                vecsb.append(Cb)
                evecsa.append(ga)
                evecsb.append(gb)
            else:
                vecsa.append(Fa)
                vecsb.append(Fb)
                evecsa.append(ga)
                evecsb.append(gb)
            edim = min(len(evecsa), maxerr)
            eva = evecsa[-edim:]
            evb = evecsb[-edim:]
            fva = vecsa[-edim:]
            fvb = vecsb[-edim:]
            B = mkB3(eva, evb, unrest)
            rhs = full.matrix((edim + 1, 1))
            rhs[-1, 0] = -1
            c = rhs / B
            subvecsa = full.matrix(Fa.shape)
            subvecsb = full.matrix(Fb.shape)
            for j in range(edim):
                subvecsa += c[j, 0] * fva[j]
                subvecsb += c[j, 0] * fvb[j]
            if saveD:
                Da = subvecsa
                Db = subvecsb
                (Fa, Fb), = two.fockab((Da, Db), hfx=hfx, filename=aotwoint)
                Fa = h + S.I @ Fa
                Fb = h + S.I @ Fb
                vecsa[i] = Da
                vecsb[i] = Db
            elif saveC:
                Ca = subvecsa
                Cb = subvecsb
                Da = dens.C1D(Ca, na) @ S
                Db = dens.C1D(Cb, nb) @ S
            else:
                Fa = subvecsa
                Fb = subvecsb
                Da = dens.C1D(Ca, na) @ S
                Db = dens.C1D(Cb, nb) @ S
    except Converged:
        print("Converged after %d iterations\n" % (i + 1,))
    except Increase:
        print("Ca Cb", Ca, Cb)
        print("Da Db", Da, Db)
        print("Na Nb", Da.tr(), Db.tr())
        print("Fa Fb", Fa, Fb)
        print("E1", (h * (Da + Db)).tr())
        print("E2", (Fa * Da + Fb * Db).tr() / 2 - (h * (Da + Db)).tr() / 2)
        print("E", E - potnuc)


if __name__ == "__main__":
    if 1:
        wrkdir = "tests/test_heh.d"
        aooneint = os.path.join(wrkdir, "AOONEINT")
        aotwoint = os.path.join(wrkdir, "AOTWOINT")
        potnuc = one.readhead(aooneint)["potnuc"]
        h = one.read("ONEHAMIL", aooneint).unpack().unblock()
        S = one.read("OVERLAP", aooneint).unpack().unblock()
        Ca = dens.cmo(h, S)
        Cb = Ca

        kwargs = dict(wrkdir=wrkdir, iters=20, threshold=1e-5)
        uroothan(Ca, Cb, 2, 1, unrest=False, **kwargs)
        diis(Ca, 2, 1, **kwargs)
        #udiis(Ca, Cb, 1, 1, **kwargs)

    if 0:
        wrkdir = "tests/test_h2.d"
        aooneint = os.path.join(wrkdir, "AOONEINT")
        aotwoint = os.path.join(wrkdir, "AOTWOINT")
        potnuc = one.readhead(aooneint)["potnuc"]
        h = one.read("ONEHAMIL", aooneint).unpack().unblock()
        S = one.read("OVERLAP", aooneint).unpack().unblock()
        Ca = dens.cmo(h, S)
        Cb = Ca

        uroothan(Ca, Cb, 1, 1, unrest=False, wrkdir=wrkdir, iters=20, threshold=1e-5)
        udiis1(Ca, Cb, 5, 5, wrkdir="tests/test_rohf.d", iter=20, threshold=1e-5)
