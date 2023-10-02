import os
import math

from daltools import one, dens
import two

from rohf import *


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
    saveD = 0
    saveC = 0
    E = 0
    aooneint = os.path.join(wrkdir, "AOONEINT")
    aotwoint = os.path.join(wrkdir, "AOTWOINT")
    vecsa = []
    vecsb = []
    evecsa = []
    evecsb = []
    Eit = []
    S = one.read("OVERLAP", aooneint).unpack().unblock()
    h = one.read("ONEHAMIL", aooneint).unpack().unblock()
    potnuc = one.readhead(aooneint)["potnuc"]

    try:
        for i in range(iters):
            Da = dens.C1D(Ca, na)
            Db = dens.C1D(Cb, nb)
            (Fa, Fb), = two.fockab((Da, Db), hfx=hfx, filename=aotwoint)
            Fa += h
            Fb += h

            E0 = E
            E = 0.5 * ((Da & (h + Fa)) + (Db & (h + Fb))) + potnuc
            Eit.append(E)

            ga = S @ Da @ Fa - Fa @ Da @ S
            gb = S @ Db @ Fb - Fb @ Db @ S

            if unrest:
                gn = -((ga & ga) + (gb & gb))
            else:
                gn = (ga + gb) & (S.I @ (ga + gb) @ S.I)

            gn = math.sqrt(2 * (ga + gb) & (S.I @ (ga + gb) @ S.I))
            print("%2d:E = %16.12f %16.5e %16.2e" % (i + 1, E, gn, E - E0))
            breakpoint()
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
                Fa = h + S.I@Fa
                Fb = h + S.I@Fb
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
        # diis(Ca, 1, 0, **kwargs)
        udiis(Ca, Cb, 2, 1, **kwargs)

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
        # udiis1(Ca, Cb, 5, 5, wrkdir="tests/test_rohf.d", iter=20, threshold=1e-5)
