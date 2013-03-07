import numpy as np
from util.unformatted import FortranBinary as FB
from util.full import matrix

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


    aofile = FB(filename) 
    aofile.find("BASTWOEL")
    J = matrix(D.shape)
    Ka = matrix(D.shape)
    Kb = matrix(D.shape)

    for rec in aofile:
        #
        #len = 8*lbuf + 4*lbuf + 4
        #lbuf is paramter interger : 600
        lbuf = (aofile.reclen-4)/12

        buf = np.array(aofile.readbuf(lbuf,'d'))
        ibuf = np.array(aofile.readbuf(4*lbuf,'b')).reshape((lbuf, 4))
        length = aofile.readbuf(1,'i')[0]
        #
        # Negative length marks end of file
        #
        if length < 0: break

        if f2py:
            J, Ka, Kb = sirfck.fckab(
                J, Ka, Kb, Da, Db, buf[:length], ibuf.T[:, :length]
                )
        else:
            for g, ig in zip(buf[:length], ibuf[:length]):
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

if __name__ == "__main__":
    from util.full import unit
    D = matrix((6, 6))
    D[0, 0] = 1.0
    Fa, Fb = fockab((D, D), "test/test_fockab.d/AOTWOINT", f2py=False)
    print Fa
    Fa, Fb = fockab((D, D), "test/test_fockab.d/AOTWOINT", f2py=True)
    print Fa
