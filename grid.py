#!/usr/bin/env python
import math
import numpy as np
import daltools
from daltools.util import full

def readfirst(quad_file):
    for f in quaditer(quad_file):
        return f

def quaditer(quad_file):
    """Loop over quadratuire file and retyrn x,y,z, weight"""
    import struct
    isize = 4
    dsize = 8
    with open(quad_file, 'rb') as qf:
        while True:
            ibytes = qf.read(2*isize)
            if not ibytes: raise StopIteration
            n, nb = struct.unpack(2*'i', ibytes)
            #print n
            ibytes = qf.read(2*nb*isize)
            fbytes = qf.read(4*n*dsize)
            xyzw = struct.unpack(4*n*'d', fbytes)
            for i in range(n):
                yield xyzw[3*i], xyzw[3*i+1], xyzw[3*i+2], xyzw[3*n+i]
        
def generate_primitive(l, e, X):

    # Primitive Normalization Factors
    # For s and p all components have the same normalization


    if l == 0:
        def S(alpha, beta):
            return (math.pi/(alpha+beta))**1.5 
        def Norm(alpha):
            return S(alpha,alpha)**(-0.5)
        Na = [Norm(alpha) for alpha in e]

        # return func
        def f(*args):
            x, y, z = args
            if type(x) == float:
                from math import exp
            else:
                from numpy import exp
            import math
            dx = (x - X[0])
            dy = (y - X[1])
            dz = (z - X[2])
            dr = (dx, dy, dz)
            dr2 = dx**2 + dy**2 + dz**2
            s = [ N_*exp(-e_*dr2) for N_, e_ in zip(Na, e) ]
            sgrad = [[-2*e_*x*N_*exp(-e_*dr2) for x in dr] for N_, e_ in zip(Na, e)]
            return s, sgrad

        f.grad = None

        return f

    elif l == 1:
        def S(alpha, beta):
            return (math.pi/(alpha+beta))**1.5 * 1/(2*(alpha+beta))
        def Norm(alpha):
            return S(alpha,alpha)**(-0.5)
        Na = [Norm(alpha) for alpha in e]

        # return func
        def f(*args):
            import math
            x, y, z = args
            dx = (x - X[0])
            dy = (y - X[1])
            dz = (z - X[2])
            dr = (dx, dy, dz)
            dr2 = dx**2 + dy**2 + dz**2
            p = [ N_*c*math.exp(-e_*dr2)  for N_, e_ in zip(Na, e) for c in dr]
            pgrad = [[-2*e_*x*N_*math.exp(-e_*dr2) for x in dr] for N_, e_ in zip(Na, e)]
            return p, pgrad

        return f

    else:
        raise Exception("l>1 not implemented")

def generate_contracted(l, e, c, X):

    # Primitive Normalization Factors
    # For s and p all components have the same normalization

    if l < 2:
        def S(alpha, beta, l):
            return (math.pi/(alpha+beta))**1.5 * 1/(2*(alpha+beta))**(l)
        def Norm(alpha, l):
            return S(alpha,alpha, l)**(-0.5)
        Na = [Norm(alpha, l) for alpha in e]
        #
        # Normalize contracted
        #
        Sc = 0
        for N1, e1, c1  in zip(Na, e, c):
            for N2, e2, c2 in zip(Na, e, c):
                Sc += c1*c2*N1*N2*S(e1, e2, l)
        Nc = Sc**(-0.5)
        cc = [Nc*c_ for c_ in c]

    if l == 0:
            
        # return func
        def f(*args):
            import math
            x, y, z = args
            dx = (x - X[0])
            dy = (y - X[1])
            dz = (z - X[2])
            dr2 = dx**2 + dy**2 + dz**2
            s = 0.0
            sx = 0.0
            sy = 0.0
            sz = 0.0
            for N, e_, c_  in zip(Na, e, cc):
                stmp  = c_*N*math.exp(-e_*dr2)
                s += stmp
                stmp *= -2*e_
                sx += dx*stmp
                sy += dy*stmp
                sz += dz*stmp

            return (s,), ((sx, sy, sz),)
        return f
    elif l == 1:

        def f(*args):
            import math
            x, y, z = args
            dx = (x - X[0])
            dy = (y - X[1])
            dz = (z - X[2])
            dr2 = dx**2 + dy**2 + dz**2
            px = 0
            py = 0
            pz = 0
            pxx = 0; pxy = 0; pxz = 0 
            pyx = 0; pyy = 0; pyz = 0 
            pzx = 0; pzy = 0; pzz = 0 
    
            # sum primitive contributions
            for N, e_, c_  in zip(Na, e, cc):
                expa = c_*N*math.exp(-e_*dr2)
                px += dx*expa
                py += dy*expa
                pz += dz*expa

                pxx += expa + dx*expa*(-2*e_*dx)
                pxy +=        dx*expa*(-2*e_*dy)
                pxz +=        dx*expa*(-2*e_*dz)

                pyx +=        dy*expa*(-2*e_*dx)
                pyy += expa + dy*expa*(-2*e_*dy)
                pyz +=        dy*expa*(-2*e_*dz)

                pzx +=        dz*expa*(-2*e_*dx)
                pzy += expa + dz*expa*(-2*e_*dy)
                pzz +=        dz*expa*(-2*e_*dz)

            return (px, py, pz), ((pxx, pxy, pxz), (pyx, pyy, pyz), (pzx, pzy, pzz))
        return f
    elif l == 2:

        def f(*args):
            import math
            x, y, z = args
            dx = (x - X[0])
            dy = (y - X[1])
            dz = (z - X[2])
            dr2 = dx**2 + dy**2 + dz**2
            dxx = 0
            dxy = 0
            dxz = 0
            dyy = 0
            dyz = 0
            dzz = 0
            # sum primitive contributions
            for N, e_, c_  in zip(Na, e, c):
                expa = N*math.exp(-e_*dr2)
                dxx += c_*dx*dx*expa
                dxy += c_*dx*dy*expa
                dxz += c_*dx*dz*expa
                dyy += c_*dy*dy*expa
                dyz += c_*dy*dz*expa
                dzz += c_*dz*dz*expa
            return (dxx, dxy, dxz, dyy, dyz, dxx)
        return f
    else:
        raise Exception("l>0 not implemented")
                
def set_cfunc(bf):
    """Return a list of contracted AO function objects

    """
    from daltools import mol
    molecule = mol.readin(bf)

    cfunc = []
    for atom in molecule:
        for center in atom['center']:
            for l, lblock in enumerate(atom['basis']):
                for block in lblock:
                    #pe = np.array(block['prim'])
                    pe = block['prim']
                    cc = np.array(block['cont'])
                    for c in range(cc.shape[1]):
                        cfunc.append(generate_contracted(l, pe, cc[:, c], center) )
    return cfunc

def eval_cfunc(cfunc, x):
    res = []
    for f in cfunc:
        fx, gx = f(*x)
        res += list(fx)

    return full.init(res)


def count_electrons(bas_file, quad_file):
    from daltools import dens
    nel = 0
    cfunc = set_cfunc(bas_file)
    Da, Db = dens.Dab()
    D = Da + Db
    for x, y, z, w in quaditer(quad_file):
        f = eval_cfunc(cfunc, (x, y, z))
        n = f&(D*f)
        nel += w*n
    return nel
            

def integrator(*args, **kwargs):
    "Integrate list of args, defined as functionals of density"
    from daltools import dens
    bf=kwargs['bf']
    qf=kwargs['qf']
    cfunc = set_cfunc(bf)
    Da, Db = dens.Dab()
    res = [0 for a in args]
    for x, y, z, w in quaditer(qf):
        cx = eval_cfunc(cfunc, (x, y, z))
        rhoax = cx&(Da*cx)
        rhobx = cx&(Db*cx)
        for i, a in enumerate(args):
            res[i] += w*a(rhoax, rhobx)
    return res

class GridIntegrator():
    def __init__(self, quad_file='DALTON.QUAD'):
        points = np.array(list(quaditer(quad_file)))
        self.x = points[:,0]
        self.y = points[:,1]
        self.z = points[:,2]
        self.w = points[:,3]
        self.points = points
    def integrate(self, f):
        vf = f(self.x, self.y, self.z)
        return np.dot(self.w, vf)
                

if __name__ == "__main__":
    import sys
    sumw = 0
    sexp = 0
    a=1
    b=1
    l=0
    for x,y,z,w in quaditer(sys.argv[1]):
        #print x,y,z,w
        sumw += w
        sexp += w*math.exp(-(a+b)*(x*x+y*y+z*z))
    sexp *= ((a+b)/math.pi)**1.5 #* 2*(a+b)

    print "Sum of weights", sumw
    print "Integral", sexp
    
