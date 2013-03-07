#!/usr/bin/env python
"""A lot of cleanup to do"""

#import dens,one,two,prop,full,sirifc,math,cmath
#class gradient(full.matrix):
#   def __init__(self,nrow,ncol):
#      full.matrix.__init__(self,nrow,ncol)
#      self.nisht=0
#      self.nasht=0
#      self.norbt=0
#   def init(self,S,C,Dc,Do,Fc,Fo):
#      self.nisht=int((S&Dc)+0.5)/2
#      self.nasht=int((S&Do)+0.5)
#      self.norbt=C.cdim
#      g=2*C.T()*S*(Dc*Fc+Do*Fo)*C
#      g=g-g.T()
#      self.data=g.data
#      return self
#def grad(S,C,Dc,Do,Fc,Fo):
#   #print "C",C
#   #print "4CS",4*C.T()*S
#   #print "Dc",Dc
#   #print "Fc",Fc
#   #print "Dc*Fc",Dc*Fc
#   #print "Do",Do
#   #print "Fo",Fo
#   #print "Dc*Fc+DoFo",Dc*Fc+Do*Fo
#   g=2*C.T()*S*(Dc*Fc+Do*Fo)*C
#   #print "grad:2*g",2*g
#   g=g-g.T()
#   #print "grad:g-gt",g
#   #print "grad.Fc,Fo",Fc.lower(),Fo.lower()
#   #print "grad.g",g.lower()
#   return g
#def gradao(S,C,Dc,Do,Fc,Fo):
#   g=2*S*(Dc*Fc+Do*Fo)
#   #print "gradao:2*g",2*g
#   g=g-g.T()
#   #print "gradao:g-gt",g
#   #print "grad.Fc,Fo",Fc.lower(),Fo.lower()
#   #print "grad.g",g.lower()
#   return g
#def gradvec(S,C,Dc,Do,Fc,Fo,nisht=0,nasht=0,nocct=0):
#   g=grad(S,C,Dc,Do,Fc,Fo)
#   nisht=int((S&Dc)/2 + .5)
#   nasht=int((S&Do)   + .5)
#   norbt=C.cdim
#   nocct=nisht+nasht
#   #print nisht,nasht,nocct,norbt
#   dim=nisht*(norbt-nisht)+nasht*(norbt-nocct)
#   gv=full.matrix(dim,2)
#   ij=0
#   for i in range(nisht):
#      for j in range(nisht,nocct):
#         #print ij+1,i+1,j+1
#         gv[ij,0]=-g[i,j]
#         gv[ij,1]=-g[j,i]
#         ij+=1
#   for i in range(nisht):
#      for j in range(nocct,norbt):
#         #print ij+1,i+1,j+1
#         gv[ij,0]=-g[i,j]
#         gv[ij,1]=-g[j,i]
#         ij+=1
#   for i in range(nisht,nocct):
#      for j in range(nocct,norbt):
#         #print ij+1,i+1,j+1
#         gv[ij,0]=-g[i,j]
#         gv[ij,1]=-g[j,i]
#         ij+=1
#   return gv
#def gradmatrix(p,nisht,nasht,norbt):
#   new=full.matrix(norbt,norbt)
#   nocct=nisht+nasht
#   ij=0
#   for i in range(nisht):
#      for j in range(nisht,norbt):
#         new[i,j]=-p[ij,0]
#         new[j,i]=p[ij,0]
#         ij+=1
#   for i in range(nisht,nocct):
#      for j in range(nocct,norbt):
#         new[i,j]=-p[ij,0]
#         new[j,i]=p[ij,0]
#         ij+=1
#   return new
#def gradnorm(S,C,Dc,Do,Fc,Fo,nisht,nasht):
#   import math
#   g=grad(S,C,Dc,Do,Fc,Fo)
#   norbt=C.cdim
#   gsum=0.0
#   for i in range(nisht):
#      for j in range(nisht,norbt):
#         gsum+=g[i,j]**2
#   for i in range(nisht,nisht+nasht):
#      for j in range(nisht+nasht,norbt):
#         gsum+=g[i,j]**2
#   return math.sqrt(gsum)
#def jensen(S,Dc,Do,Fc,Fo,C=None,h=None):
#   I=full.unit(S.rdim)
#   corr=S*Do*(Fc-Fo)*((Dc+Do)*S-I)
#   F=Fc+corr+corr.T()
#   return F
#def jensen_scaled(S,Dc,Do,Fc,Fo):
#   I=full.unit(S.rdim)
#   corr=S*Dc*(Fc-Fo)*((Dc+Do)*S-I)
#   F=Fc
#   for i in range(nisht):
#      for j in range(nisht+nasht,F.rdim):
#         F[i,j]*=2
#         F[j,i]*=2
#   return F+corr+corr.T()
#def filatov(S,Dc,Do,Fc,Fo):
#   I=full.unit(S.rdim)
#   f=0.5
#   beta=1./(1-f)
#   PS = Do*S
#   SG = S*Dc/2 - (1/beta)*I + (1/(2*beta))*S*Do
#   F=Fc
#   corr=beta*SG*(Fc-Fo)*PS
#   F+=corr+corr.T()
#   return F
#def okazaki(S,Dc,Do,Fc,Fo):
#   Dc=0.5*Dc
#   Fo=0.5*Fo
#   I=full.unit(S.rdim)
#   F=(I-S*Do)*Fc*(I-Do*S) + (I-S*Dc)*Fo*(I-Dc*S)\
#    +S*(Dc*(Fc-Fo)*Do + Do*(Fc-Fo)*Dc)*S
#   return F
#   

class Converged(Exception):
    """To be removed"""

    def __init__(self, value):
        Exception.__init__(self)
        self.value = value

    def __str__(self):
        return repr(self.value)

#class BackstepFail(Exception):
#   def __init__(self,value):
#      self.value=value
#   def __str__(self):
#      return repr(self.value)

from math import sqrt
from util import full
from daltools import one, dens
from dalmisc.fockab import fockab

def energy(Da, Db, h1, Fa, Fb):    
    e1 = h1&(Da+Db)
    e2 = 0.5*((Da&Fa) + (Db&Fb))
    return e1 + e2
    
def uroothan(Ca, Cb, na, nb, hfx=1, iters=10, threshold=1e-6, unrest=False):
    """Open-shell Roothan iterations, restricted or unrestricted"""
    if (unrest):
        print "Unrestricted HF Na=%d Nb=%d\n" % (na, nb)
    else:
        print "Restricted RHF Na=%d Nb=%d\n" % (na, nb)
    E0 = 0.0
    h = one.read('ONEHAMIL','AOONEINT').unpack().unblock()
    S = one.read('OVERLAP','AOONEINT').unpack().unblock()
    I = full.unit(Ca.rdim)
    Si = S.I
    h = Si*h
    potnuc = one.readhead("AOONEINT")["potnuc"]
    iterinf = []
    try:
        for i in range(iters):
            Da = dens.C1D(Ca, na)
            Db = dens.C1D(Cb, nb)
            Fa, Fb = fockab(Da, Db, hfx=hfx)
            Fa += h
            Fb += h
            E = 0.5*((Da&(h+Fa)) + (Db&(h+Fb))) + potnuc
            ga = Da*Fa - Fa*Da
            gb = Db*Fb - Fb*Db
            if (unrest):
                g2 = -(ga*ga+gb*gb)
            else:
                g = ga+gb
                g2 = -(ga+gb)*(ga+gb)/2
            gn = sqrt(g2.tr())
            iterinf.append((E, gn))
            print "%2d:E=%16.12f %16.5e %16.2e" % (i+1, E, gn, E-E0)
            if  (gn < threshold): raise Converged(gn)
            if unrest:
                Ca = dens.cmo(Fa, S)
                Cb = dens.cmo(Fb, S)
            else:
                D = Da+Db
                Ds = Da-Db
                Fs = Fa-Fb
                ID = I-D
                F = ((Fa+Fb) + Ds*Fs*ID + ID*Fs*Ds)/2
                Ca = dens.cmo(F, S)
                Cb = Ca
    except Converged:
        print "-Converged-"
    if (unrest):
        return (Ca, Cb)
    else:
        return Ca

#def mkB(vecs):
#   B=full.matrix(len(vecs)+1,len(vecs)+1)
#   for i in range(len(vecs)):
#      for j in range(len(vecs)):
#         B[i,j]=vecs[i]&vecs[j]
#      B[i,len(vecs)]=1
#      B[len(vecs),i]=1
#   return B
#def mkB2(vecs):
#   B=full.matrix(len(vecs)+1,len(vecs)+1)
#   for i in range(len(vecs)):
#      for j in range(len(vecs)):
#         B[i,j]=-(vecs[i]*vecs[j]).tr()
#      B[i,len(vecs)]=-1
#      B[len(vecs),i]=-1
#   return B
#def mkB3(avecs,bvecs,unrest):
#   dim=len(avecs)
#   B=full.matrix(dim+1,dim+1)
#   for i in range(dim):
#      for j in range(dim):
#         if unrest:
#            B[i,j]=-(avecs[i]*avecs[j]).tr() - (bvecs[i]*bvecs[j]).tr()
#         else:
#            vecsi=avecs[i]+bvecs[i]
#            vecsj=avecs[j]+bvecs[j]
#            B[i,j]=(vecsi*vecsj.T()).tr()
#      B[i,dim]=-1
#      B[dim,i]=-1
#   return B
#
#def Eg(C,na,nb,hfx=1):
#   import one
#   potnuc=one.info()[0]["potnuc"]
#   S=one.read('OVERLAP','AOONEINT').unpack().unblock()
#   Si=S.inv()
#   h=Si*one.read('ONEHAMIL','AOONEINT').unpack().unblock()
#   Da=C[:,:na]*C[:,:na].T()*S
#   Db=C[:,:nb]*C[:,:nb].T()*S
#   Fa,Fb=fab(Da,Db,Si=Si)
#   Fa+=h
#   Fb+=h
#   E=0.5*((Da*(h+Fa)) + (Db*(h+Fb))).tr() + potnuc
#   g=Da*Fa-Fa*Da + Db*Fb-Fb*Db
#   g2=-g*g/2
#   gn=math.sqrt(g2.tr())
#   return E,g,gn
#
#def qnr(C,na,nb,iter=10,hfx=1,threshold=1e-6):
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
#            H += v*v.T()/(v&dg)
#            H -= Hg*Hg.T()/gHg
#            #
#            # BFGS
#            #
#            u=v/(v&dg) - Hg/gHg
#            H+=gHg*u*u.T()
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
#def hinit(n,type="unit"):
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
#def diis(C,nisht,nasht,iter=10,fock=jensen,hfx=1,threshold=1e-6,maxerr=2):
#   import one
#   C1=C
#   E=0
#   potnuc=one.info()[0]["potnuc"]
#   vecs=[]
#   evecs=[]
#   h=one.read('ONEHAMIL','AOONEINT').unpack().unblock()
#   S=one.read('OVERLAP','AOONEINT').unpack().unblock()
#   I=full.unit(S.rdim)
#   try:
#      for i in range(iter):
#         print "--- Iteration %d"%(i+1)
#         Di,Da=dens.C2D(C,nisht,nasht)
#         D=Di+Da
#         Fc=h+two.fock(D,hfx)
#         Fo=two.fock(Da,hfc=0)+Fc
#         F=fock(S,Di,Da,Fc,Fo,C,h)
#         #print "F",(C.T()*F*C).pack()
#         E0=E
#         E=((h+Fc)&D)/2+((Fo-Fc)&Da)/2 + potnuc
#         g=grad(S,C1,Di,Da,Fc,Fo)
#         gao=gradao(S,C,Di,Da,Fc,Fo)
#         gn=gradvec(S,C,Di,Da,Fc,Fo,nisht,nasht)[:,0].norm2()
#         print "%2d:E=%16.12f %16.5e %16.2e"%(i+1,E,gn,E-E0)
#         if  (gn  < threshold): raise "converged"
#         vecs.append(F)
#         evecs.append(g)
#         #evecs.append(gao)
#         #print "vecs",vecs
#         #print "evecs",evecs
#         edim=min(len(evecs),maxerr)
#         ev=evecs[-edim:]
#         print "ev",ev
#         fv=vecs[-edim:]
#         B=mkB(ev)
#         print "B",B
#         rhs=full.matrix(edim+1,1)
#         rhs[-1,0]=1
#         #print "B",B
#         #print "B/2",B/2
#         print "rhs",rhs
#         c=rhs/B
#         print "c",c
#         subevecs=full.matrix(g.rdim,g.cdim)
#         subvecs=full.matrix(F.rdim,F.cdim)
#         for i in range(edim):
#            #print "diis average %d(%d)"%(i+1,edim),c[i,0],fv[i]
#            subevecs+=c[i,0]*ev[i]
#            subvecs+=c[i,0]*fv[i]
#         if 0:
#            #TEST
#            # sum previous
#            # add current correction
#             print "last fv",(fv[edim-1]-h)
#         #evecs[edim]=
#         update=-subevecs
#         upd=update.lower()
#         upd.anti=0
#         update=upd.unpack()
#         #print 'update',C.T()*update*C
#         F=subvecs#+update
#         #print "F(ao)",F.lower()
#         #print "F(mo)",(C.T()*F*C).lower()
#         #print "F(eig)",(C.T()*F*C).eig()
#         C=dens.cmo(F,S)
#         print "C",C
#   except "converged":
#      print "-Converged-"
#   except "stop":
#      print "-STOP-"
#
#def fab(Da,Db,Si=None,hfx=1):
#   if Si: # mixed representation in/out
#      Fs=Si*two.fock((Da+Db)*Si,hfx=hfx)
#      Ft=Si*two.fock((Da-Db)*Si,hfx=hfx,hfc=0)
#   else:
#      Fs=two.fock(Da+Db,hfx=hfx)
#      Ft=two.fock(Da-Db,hfx=hfx,hfc=0)
#   Fa=Fs+Ft
#   Fb=Fs-Ft
#   #print ((Da*Fa).tr() + (Db*Fb).tr())/2
#   return Fa,Fb
#
#def Feff(Da,Db,Fa,Fb):
#   I=full.unit(Da.rdim)
#   D=Da+Db
#   Ds=Da-Db
#   ID=I-D
#   Fs=Fa-Fb
#   F=((Fa+Fb) + Ds*Fs*ID + ID*Fs*Ds)/2
#   return F
#         
#
#def udiis1(Ca,Cb,na,nb,iter=10,fock=jensen,hfx=1,threshold=1e-6,maxerr=2,unrest=0):
#   #class Converged(Exception)
#   saveD=1
#   saveC=0
#   import one
#   E=0
#   potnuc=one.info()[0]["potnuc"]
#   vecs=[]
#   vecsa=[]
#   vecsb=[]
#   evecs=[]
#   evecsa=[]
#   evecsb=[]
#   Eit=[]
#   S=one.read('OVERLAP','AOONEINT').unpack().unblock()
#   Si=1/S
#   Si=Ca*Ca.T()
#   #print Si*S
#   h=Si*one.read('ONEHAMIL','AOONEINT').unpack().unblock()
#   I=full.unit(S.rdim)
#   Da=dens.C1D(Ca,na)*S
#   Db=dens.C1D(Cb,nb)*S
#   C0=Ca[:,:]
#   p0=0*S
#   #print "D(mo)",Ca.T()*S*Da*Ca,Cb.T()*S*Db*Cb
#
#   try:
#      for i in range(iter):
#         #print "--- Iteration %d"%(i+1)
#         #Da=dens.C1D(Ca,na)*S
#         #Db=dens.C1D(Cb,nb)*S
#         #print Da.tr(),Db.tr()
#         Fa,Fb=fab(Da,Db,Si=Si,hfx=hfx)
#         Fa=h+Fa
#         Fb=h+Fb
#         E0=E
#         E=(Da*(h+Fa) + Db*(h+Fb)).tr()/2 + potnuc
#         Eit.append(E)
#         ga=Da*Fa-Fa*Da
#         gb=Db*Fb-Fb*Db
#         if unrest:
#            g2=-(ga*ga+gb*gb)
#         else:
#            g2=-(ga+gb)*(ga+gb)/2
#         gn=math.sqrt(g2.tr())
#         print "%2d:E=%16.12f %16.5e %16.2e"%(i+1,E,gn,E-E0)
#         #print "D(mo)",Ca.T()*S*Da*Ca,Cb.T()*S*Db*Cb
#         if  (gn  < threshold): raise Converged(gn)
#         #if  (E  > E0): raise Exception('Energy increase')
#         if unrest:
#            Ca=dens.cmo(Fa)
#            Cb=dens.cmo(Fb)
#            #Ca=Ca*Ua
#            #Cb=Cb*Ub
#         else:
#            Ca=dens.cmo(Feff(Da,Db,Fa,Fb))
#            Cb=Ca[:,:]
#            #print Ca.T()*S*Ca
#         Da=dens.C1D(Ca,na)*S
#         Db=dens.C1D(Cb,nb)*S
#         if saveD:
#            vecsa.append(Da)
#            vecsb.append(Db)
#            evecsa.append(ga*Da - Da*ga)
#            evecsb.append(gb*Db - Db*gb)
#         elif saveC:
#            vecsa.append(Ca)
#            vecsb.append(Cb)
#            evecsa.append(ga)
#            evecsb.append(gb)
#         else:
#            vecsa.append(Fa)
#            vecsb.append(Fb)
#            evecsa.append(ga)
#            evecsb.append(gb)
#         edim=min(len(evecsa),maxerr)
#         eva=evecsa[-edim:]
#         evb=evecsb[-edim:]
#         fva=vecsa[-edim:]
#         fvb=vecsb[-edim:]
#         B=mkB3(eva,evb,unrest)
#         rhs=full.matrix(edim+1,1)
#         rhs[-1,0]=-1
#         #print "rhs",rhs
#         c=rhs/B
#         #print "c",c
#         subvecsa=full.matrix(Fa.rdim,Fa.cdim)
#         subvecsb=full.matrix(Fb.rdim,Fb.cdim)
#         for j in range(edim):
#            subvecsa+=c[j,0]*fva[j]
#            subvecsb+=c[j,0]*fvb[j]
#         if saveD:
#            Da=subvecsa
#            #print Da*Da-Da
#            Db=subvecsb
#            #print Da.tr(),Db.tr()
#            #Fa,Fb=fab(Da,Db,Si,hfx)
#            #Fa=h+Fa
#            #Fb=h+Fb
#            #vecsa[i]=Da
#            #vecsb[i]=Db
#            #print "vecsa",vecsa
#         elif saveC:
#            Ca=subvecsa
#            Cb=subvecsb
#            Da=dens.C1D(Ca,na)*S
#            Db=dens.C1D(Cb,nb)*S
#         else:
#            Fa=subvecsa
#            Fb=subvecsb
#            Da=dens.C1D(Ca,na)*S
#            Db=dens.C1D(Cb,nb)*S
#   except Converged:
#      print "Converged after %d iterations\n"%(i+1,)
##  except Increase:
##        print "Ca Cb",Ca,Cb
##        print "Da Db",Da,Db
##        print "Na Nb",Da.tr(), Db.tr()
##        print "Fa Fb",Fa,Fb
##        print "E1",(h*(Da+Db)).tr()
##        print "E2",(Fa*Da+Fb*Db).tr()/2-(h*(Da+Db)).tr()/2
##        print "E",E-potnuc
#   except "stop":
#      print "-STOP-"
##  for j in range(len(Eit)):
##     print "%d %10.6f %10.6f"%(j+1,Eit[j],Eit[j]-Eit[-1])


if __name__ == "__main__":
    pass
