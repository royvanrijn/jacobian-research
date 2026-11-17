from pathlib import Path
import pickle
from flint import fmpq,fmpq_poly
import exact_core as ec
S=pickle.load(open('hne0_branch2_deg35.pkl','rb'));W=(5,6,5);Q=(6,5,1,0,6,5);U=ec.K(fmpq_poly([0,1]));F=fmpq_poly([26,0,3,3,-1,1]);PHI=fmpq_poly([fmpq(-9725570295901,12623962),fmpq(-1170753213563,971074),fmpq(-387111042229,12623962),fmpq(1578225240619,12623962),fmpq(-469713794365,6311981)])
def ev(p,x,mod):
 o=fmpq_poly([0])
 for i in range(len(p)-1,-1,-1):o=(o*x+fmpq_poly([p[i]]))%mod
 return o
H=ec.MOD;G=fmpq_poly([H[7*k] for k in range(6)]);assert ev(G,PHI,F).is_zero()
def kk(c):
 a=[fmpq(0)]*35
 for j,(n,d) in c.items():a[j]=fmpq(n,d)
 return ec.K(fmpq_poly(a))
def tw(k):
 p=k.p;assert all(not p[i] or i%7==0 for i in range(len(p)));v=fmpq_poly([p[7*j] if 7*j<len(p) else 0 for j in range(5)]);return ev(v,PHI,F)
out=[]
for i,e in enumerate(S):
 ee={}
 for m,c in e.items():
  z=tw(kk(c)*U**(sum(a*b for a,b in zip(W,m))-Q[i]));d={j:(int(z[j].p),int(z[j].q)) for j in range(len(z)) if z[j]}
  if d:ee[m]=d
 out.append(ee)
pickle.dump(out,open('hne0_branch2_polred.pkl','wb'),protocol=4);print([len(e) for e in out])
