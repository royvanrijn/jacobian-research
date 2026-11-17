import pickle
from flint import fmpq,fmpq_poly
import exact_core as ec
S=pickle.load(open('hne0_deg35.pkl','rb'));W=(5,6,5);Q=(6,5,1,0,6,5);U=ec.K(fmpq_poly([0,1]));H=ec.MOD
def kk(c):
 a=[fmpq(0)]*35
 for j,(n,d) in c.items():a[j]=fmpq(n,d)
 return ec.K(fmpq_poly(a))
def tv(k):
 p=k.p;assert all(not p[i] or i%7==0 for i in range(len(p)));return fmpq_poly([p[7*j] if 7*j<len(p) else 0 for j in range(5)])
out=[]
for i,e in enumerate(S):
 ee={}
 for m,c in e.items():
  z=tv(kk(c)*U**(sum(a*b for a,b in zip(W,m))-Q[i]));ee[m]={j:(int(z[j].p),int(z[j].q)) for j in range(len(z)) if z[j]}
 out.append(ee)
pickle.dump(out,open('hne0_v.pkl','wb'))
def stat(X):
 vals=[abs(z) for e in X for c in e.values() for nd in c.values() for z in nd if z]
 return max(len(str(x)) for x in vals),sum(len(str(x)) for x in vals)/len(vals)
print('v stat max/avg digits',stat(out));print('w stat',stat(pickle.load(open('hne0_polred.pkl','rb'))))
print('G coeff stats',[(len(str(abs(int(H[7*k].p)))),len(str(abs(int(H[7*k].q))))) for k in range(6)])
