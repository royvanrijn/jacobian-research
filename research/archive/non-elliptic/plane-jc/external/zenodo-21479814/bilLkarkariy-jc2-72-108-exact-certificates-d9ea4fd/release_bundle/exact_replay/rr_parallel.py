import pickle,sys,os
from flint import fmpz
from multiprocessing import Pool
base,M,xs=pickle.load(open('crt_cache.pkl','rb'))
def rr_one(x):
 x%=M;B=(M//2).isqrt();r0,r1=M,x;t0,t1=fmpz(0),fmpz(1)
 while abs(r1)>B:q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
 a,b=r1,t1
 if b<0:a,b=-a,-b
 if not b or abs(a)>B or b>B or a.gcd(b)!=1 or (x*b-a)%M:return None
 return int(a),int(b)
with Pool(min(24,os.cpu_count() or 4)) as P:res=list(P.imap(rr_one,xs,chunksize=8))
pickle.dump((base,res),open('rr_cache.pkl','wb'),protocol=4)
print('reconstructed',sum(x is not None for x in res),'missing',sum(x is None for x in res))
