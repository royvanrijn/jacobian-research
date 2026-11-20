#!/usr/bin/env python3
"""Complete actual a=2 finite-pair replay, independent recurrence and prime.

For each m the normalized physical moment is
T_m(t)=sum binom(m,k)(2m+1)_(d*k)t^k. BOTH actual shapes are
constructed separately: upper2N+1, previous2N-1. No fixed-shape
discriminant identity is used. A complete finite PASS is NOT an infinite
all-terminal theorem until an independently justified tail is supplied.
"""
if not __debug__:raise SystemExit('Unoptimized Python required for exact checks')
from pathlib import Path
from math import isqrt,comb,factorial
from hashlib import sha256
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,json,struct,time
import flint
from flint import nmod_poly
PRIMES=(1000033,1000037,1000039,1000003)

def tables(d,limit,p):
 assert 1<=d<=limit<p and p>1 and all(p%j for j in range(2,isqrt(p)+1))
 inv=[0]*(limit+1);inv[1]=1
 for j in range(2,limit+1):inv[j]=-(p//j)*inv[p%j]%p
 rise=[0]*(limit+1);q=1
 for j in range(1,d+1):q=q*j%p
 rise[1]=q
 for t in range(1,limit-d+1):
  q=q*(t+d)%p*inv[t]%p;rise[t+1]=q
 return inv,rise

def coefficients(m,d,p,tab):
 inv,rise=tab;c=[1];v=1
 for k in range(m):
  v=v*(m-k)%p*inv[k+1]%p*rise[2*m+1+d*k]%p
  c.append(v)
 return c

def worker(args):
 d,last=args;start=time.monotonic();limit=(d+2)*last
 cache={};counts={};fallback=[];h=sha256();checks=0
 for N in range(2,last+1):
  attempts=[]
  for p in PRIMES:
   if p not in cache:cache[p]=tables(d,limit,p)
   a=coefficients(N,d,p,cache[p]);b=coefficients(N-1,d,p,cache[p])
   assert len(a)==N+1 and len(b)==N and a[0]==b[0]==1 and a[-1]and b[-1]
   f=nmod_poly(a,p);g=nmod_poly(b,p)
   assert f.degree()==N and g.degree()==N-1
   q,s,t=f.xgcd(g);assert s*f+t*g==q;checks+=1
   attempts.append([p,q.degree()])
   if q.degree()==0:
    assert q==nmod_poly([1],p);counts[str(p)]=counts.get(str(p),0)+1
    h.update(struct.pack('<III',d,N,p))
    if len(attempts)>1:fallback.append({'terminal':N,'attempts':attempts})
    break
  else:raise AssertionError(('Unresolved physical a2 finite pair',d,N,attempts))
 return {'gap':d,'first_terminal':2,'last_terminal':last,'tested':last-1,
         'verified_full_Bezout_identities':checks,'successful_primes':counts,
         'fallbacks':fallback,'maximum_upper_degree':last,'maximum_previous_degree':last-1,
         'case_stream_sha256':h.hexdigest(),'seconds':time.monotonic()-start}

def controls():
 count=0
 for d in [1,2,6,10,45]:
  p=PRIMES[0];t=tables(d,10*(d+2),p)
  for N in range(2,10):
   for m in [N,N-1]:
    c=coefficients(m,d,p,t)
    literal=[comb(m,k)*factorial(2*m+d*k)//factorial(2*m)%p for k in range(m+1)]
    assert c==literal;count+=1
 # Explicit negative for replacing upper Gamma2N+1 by previous Gamma2N-1.
 N=3;d=2;c=coefficients(N,d,PRIMES[0],tables(d,15,PRIMES[0]))
 assert c[1]==168 and comb(N,1)*(2*N-1)*(2*N)==90!=c[1]
 p=PRIMES[0];f=nmod_poly([1,1],p);g=f*f;q,s,t=f.xgcd(g)
 assert q.degree()==1 and s*f+t*g==q
 assert (s+1)*f+t*g!=q
 return count

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=4)
 ap.add_argument('--output',type=Path,default=Path('/tmp/independent_lower_two_practical_finite_box.json'))
 ap.add_argument('--benchmark-last',type=int)
 a=ap.parse_args();start=time.monotonic();cc=controls()
 expected=[(d,max(4000,128*(d+2))-1)for d in range(6,46)]
 assert len(expected)==40 and sum(last-1 for d,last in expected)==176816
 assert max((d+2)*last for d,last in expected)==282705
 domain=expected if a.benchmark_last is None else[(d,min(last,a.benchmark_last))for d,last in expected]
 rows=[]
 with ProcessPoolExecutor(max_workers=a.jobs)as pool:
  fs={pool.submit(worker,c):c for c in domain}
  for f in as_completed(fs):
   v=f.result();rows.append(v);print(json.dumps(v),flush=True)
   a.output.with_suffix('.running.json').write_text(json.dumps({'status':'RUNNING','completed_gaps':sorted(rows,key=lambda r:r['gap'])},indent=2)+'\n')
 rows.sort(key=lambda r:r['gap']);total=sum(r['tested']for r in rows);full=a.benchmark_last is None
 if full:
  assert [(r['gap'],r['last_terminal'])for r in rows]==expected and total==176816
 out={'status':'PASS','complete_finite_box':full,'tested_cases':total,
      'verified_full_Bezout_identities':sum(r['verified_full_Bezout_identities']for r in rows),
      'physical_shapes':['2N+1','2N-1'],'maximum_factorial_argument':282705,
      'coefficient_construction':'Separate binomial/sliding-rising recurrences c(k+1)=c(k)*(m-k)/(k+1)*(2m+1+d*k)_(d); no factorial tables',
      'gaps':rows,'literal_integer_coefficient_controls':cc,
      'versions':{'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__},
      'seconds':time.monotonic()-start,'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
      'scope':'Only the completed ACTUAL lower-exponent-two finite box. No all-gap/all-terminal conclusion without the separately audited infinite analytical tail and complementary small-gap proofs.'}
 a.output.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','complete_finite_box':full,'tested_cases':total,'verified_full_Bezout_identities':out['verified_full_Bezout_identities'],'seconds':out['seconds']}),flush=True)
if __name__=='__main__':main()
