#!/usr/bin/env python3
"""Exact finite relative-gap box ONLY; a uniform analytic theorem is separate."""
if not __debug__:raise SystemExit('Exact verification requires unoptimized Python')
from pathlib import Path
from math import comb,factorial,isqrt
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
import argparse,json,struct,time,sys,traceback
import flint
from flint import nmod_poly
PRIMARY=1000003;ALTERNATIVE=1000033
for p in[PRIMARY,ALTERNATIVE]:
 assert all(p%j for j in range(2,isqrt(p)+1))
CASES=[(a,d,max(4000,128*(a+d))-1)for d in range(9,64)for a in range(3,d//3+1)]
assert len(CASES)==532 and sum(n-1 for a,d,n in CASES)==3694040
assert max(n for a,d,n in CASES)==10751
assert max((a+d)*n for a,d,n in CASES)==903084<PRIMARY
def tables(p,limit):
 assert limit<p
 f=[1]*(limit+1)
 for j in range(1,limit+1):f[j]=f[j-1]*j%p
 iv=[1]*(limit+1);iv[limit]=pow(f[limit],-1,p)
 for j in range(limit,0,-1):iv[j-1]=iv[j]*j%p
 return f,iv
def polynomial(m,a,d,p,tab):
 f,iv=tab;pref=f[m]*iv[a*m]%p
 c=[pref*f[a*m+d*k]*iv[k]*iv[m-k]%p for k in range(m+1)]
 assert c[0]==1 and c[-1]
 out=nmod_poly(c,p);assert out.degree()==m
 return out
def worker(case):
 a,d,last=case;p=PRIMARY;limit=(a+d)*last
 tab=tables(p,limit);other=None;previous=polynomial(1,a,d,p,tab)
 stream=sha256();fallbacks=[]
 for N in range(2,last+1):
  current=polynomial(N,a,d,p,tab)
  assert previous.degree()==N-1 and current.degree()==N
  g=previous.gcd(current);used=p
  if g.degree()!=0:
   if other is None:other=tables(ALTERNATIVE,limit)
   lo=polynomial(N-1,a,d,ALTERNATIVE,other);hi=polynomial(N,a,d,ALTERNATIVE,other)
   gg,s,t=hi.xgcd(lo)
   assert gg==nmod_poly([1],ALTERNATIVE)and s*hi+t*lo==gg,('BOTH_FIELDS_FAIL',a,d,N,g.degree(),gg.degree())
   used=ALTERNATIVE;fallbacks.append({'terminal':N,'primary_gcd_degree':g.degree(),'success_prime':used,'full_Bezout':True})
  else:assert int(g[0])==1
  stream.update(struct.pack('<IIIII',a,d,N,used,0))
  previous=current
 return{'lower_exponent':a,'gap':d,'first_terminal':2,'last_terminal':last,'tested':last-1,'max_factorial':limit,'primary_prime':p,
 'all_actual_full_degree_pairs_certified':True,'fallbacks':fallbacks,'case_stream_sha256':stream.hexdigest()}
def controls():
 for a,d in[(3,9),(7,24),(21,63)]:
  for p in[PRIMARY,ALTERNATIVE]:
   tab=tables(p,(a+d)*9)
   for n in[1,2,5,8]:
    f=polynomial(n,a,d,p,tab)
    assert list(map(int,f))==[comb(n,k)*(factorial(a*n+d*k)//factorial(a*n))%p for k in range(n+1)]
 assert nmod_poly([1,1],PRIMARY).gcd(nmod_poly([1,2,1],PRIMARY)).degree()==1
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--jobs',type=int,default=4);ap.add_argument('--output',type=Path,default=Path('/tmp/root_relative_gap_finite_box.json'));a=ap.parse_args()
 assert 1<=a.jobs<=16;controls();started=time.monotonic();rows=[]
 running=a.output.with_suffix('.running.json');journal=a.output.with_suffix('.jsonl')
 sourcehash=sha256(Path(__file__).read_bytes()).hexdigest()
 def save(complete=False):
  result={'status':'PASS'if complete else'RUNNING','scope':'ONLY the actual finite box9<=d<=63,3<=a<=floor(d/3),2<=N<max(4000,128(a+d)); a uniform all-order analytic theorem for d>=3a is a separate UNPROVED premise. No all-binomial conclusion.',
    'expected_support_pairs':532,'expected_cases':3694040,'completed_support_pairs':len(rows),'completed_cases':sum(v['tested']for v in rows),
    'max_terminal':10751,'max_factorial_argument':903084,'source_sha256':sourcehash,
    'actual_polynomials':'P_m(t)=sum binom(m,k)*(a*m+1)_(d*k)*t^k, consecutive m=N-1,N, each full degree retained',
    'primary_prime':PRIMARY,'alternative_prime':ALTERNATIVE,'rows':sorted(rows,key=lambda v:(v['gap'],v['lower_exponent'])),
    'versions':{'python':sys.version.split()[0],'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__},
    'seconds':time.monotonic()-started}
  running.write_text(json.dumps(result,indent=2)+'\n')
  if complete:a.output.write_text(json.dumps(result,indent=2)+'\n')
 save()
 with journal.open('w')as log:
  with ProcessPoolExecutor(max_workers=a.jobs)as pool:
   futures={pool.submit(worker,c):c for c in CASES}
   try:
    for future in as_completed(futures):
     row=future.result();rows.append(row);log.write(json.dumps(row)+'\n');log.flush();save()
     print('COMPLETE',len(rows),'of532 support pairs;',sum(v['tested']for v in rows),'of3694040 actual tests',flush=True)
   except BaseException as e:
    for f in futures:f.cancel()
    running.with_suffix('.failure.json').write_text(json.dumps({'status':'FAIL','case':futures.get(future),'traceback':traceback.format_exc(),'completed_supports':len(rows)},indent=2)+'\n')
    raise
 assert len(rows)==532 and sum(v['tested']for v in rows)==3694040
 save(True);print('PASS COMPLETE3694040 ACTUAL FINITE PAIRS; ANALYTIC UNIFORM TAIL IS SEPARATE',flush=True)
if __name__=='__main__':main()
