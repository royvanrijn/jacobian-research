"""A bounded exact norm elimination of w on the delta-zero branch.
A constant-leading cubic gives a finite free rank3 w presentation globally.
No coefficient denominator or exceptional chart is discarded.
"""
from pathlib import Path
from math import gcd
from functools import reduce
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def worker():
 from flint import fmpz_mpoly_ctx
 src=json.loads((HERE/'delta-zero-fivevar-input.json').read_text());ctx=fmpz_mpoly_ctx.get(src['variables']);r,u,v,z,w=ctx.gens();eq={row['old_coefficient']:ctx.from_dict({tuple(ex):int(c) for ex,c in row['terms']}) for row in src['equations']}
 def wc(F,n):return ctx.from_dict({ex[:4]+(0,):int(c) for ex,c in F.to_dict().items() if ex[4]==n})
 cubic=eq[3];leading=wc(cubic,3);assert leading==3072
 rows=[];quadratics=[]
 for index in [1,2,4,5]:
  F=eq[index];n=int(F.degrees()[4]);assert n in [2,3]
  if n==3:F=3072*F-wc(F,3)*cubic
  content=reduce(gcd,(abs(int(c)) for c in F.coeffs()),0);F=F//content;assert F.degrees()[4]<=2
  quadratics.append(F);emit({'stage':'quadratic','coefficient':index,'terms':len(F.to_dict()),'degrees':list(F.degrees()),'seconds':time.monotonic()-START})
  N=cubic.resultant(F,'w');td=N.to_dict();assert all(ex[4]==0 for ex in td);content=reduce(gcd,(abs(int(c)) for c in N.coeffs()),0);N=N//content
  rows.append({'coefficient':index,'quadratic_terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(F.to_dict().items())],'norm_content':str(content),'norm_terms':[[list(map(int,ex[:4])),str(c)] for ex,c in sorted(N.to_dict().items())],'degrees':list(map(int,N.degrees()[:4])),'total_degree':int(max(sum(ex) for ex in N.to_dict()))})
  emit({'stage':'norm','coefficient':index,'terms':len(N.to_dict()),'degrees':list(N.degrees()),'seconds':time.monotonic()-START})
  (HERE/'delta-zero-norm-checkpoint.json').write_text(json.dumps({'status':'UNVERIFIED_NORM_CANDIDATE','norms':rows},separators=(',',':'))+'\n')
 out={'status':'UNVERIFIED_NORM_CANDIDATE','variables':['r','u','v','z'],'source_sha256':hashlib.sha256((HERE/'delta-zero-fivevar-input.json').read_bytes()).hexdigest(),'constant_cubic_leading_coefficient':3072,'norms':rows,'pure_equation_terms':[[list(map(int,ex[:4])),str(c)] for ex,c in sorted(eq[7].to_dict().items())],'boundary':'All norms are necessary without division, since the cubic leading coefficient is a nonzero constant. Individual norms may have common zeros using different w roots. No characteristic-zero exclusion follows without replayed identities or an exact integer witness.'}
 (HERE/'delta-zero-norm-input.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
if __name__=='__main__':
 if '--worker' in sys.argv:worker()
 else:
  def limits():
   resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(50*1024**2,50*1024**2))
  receipt={'cap_seconds':60,'address_cap_GiB':4,'per_file_cap_MiB':50,'theorem_exclusion':False,'motivation':'Eliminate one variable through a constant-leading cubic and four quadratic norms. This is a new four-variable certificate target, not a restart of the capped five-variable rational lift.'}
  with (HERE/'delta-zero-norm-build.log').open('wb') as fd:
   try:
    proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limits);receipt.update(outcome='completed',returncode=proc.returncode)
   except subprocess.TimeoutExpired:receipt.update(outcome='timeout',returncode=None)
  receipt['seconds']=time.monotonic()-START
  for fn in ['delta-zero-norm-input.json','delta-zero-norm-checkpoint.json']:
   p=HERE/fn
   if p.exists():receipt[fn]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
  (HERE/'delta-zero-norm-build-run.json').write_text(json.dumps(receipt,indent=2)+'\n');emit(receipt);print((HERE/'delta-zero-norm-build.log').read_text()[-2500:])
