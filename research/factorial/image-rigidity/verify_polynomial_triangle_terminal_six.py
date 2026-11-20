#!/usr/bin/env python3
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions required')
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,isqrt
from hashlib import sha256
import argparse,json,struct,copy
HERE=Path(__file__).resolve().parent
PROOF='POLYNOMIAL_TRIANGLE_AND_TERMINAL_SIX.md';INPUTS='polynomial_triangle_terminal_six_inputs.json';LEDGER='polynomial_triangle_terminal_six_source_ledger.json';RECEIPT='polynomial_triangle_terminal_six_verification.json'
def digest(p):return sha256(p.read_bytes()).hexdigest()
def cutoff(N):
 assert N>=2;e=(4*N*(N-1)-1).bit_length();return e,2*(N-1)*e

def domain(inp):
 assert set(inp)=={'terminal','prime','binary_exponent','gap_cutoff','strict_lower_factor','expected_cases','maximum_factorial_argument','case_stream_sha256'}
 N,p=inp['terminal'],inp['prime'];assert N==6 and p==1000003
 assert p>1 and all(p%j for j in range(2,isqrt(p)+1))
 e,D=cutoff(N);assert(e,D,2*N-3)==(inp['binary_exponent'],inp['gap_cutoff'],inp['strict_lower_factor'])==(7,70,9)
 pairs=[(a,d)for d in range(2,D)for a in range(1,d*(2*N-3))]
 assert len(pairs)==inp['expected_cases']==21658
 assert max((a+d)*N for a,d in pairs)==inp['maximum_factorial_argument']==4134<p
 return N,p,pairs

def tables(p,L):
 f=[1]*(L+1)
 for j in range(1,L+1):f[j]=f[j-1]*j%p
 g=[1]*(L+1);g[L]=pow(f[L],-1,p)
 for j in range(L,0,-1):g[j-1]=g[j]*j%p
 return f,g

def polynomial(a,d,n,p,tab):
 f,g=tab;v=[comb(n,k)*f[a*n+d*k]*g[a*n]%p for k in range(n+1)]
 assert len(v)==n+1 and v[0]==1 and v[-1]!=0;return v

def trim(v):
 while v and v[-1]==0:v.pop()
 return v

def gcd(u,v,p):
 u,v=u[:],v[:]
 while v:
  q=u[:];inv=pow(v[-1],-1,p)
  while len(q)>=len(v):
   c=q[-1]*inv%p;off=len(q)-len(v)
   for j,b in enumerate(v):q[off+j]=(q[off+j]-c*b)%p
   trim(q)
  u,v=v,q
 inv=pow(u[-1],-1,p);return [c*inv%p for c in u]

def ev(v,t):
 out=Q(0)
 for c in reversed(v):out=out*t+c
 return out

def root_path(arg):
 if arg:return arg.resolve()
 for p in[HERE,*HERE.parents]:
  if(p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Use --repo-root for relocated replay')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--write-receipt',action='store_true');args=ap.parse_args();root=root_path(args.repo_root)
 source=json.loads((HERE/LEDGER).read_text());mandatory={('packet',PROOF),('packet',Path(__file__).name),('packet',INPUTS)}|{('repository','research/factorial/image-rigidity/'+n)for n in['lower_exponent_initial_and_boundary.tex','GAMMA_PRODUCT_ROOT_MODULUS.tex','ALL_UNIVARIATE_TWO_TERM_SFC.tex','UNIFORM_PADE_PARABOLA_SFC.tex']}
 keys={(v['location'],v['path'])for v in source['files']};assert mandatory<=keys and len(keys)==len(source['files']);pins=[]
 for item in source['files']:
  p=(root if item['location']=='repository'else HERE)/item['path'];h=digest(p);assert h==item['sha256'],str(p);pins.append({'path':item['path'],'sha256':h})
 inp=json.loads((HERE/INPUTS).read_text());N,p,pairs=domain(inp);tab=tables(p,4134)
 literal=0
 for a in[1,2,7,20]:
  for d in[2,9,69]:
   for m in[5,6]:
    v=polynomial(a,d,m,p,tab);w=[comb(m,k)*(factorial(a*m+d*k)//factorial(a*m))%p for k in range(m+1)];assert v==w;literal+=1
 rational=0
 for a,d in[(0,1),(1,2),(2,3),(7,4)]:
  for m in[5,6]:
   A,B=Q(-2,3),Q(5,7);v=[comb(m,k)*(factorial(a*m+d*k)//factorial(a*m))for k in range(m+1)]
   actual=sum(comb(m,k)*A**(m-k)*B**k*factorial(a*m+d*k)for k in range(m+1));assert actual==A**m*factorial(a*m)*ev(v,B/A);rational+=1
 stream=sha256()
 for a,d in pairs:
  u,v=polynomial(a,d,5,p,tab),polynomial(a,d,6,p,tab);g=gcd(u,v,p);assert g==[1],(a,d,g)
  stream.update(struct.pack('<IIII',a,d,p,0))
 assert stream.hexdigest()==inp['case_stream_sha256']
 # Exact integer-cutoff controls, not moment scans at later terminals.
 for n in range(2,129):
  e,D=cutoff(n);assert 2**(e-1)<4*n*(n-1)<=2**e
  assert D==2*(n-1)*e
  count=sum(d*(2*n-3)-1 for d in range(2,D));assert count==(2*n-3)*(D*(D-1)//2-1)-(D-2)
 # Malformed domains and an actual full-degree common-root corruption fail.
 for key,value in [('terminal',7),('gap_cutoff',69),('prime',1000005)]:
  bad=copy.deepcopy(inp);bad[key]=value
  try:domain(bad)
  except AssertionError:pass
  else:raise AssertionError('Malformed domain was accepted')
 u,v=polynomial(1,2,5,p,tab),polynomial(1,2,6,p,tab);u[0]=-sum(u[1:])%p;v[0]=-sum(v[1:])%p
 assert len(u)==6 and len(v)==7 and len(gcd(u,v,p))>1
 wrong=[comb(6,k)*(factorial(5+2*k)//factorial(5))%p for k in range(7)]
 actual=polynomial(1,2,6,p,tab)
 assert wrong!=actual and wrong[1]==252 and actual[1]==336 # Incorrect previous shared Gamma6 versus actual Gamma7.
 result={'status':'PASS','scope':'Complete ACTUAL terminal-six finite triangle only; the all-exponent theorem and all-N>=7 equivalence require the written complement arguments. No all-order Strong Factorial closure.','terminal':6,'tested_cases':len(pairs),'gap_cutoff':70,'binary_exponent':7,'prime':p,'maximum_factorial_argument':4134,'actual_shapes':['5*a+1','6*a+1'],'all_full_degree_gcds_equal_one':True,'case_stream_sha256':stream.hexdigest(),'literal_factorial_polynomial_controls':literal,'rational_moment_controls':rational,'integer_cutoff_controls':127,'negative_controls':['three_malformed_domains','full_degree_common_root_corruption','incorrect_shared_Gamma_shape'],'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__).resolve()),'inputs_sha256':digest(HERE/INPUTS),'ledger_sha256':digest(HERE/LEDGER),'verified_source_pins':pins}
 out=args.output or HERE/RECEIPT;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print('PASS21658 actual degree5/6 gcds;',literal,'literal polynomial controls;',rational,'rational moment controls;127integer cutoffs');print(out)
if __name__=='__main__':main()
