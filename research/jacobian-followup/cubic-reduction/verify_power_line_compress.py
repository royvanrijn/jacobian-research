#!/usr/bin/env python3
"""Independent Fraction/gcd certificate against pair-line single-cube compression."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict
import hashlib,json,time
from verify import add,mul,variable
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean(p):return {k:v for k,v in p.items() if v}
def pscale(p,c):return clean({e:c*v for e,v in p.items()})
def lin(f):
 out={}
 for i,c in enumerate(f):
  if c:out=add(out,variable(i,18),Q(c))
 return out
# Dense univariate polynomials in ascending order; no CAS or floating point.
def trim(a):
 a=list(a)
 while a and a[-1]==0:a.pop()
 return a
def rem(a,b):
 a=trim(a);b=trim(b);assert b
 while len(a)>=len(b):
  q=a[-1]/b[-1];k=len(a)-len(b)
  for i,c in enumerate(b):a[k+i]-=q*c
  a=trim(a)
 return a
def pgcd(a,b):
 a=trim(a);b=trim(b)
 while b:a,b=b,rem(a,b)
 return [c/a[-1] for c in a] if a else []
def derivative(a):return trim([Q(i)*a[i] for i in range(1,len(a))])
def quotient(a,b):
 a=trim(a);b=trim(b);assert b;out=[Q(0)]*max(0,len(a)-len(b)+1)
 while len(a)>=len(b):
  q=a[-1]/b[-1];k=len(a)-len(b);out[k]+=q
  for i,c in enumerate(b):a[k+i]-=q*c
  a=trim(a)
 assert not a
 return trim(out)
def constraint(a,b):
 ratio=None
 for e in set(a)|set(b):
  x=a.get(e,Q(0));y=b.get(e,Q(0))
  if not x or not y:return False,None
  value=-x/y
  if ratio is not None and value!=ratio:return False,None
  ratio=value
 return True,ratio

def main():
 source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());n=d["size"];assert n==89
 assert d["source_sha256"]==sha(HERE/"power_basis_anneal.py")
 B=[[Q(c) for c in row] for row in d["B"]];lins=[lin(f) for f in d["forms"]];squares=[mul(l,l) for l in lins];cubes=[mul(s,l) for s,l in zip(squares,lins)]
 basis={}
 for j,cube in enumerate(cubes):
  p=dict(cube);coords={j:Q(1)}
  for e,(b,bc) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,b,-c);coords=add(coords,bc,-c)
  assert p,"The saved89 cubes must be linearly independent"
  pivot=min(p);c=p[pivot];basis[pivot]=(pscale(p,1/c),pscale(coords,1/c))
 def project(p):
  p=dict(p);coords={}
  for e,(b,bc) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,b,-c);coords=add(coords,bc,c)
  return p,coords
 for j,cube in enumerate(cubes):
  residual,coords=project(cube);assert not residual and coords=={j:Q(1)}
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);scale=next(c for c in v if c)
  groups[tuple(c/scale for c in v)].append(j)
 eligible=[(j,k) for key,ids in groups.items() if key[17]==0 for j,k in combinations(ids,2)]
 support=set().union(*(set(p) for p in cubes));started=time.monotonic()
 stats=dict(pairs=0,outside_rejected=0,residual_rejected=0,fixed_points=0,full_binary_pairs=0,eligible_output_pairs=len(eligible),gcd_exclusions=0)
 fixed_certificates=[];families=[]
 for p,q in combinations(range(n),2):
  stats["pairs"]+=1;E1=pscale(mul(squares[p],lins[q]),Q(3));E2=pscale(mul(lins[p],squares[q]),Q(3))
  a={e:c for e,c in E1.items() if e not in support};b={e:c for e,c in E2.items() if e not in support}
  ok,t0=constraint(a,b)
  if not ok:stats["outside_rejected"]+=1;continue
  r1,c1=project(E1);r2,c2=project(E2);ok,t1=constraint(r1,r2)
  if not ok or (t0 is not None and t1 is not None and t0!=t1):stats["residual_rejected"]+=1;continue
  fixed=t0 if t0 is not None else t1
  coefficients=[[Q(i==p),c1.get(i,Q(0)),c2.get(i,Q(0)),Q(i==q)] for i in range(n)]
  if fixed is not None:
   stats["fixed_points"]+=1
   values=[sum((c*fixed**a for a,c in enumerate(row)),Q(0)) for row in coefficients]
   # Cancellation needs two proportional coefficient columns with exactly these coordinate ratios.
   for j,k in eligible:
    v=next(o for o in range(18) if B[o][j])
    assert not values[j] or B[v][k]*values[j]!=B[v][j]*values[k],(p,q,j,k,fixed)
   fixed_certificates.append(dict(pair=[p,q],parameter=str(fixed)))
  else:
   stats["full_binary_pairs"]+=1
   for j,k in eligible:
    aj=trim(coefficients[j]);ak=trim(coefficients[k])
    if not aj:continue
    v=next(o for o in range(18) if B[o][j])
    ratio=trim([B[v][k]*a-B[v][j]*b for a,b in zip(coefficients[j],coefficients[k])])
    assert ratio,"An identically satisfied ratio with aj nonzero would permit compression"
    # All complex roots of ratio must be zeros of aj (or the excluded parameter0).
    radical=quotient(ratio,pgcd(ratio,derivative(ratio)))
    forbidden=[Q(0)]+aj
    assert len(pgcd(radical,forbidden))==len(radical),(p,q,j,k,ratio,aj)
    stats["gcd_exclusions"]+=1
   families.append([p,q])
 receipt=dict(status="pass",source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/"verify.py"),
  input_sha256=sha(source),primary_source_sha256=sha(HERE/"power_line_compress.py"),
  primary_checkpoint_sha256=sha(HERE/"power-line-compress/checkpoint.json"),
  **stats,seconds=round(time.monotonic()-started,3),fixed_certificates=fixed_certificates,full_binary_pairs_list=families,
  independent_cube_rank=89,all_pairs_covered=True,
  scope="No two or more saved89 cubes can be removed and replaced by one cube whose linear form lies in the source span of any two saved forms, while keeping all other selected forms and re-solving their coefficients. Characteristic-zero extensions included. Does not exclude two new cubes or forms outside these pair spans, nor all smaller decompositions.")
 primary=json.loads((HERE/"power-line-compress/checkpoint.json").read_text())
 assert primary["status"]=="completed-no-compression" and primary["source_sha256"]==sha(HERE/"power_line_compress.py")
 for key in ["pairs","outside_rejected","residual_rejected","fixed_points"]:assert primary[key]==stats[key]
 (HERE/"power-line-compress/validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
 print(json.dumps({k:v for k,v in receipt.items() if k not in ["fixed_certificates","full_binary_pairs_list"]}))
if __name__=="__main__":main()
