"""Reproducible exact finite controls for EisJointCosetDistribution.

Run from any directory, without Python -O. The proof is analytic; this
checks selected rational identities, every invertible matrix modulo 49,
and deliberately incorrect coefficients. It does not verify the full theorem.
No file is written unless --output is explicitly supplied.
"""
from fractions import Fraction as F
from collections import defaultdict
from math import isqrt
from pathlib import Path
import hashlib,json,sys
if not __debug__:
    raise RuntimeError('This finite diagnostic requires assertions; run without -O.')
import argparse
parser=argparse.ArgumentParser(description='Exact finite controls for the joint quaternionic coset distribution identities; not an all-target theorem verifier.')
parser.add_argument('--output',type=Path,help='Optional JSON output; otherwise print only.')
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[3]
source=ROOT/'research/arithmetic/EISENSTEIN_BOUNDED_PRIME_FACTORS.md'
source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
checker_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
NMAX=20000
r=[0]*(NMAX+1)
for d in range(1,NMAX+1):
    ch=0 if d%3==0 else (1 if d%3==1 else -1)
    if ch:
        for a in range(d,NMAX+1,d): r[a]+=6*ch

def val(a,p):
    e=0
    while a%p==0:
        a//=p;e+=1
    return e

def W(a,b):
    t=F(a-b,a+b)
    return (1-t*t)*(1+t/3)

def count(N, *,bad=(), good=(), exact=(), local=(), orient=None):
    ans=F(0)
    for a in range(1,N):
        b=N-a
        if not r[a] or not r[b]:continue
        if any(a%p and b%p for p in bad):continue
        if any(a%p==0 or b%p==0 for p in good):continue
        if any(val(a if side==0 else b,p)!=1 for p,side in exact):continue
        if orient=='left' and not (val(a,3)==1 and val(b,3)==0):continue
        if orient=='right' and not (val(a,3)==0 and val(b,3)==1):continue
        factor=F(r[a]*r[b])
        for p in local:
            va,vb=val(a,p),val(b,p)
            factor*=F(va,va+1)+F(vb,vb+1)
        ans+=factor*W(a,b)
    return ans

def ac(j):
    # Direct coefficient convolution of (1-T)^3 and sum (-T)^j.
    return sum(c*((-1)**(j-i)) for i,c in enumerate([1,-3,3,-1]) if i<=j)
def bc(e):return 2*((-1)**(e-1))

tested=[]
for p in (7,13):
  for e in range(1,5 if p==7 else 4):
    for base in (2,4):
      for controls in ({},{'exact':((19,0),)},{'local':(13 if p==7 else 7,)},{'bad':(5,)}, {'orient':'left'} if base%3==1 else {'orient':'right'}):
        N=base*p**e
        if N>NMAX:continue
        left=count(N,**controls)-count(N,good=(p,),**controls)
        right=-sum(ac(j)*count(N//p**j,**controls) for j in range(1,e+1))-bc(e)*count(base,local=tuple(controls.get('local',()))+(p,),**{k:v for k,v in controls.items() if k!='local'})
        assert left==right,(p,e,base,controls,left,right)
        tested.append((p,e,base,controls))
# Carry another singular bad test through the split recurrence.
for e in (1,2):
  for base in (26,52):
    N=base*7**e
    controls={'bad':(13,), 'orient':'left' if base%3==1 else None}
    left=count(N,**controls)-count(N,good=(7,),**controls)
    right=-sum(ac(j)*count(N//7**j,**controls) for j in range(1,e+1))-bc(e)*count(base,local=(7,),**controls)
    assert left==right,(N,left,right)
    tested.append((7,e,base,controls))

def terminal_expansion(N,ps,local=()):
    if not ps:return count(N,local=local,orient='left')
    p,*rest=ps;e=val(N,p)
    if p%3==2:
        return F(0) if e==1 else terminal_expansion(N//(p*p),rest,local)
    ans=F(0)
    for j in range(1,e+1):ans-=ac(j)*terminal_expansion(N//p**j,rest,local)
    ans-=bc(e)*terminal_expansion(N//p**e,rest,local+(p,))
    return ans
multi=[]
for N,ps in ((4*7*13,(7,13)),(4*7**2*13,(7,13)),(4*7*13**2,(13,7)),(4*7*5,(7,5)),(4*7*5**2,(5,7))):
    actual=count(N,bad=ps,orient='left')
    expanded=terminal_expansion(N,list(ps))
    assert actual==expanded,(N,ps,actual,expanded)
    multi.append((N,ps,str(actual)))

scalar_checks=0
for p in (7,13,19,31):
  sig=lambda e:F(p**(e+1)-1,p-1)
  for e in range(1,61):
    bad=-sum(ac(j)*sig(e-j) for j in range(1,e+1))-bc(e)*F(2,p+1)
    expected=sig(e)-F((p-1)**2*p**(e-1),p+1)
    assert bad==expected
    scalar_checks+=1

# Exhaust the actual prime-square local indicator at p=7, on every invertible matrix.
p=7;q=p*p
squares={a*a%p for a in range(1,p)}
total=[0,0];events=[0,0];formula_checks=0
for a in range(q):
  E1=int(a%p==0);E2=int(a==0)
  for d in range(q):
    F1=int(d%p==0);F2=int(d==0)
    test=int((a*d)%p==0 and (a*d)%q!=0)
    formula=E1+F1-E2-F2-2*E1*F1+E2*F1+E1*F2
    assert test==formula,(a,d)
    for b in range(q):
      for c in range(q):
        det=(a*d-b*c)%p
        if det==0:continue
        cl=int(det not in squares)
        total[cl]+=1;events[cl]+=test;formula_checks+=1
assert F(events[0],total[0])==F(events[1],total[1])==F(2*(p-1)**2,p*p*(p+1))
negative_actual=count(364,bad=(7,13),orient='left')
negative_perturbed=3*count(52,bad=(13,),orient='left')-2*count(52,bad=(13,),local=(7,),orient='left')
assert negative_actual!=negative_perturbed
control_a=control_d=7
control_b=control_c=1
control_E1=int(control_a%7==0); control_F1=int(control_d%7==0)
control_E2=int(control_a%49==0); control_F2=int(control_d%49==0)
control_wrong=control_E1+control_F1-control_E2-control_F2+control_E2*control_F1+control_E1*control_F2
control_truth=int(control_a*control_d%7==0 and control_a*control_d%49!=0)
assert (control_a*control_d-control_b*control_c)%7!=0
assert control_wrong!=control_truth
assert (r[1],r[3],r[7],r[5],r[25])==(6,6,12,0,6)
result={'coverage':'Finite exact recurrence, local density and negative controls only; not the spectral theorem or all-target estimate.', 'checker_sha256':checker_hash,'source_sha256':source_hash,'python':sys.version.split()[0], 'single_split_recurrence_cases':len(tested),'asymmetric_weight':'(1-t^2)(1+t/3), exact rational arithmetic, endpoints omitted as zero','multiple_singular_cases':multi,'scalar_density_checks':scalar_checks,'scalar_max_exponent':60,'prime_square':{'p':p,'invertible_matrices':formula_checks,'class_totals':total,'class_events':events,'exact_density':str(F(events[0],total[0]))},'negative_controls':{'perturbed_a1_rejected':True,'perturbed_C6_rejected':True},'outcome':'PASS'}
if args.output is not None:
    args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
