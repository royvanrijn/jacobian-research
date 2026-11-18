#!/usr/bin/env python3
"""Portable exact finite local controls for critical semilinear switching.

Run without Python -O. Output is printed unless --output is supplied.
The global distribution estimate BVcrit is not verified by these controls.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,sys
if not __debug__:
    raise RuntimeError("Exact controls require assertions; run without Python -O.")
parser=argparse.ArgumentParser(description="Exact finite switching local controls; not a sieve distribution verifier.")
parser.add_argument("--output",type=Path,help="Optional JSON receipt; otherwise print only.")
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[3]
def A(e):return (e+1)**2
def B(e):return e+1
def H(p,v=0):
    x=F(1,p);q=1-x
    return 1/q+2*x/((v+1)*q*q)+x*(1+x)/((v+1)**2*q**3)
def K(p,r,j):
    x=F(1,p);q=1-x
    if r!=j:return q*A(min(r,j))
    tail=(1+x)/q**3-sum(A(i)*x**i for i in range(r+1))
    return (1-2*x)*A(r)+q*x**(-r)*tail

def val(n,p):
    e=0
    while n and n%p==0:e+=1;n//=p
    return e

def direct_K(p,r,j):
    shift=min(r,j);h=2;mod=p**h;x=F(1,p);q=1-x
    average0=q*x**(-h)*(H(p,shift)-sum(F(A(shift+i),A(shift))*x**i for i in range(h)))
    answer=F(0)
    for s in range(1,mod):
        if s%p==0:continue
        value=(p**(r-shift)-p**(j-shift)*s)%mod
        answer+=F(A(shift+val(value,p)),A(shift)) if value else average0
    return A(shift)*answer/mod

def integral_I(p,r):
    x=F(1,p);q=1-x;ha=(1+x)/q**3;hb=1/q**2
    ta=ha-sum(A(i)*x**i for i in range(r+1));tb=hb-sum(B(i)*x**i for i in range(r+1))
    return q*sum(A(i)*B(i)*x**i for i in range(r))+x**r*(1-2*x)*A(r)*B(r)+q*(B(r)*ta+A(r)*tb)

def euler_K(p,r):
    x=F(1,p);q=1-x
    prefix=sum(B(j)*x**j*K(p,r,j) for j in range(r+1))
    tail=1/q**2-sum(B(j)*x**j for j in range(r+1))
    return prefix+q*A(r)*tail
nd=ne=nh=0
for p in [13,37,61]:
    x=F(1,p);q=1-x
    for r in range(13):
        for j in range(13):
            direct=direct_K(p,r,j);formula=K(p,r,j);assert direct==formula
            # Primitive normalized forms have 2 roots if r=j, else one root.
            c00=1-(2 if r==j else 1)*x
            actual=direct/A(min(r,j));T=actual-c00
            assert T>=0 and c00*(1+T)<=actual
            nd+=1;nh+=1
    for r in range(51):
        assert euler_K(p,r)==integral_I(p,r)
        assert integral_I(p,r)>=integral_I(p,0)
        ne+=1
    assert K(p,0,0)>=q
# Deliberate mistakes lose an exact content factor or treat E as the expectation.
p=13;r=j=1;lost_content=K(p,r,j)/A(1)
assert lost_content!=K(p,r,j)
x=F(1,p);actual=K(p,0,0);c00=1-2*x;E=1+actual-c00
assert E>actual and c00*E<=actual
paths=['research/arithmetic/ASYMMETRIC_CRITICAL_SWITCHING.md','research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md','research/arithmetic/ASYMMETRIC_SEMILINEAR_SIEVE.md','research/arithmetic/sources/henriot-erratum-2014.pdf']
result={'scope':'Exact local K, content and corrected-Henriot inequality diagnostics; BVcrit remains unproved','direct_K_cases':nd,'henriot_local_inequality_cases':nh,'all_valuation_euler_cases':ne,'content_negative_control':{'p':13,'r':1,'j':1,'true_K':str(K(13,1,1)),'missing_Uh_K':str(lost_content)},'expectation_negative_control':{'actual':str(actual),'uncorrected_interpretation_E':str(E)},'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}
result['checker']={'path':Path(__file__).relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version.split()[0]}
result['result']='PASS'
serialized=json.dumps(result,indent=2)+'\n'
if args.output is not None:
    args.output.write_text(serialized)
print(serialized,end='')
