"""Finite normalization checks for the written all-order counterexample proof.

These checks do not replace the symbolic Haar/Cartier/Dirichlet argument.
The moment calculation uses the ORIGINAL polynomial, independently of the
exploratory reciprocal-transvectant calculation.
"""
from pathlib import Path
from math import comb, factorial
from fractions import Fraction
import hashlib, json, time
import sympy as S

ROOT = Path(__file__).resolve().parents[3]
x, y, z, w, t = S.symbols('x y z w t')

def partitions(total, limit=None):
    if total == 0:
        yield ()
        return
    for q in range(min(total, limit or total), 1-1, -1):
        for rest in partitions(total-q, q):
            yield (q,)+rest

def check(n):
    begin=time.monotonic()
    F = x**8+14*x**4*y**4+y**8 if n==8 else x*y*(x**10+11*x**5*y**5-y**10)
    G = x**2*y**2*(x**4-y**4)**2 if n==8 else x**20+y**20-228*x**15*y**5+494*x**10*y**10+228*x**5*y**15
    v=S.Poly(G,x,y).total_degree()
    C=S.expand((S.diff(F,x)*S.diff(G,y)-S.diff(F,y)*S.diff(G,x))/(n*v))
    relation=C**2-F**3*G/36+3*G**3 if n==8 else C**2-12*F**5-G**3/144
    assert S.Poly(S.expand(relation),x,y).is_zero
    assert S.expand(F.subs({x:-y,y:x}, simultaneous=True)-F)==0
    # At y=w=1, bidegrees recover every discarded exponent uniquely.
    m=2; b=n//2; a=b-1; A=a*m; B=b*m; s=b*(m-1)
    j=(n-1)*B+b; p=b*(n-1)*m-a; N=p-1
    assert S.isprime(p)
    H=S.Poly(F.subs({x:z+t*x,y:1+t}, simultaneous=True),t)
    Hq=[S.Poly(H.nth(q),x,z) for q in range(s+1)]
    Fb=S.Poly(F.subs({x:z,y:1}, simultaneous=True),x,z)
    U=S.Poly(F.subs(y,1),x,z)
    J=S.Poly(0,x,z)
    # [t^s] F(bottom+t*top)^B / F(bottom), with no rational division.
    for part in partitions(s):
        count={q:part.count(q) for q in set(part)}; k=len(part)
        coefficient=factorial(B)//factorial(B-k)
        for multiplicity in count.values(): coefficient//=factorial(multiplicity)
        summand=coefficient*Fb**(B-k-1)
        for q,multiplicity in count.items():summand*=Hq[q]**multiplicity
        J+=summand
    R=U**A*J
    total=S.Integer(0)
    for (ix,iz),coefficient in R.terms():
        if ix+iz==N:
            total+=coefficient*(-1)**iz*factorial(ix)*factorial(iz)
    moment=S.cancel(total/factorial(N+1))
    expected={8:S.Rational(-865594153113747456,298435517190385025),12:S.Rational(41958690403313837610721588134765625,1532800285927159272199695730537217233992065276406)}
    assert moment==expected[n]
    # The general rank-one congruence is checked on this original polynomial.
    scaled=S.cancel(p*moment)
    residue=int(S.numer(scaled)%p)*pow(int(S.denom(scaled)%p),-1,p)%p
    h=(p-1)//b; c=14 if n==8 else 11; sign=1 if n==8 else -1
    Hasse=sum(comb(h,i)*comb(h-i,i)*pow(c,h-2*i,p)*pow(sign,i,p) for i in range(h//2+1))%p
    k=(p-1)//2
    elliptic=(comb(k,k//2)*pow(-pow(12,-1,p),k//2,p))%p if n==8 else (comb(k,(p-1)//3)*pow(12,k-(p-1)//3,p)*pow(pow(144,-1,p),(p-1)//3,p))%p
    assert Hasse==elliptic!=0
    assert residue==comb(n*B,j)*Hasse%p!=0
    return dict(n=n,m=m,p=p,moment=str(moment),p_times_moment_mod_p=residue,Hasse=Hasse,elliptic=elliptic,original_zero_weight_terms=len(R.terms()),Klein_relation='PASS',self_duality='PASS',independent_formula_agreement='PASS',seconds=time.monotonic()-begin)

if __name__=='__main__':
    results=[check(8),check(12)]
    output={'checks':results,'scope':'Finite exact normalization checks only. Infinite fixed-marker tail is proved in the canonical written argument, not by these samples.','external_review':False,'formal_verification':False}
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    if args.output:args.output.write_text(json.dumps(output,indent=2)+'\n')
    for row in results:print(json.dumps(row),flush=True)
