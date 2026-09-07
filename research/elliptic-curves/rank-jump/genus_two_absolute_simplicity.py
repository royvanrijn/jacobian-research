#!/usr/bin/env python3
"""One Frobenius certificate excluding every geometric elliptic factor."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,cyclotomic_polynomial,euler_phi
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'GENUS_TWO_RANK_AND_SIMPLICITY_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_native_genus_two_lift_gate_v1.json'
OUTPUT=r.OUT/'rank_jump_genus_two_absolute_simplicity_v1.json'


def independent_count(coeffs,p,d):
    """Integer-pair arithmetic in Fp[w]/(w^2-d), using the norm character."""
    count=2 # nonzero rational leading coefficient is square over Fp^2
    for a in range(p):
        for b in range(p):
            u=v=0
            for c in reversed(coeffs):u,v=(u*a+d*v*b+c)%p,(u*b+v*a)%p
            norm=(u*u-d*v*v)%p
            chi=pow(norm,(p-1)//2,p)
            count+=1+(0 if chi==0 else 1 if chi==1 else -1)
    return count


def compute():
    inp=r.read(INPUT);p=131
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'x');pol=R(inp['sextic_coefficients']);F=GF(p);rp=PolynomialRing(F,'x');q=rp(pol)
    assert q.degree()==6 and q.is_squarefree()
    n1=sum(1 if q(a)==0 else 2 if q(a).is_square() else 0 for a in F)
    n1+=2 if q.leading_coefficient().is_square() else 0
    assert n1==next(x['genus_two_points'] for x in inp['finite_field_lift_checks'] if x['prime']==p)
    d=next(k for k in range(2,p) if not F(k).is_square())
    K=GF(p*p,'w',modulus=rp.gen()**2-F(d));rpk=PolynomialRing(K,'x');qk=rpk(q)
    n2=sum(1 if qk(a)==0 else 2 if qk(a).is_square() else 0 for a in K)+2
    assert n2==independent_count([int(x) for x in q.list()],p,d)
    a1=p+1-n1;s2=p*p+1-n2;a2=(a1*a1-s2)//2
    assert 2*a2==a1*a1-s2
    x=R.gen();P=x**4-a1*x**3+a2*x**2-p*a1*x+p*p
    Z=PolynomialRing(QQ,'z');z=Z.gen();T=PolynomialRing(Z,'x');xx=T.gen()
    ratio=Z(T(P.list()).resultant(sum(P[i]*(z*xx)**i for i in range(5)))).monic()
    diagonal=(z-1)**4;off,remainder=ratio.quo_rem(diagonal)
    assert remainder==0 and off.degree()==12 and off(1)!=0
    orders=[n for n in range(2,289) if euler_phi(n)<=12]
    cyclo_hits=[n for n in orders if off.gcd(Z(cyclotomic_polynomial(n))).degree()>0]
    passed=P.is_irreducible() and not cyclo_hits
    return {'schema':'rank-jump.genus-two-absolute-simplicity.v1','status':'PASS' if passed else 'UNKNOWN',
        'prime':p,'quadratic_nonresidue':d,'points_Fp':int(n1),'points_Fp2':int(n2),
        'independent_integer_pair_count_agrees':True,'Frobenius_polynomial':list(map(str,P.list())),
        'Frobenius_polynomial_irreducible':bool(P.is_irreducible()),
        'off_diagonal_eigenvalue_ratio_polynomial':list(map(str,off.list())),
        'cyclotomic_orders_checked':orders,'cyclotomic_factors_found':cyclo_hits,
        'absolute_simplicity':bool(passed),'geometric_elliptic_quotient_exists':False if passed else 'UNKNOWN',
        'bindings':{str(path.relative_to(r.ROOT)):r.digest(path.read_bytes()) for path in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':'One exact good-reduction certificate and the Frobenius/conjugate proof in the note. This proves geometric simplicity, not the rational rank, Selmer dimension, or rational points of the genus-two curve.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print(d['status'],'Frobenius',d['Frobenius_polynomial'],'cyclotomic hits',d['cyclotomic_factors_found'])
