#!/usr/bin/env sage-python
"""Uniform injectivity of the class1 smooth-fibre branch-divisor map."""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,lcm,gcd,prime_range
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
EQ=ART/'x1092_class1_realization_equation_v1.json'
CARRIERS=ART/'curve302_class1_fibre_carriers_v1.json'
OUTPUT=ART/'curve302_class1_branch_veronese_v1.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compute(check=False):
    eq,data=map(read,[EQ,CARRIERS]);R=PolynomialRing(QQ,'u');K=R.fraction_field()
    cs=[K(R(c['numerator']))/R(c['denominator'])for c in eq['quartic_coefficients']]
    denominator=lcm([c.denominator()for c in cs]);raw=[R(denominator*c)for c in cs]
    common=gcd(raw);ps=[q//common for q in raw]
    assert len(ps)==5 and all(p.degree()==4 for p in ps)
    M=matrix(QQ,[p.list()for p in ps]);assert M.dimensions()==(5,5)
    witness=None
    for pp in prime_range(3,2001):
        p=int(pp)
        if any(c.denominator()%p==0 for c in M.list()):continue
        value=int(M.change_ring(GF(p)).det())
        if value:witness={'prime':p,'determinant_mod_prime':value};break
    assert witness is not None
    for row in data['carriers']:
        u=QQ(row['u']);v=[p(u)for p in ps];q=list(map(QQ,row['branch_quartic']))
        k=next(i for i in range(5)if v[i]);ratio=q[k]/v[k]
        assert ratio and q==[ratio*x for x in v]
    if check:
        assert ZZ(witness['prime']).is_prime(proof=True)
        # Exact rational invertibility is independent of modular determinant.
        inverse=M.inverse();assert inverse*M==matrix.identity(QQ,5)
        assert all(denominator*c==common*p for c,p in zip(cs,ps))
    result={'schema':'curve302.class1-branch-veronese.v1',
            'bindings':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),EQ,CARRIERS]},
            'common_denominator':list(map(str,denominator.list())),'removed_polynomial_content':list(map(str,common.list())),
            'coefficient_matrix':[[str(x)for x in row]for row in M],'finite_certificate':witness,
            'matrix_rank':5,'projective_map_degree':4,'parameter_map':'M*(U^4,U^3*V,U^2*V^2,U*V^3,V^4)^t, where u=V/U',
            'status':'PASS_VERONESE_BRANCH_EMBEDDING',
            'theorem':'Different parameters have different projective binary quartics. Thus distinct smooth class1 fibres never define the same quadratic extension of the A-base, even up to a constant twist; every joint normalization of two such smooth covers has genus at least2.',
            'boundary':'Smooth fibres in this explicit family only. Singular-fibre normalizations and other multisections are outside this theorem.',
            'limits':{'wall_seconds':60,'workers':1,'coefficient_matrix_size':5,'new_parameters':0,'searches':0}}
    print('PASS branch map is a degree4 Veronese embedding; prime',witness['prime'],flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    result=compute(args.check)
    if args.check:assert read(OUTPUT)==result
    else:
        with OUTPUT.open('x')as out:json.dump(result,out,indent=2,sort_keys=True);out.write('\n')
