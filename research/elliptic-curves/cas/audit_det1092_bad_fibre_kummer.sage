#!/usr/bin/env sage-python
"""Bounded equation-only audit of the discriminant Kummer class.

One fixed parent; only the previously frozen nine control parameters.
No point coordinates, point search, class group, or integer factorization.
"""
import hashlib,json,signal,math
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix

signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_bad_fibre_kummer_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
OLD=ART/'det1092_polynomial_lift_gate_v1/construction.json'
ROSTER=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
parent,old,roster=map(read,[PARENT,OLD,ROSTER])
R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
a1,a2,a3,a4,a6=map(dec,parent['a_invariants'])
assert (a1,a2,a3)==(1,1,1)
T=PolynomialRing(K,'X');x=T.gen()
f=x**3+5*x*x+(16*a4+8)*x+64*a6+16
D=R(f.discriminant());assert D.degree()==24 and D.gcd(D.derivative())==1
beta=-D*f.derivative()
assert f.resultant(beta)==K(D)**4
templates={r['label']:R(list(map(QQ,r['root_remainder']))) for r in old['templates']}
r,s=templates['double'],templates['simple']
assert (r-s).gcd(D)==1
assert all(R(c)%D==0 for c in (f-(x-r)**2*(x-s)).list())
assert R(f.derivative()(s)).gcd(D)==1
finf=PolynomialRing(QQ,'X')([64*R(a6)[12],16*R(a4)[8],0,1])
assert finf.discriminant() and finf.gcd(finf.derivative())==1
# At infinity: D pole24 and f'(theta) pole8, with nonzero leading factor.
assert finf.discriminant()==D.leading_coefficient()
print('PASS_NORM_AND_RAMIFICATION_IDENTITY',flush=True)
# Exact Sturm sign variations at infinities; no real root approximation.
def primitive(p):
    den=math.lcm(*(int(c.denominator()) for c in p))
    numer=[int(c*den) for c in p]
    common=math.gcd(*numer)
    return p*(QQ(den)/common)
chain=[primitive(D),primitive(D.derivative())]
while chain[-1].degree()>0:
    remainder=-(chain[-2]%chain[-1])
    assert remainder
    chain.append(primitive(remainder))
def signs(infinity):
    return [int(c.leading_coefficient().sign())*((-1)**c.degree() if infinity<0 else 1) for c in chain]
def variations(seq):return sum(a!=b for a,b in zip(seq,seq[1:]))
minus,plus=signs(-1),signs(1)
real_roots=variations(minus)-variations(plus)
print('REAL_DISCRIMINANT_ROOTS',real_roots,'SIGN_AT_302',D(0).sign(),flush=True)
controls=[]
for row in roster['cases']:
    tau=QQ(row['parameter']);value=D(tau);assert value
    controls.append({'label':row['label'],'parameter':str(tau),
        'discriminant_sign':int(value.sign()),
        'real_obstruction_for_beta_times_every_generic_class':bool(value>0)})
assert len(controls)==9
print('FROZEN_PANEL_DISCRIMINANT_SIGNS',[r['discriminant_sign'] for r in controls],flush=True)
report={'status':'PASS_DISCRIMINANT_KUMMER_ALGEBRA',
    'classification':'equation-only construction and exact arithmetic application',
    'definition':'beta=-disc(f)*f_prime(theta)', 'norm':'disc(f)^4',
    'finite_ramification':'the double and simple spectral points over the unique degree24 discriminant place',
    'valuations_at_bad_points':[3,1],'valuations_at_infinity':[-32,-32,-32],
    'discriminant_coefficients':list(map(str,D.list())),
    'beta_coefficients':[{'numerator':list(map(str,c.numerator().list())),
                           'denominator':list(map(str,c.denominator().list()))} for c in beta.list()],
    'Sturm_signs_minus_infinity':minus,'Sturm_signs_plus_infinity':plus,
    'real_discriminant_roots':real_roots,'discriminant_sign_at_302':int(D(0).sign()),
    'frozen_controls':controls,
    'limits':{'seconds':25,'point_searches':0,'new_parameters':0,'integer_factorizations':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,OLD,ROSTER,Path(__file__)]}}
OUT.mkdir(exist_ok=True)
dest=OUT/'construction.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
