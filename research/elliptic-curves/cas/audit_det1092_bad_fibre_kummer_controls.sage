#!/usr/bin/env sage-python
"""Fixed nine-fibre test of beta; no point inputs or prime factorization.

If disc(f) has odd valuation away from coefficient denominators,2,3,c4,
the reduction is nodal with odd discriminant valuation. Local Kummer
images are unramified there, while beta has odd valuation in both factors.
"""
import hashlib,json,signal,math
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,prod
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_bad_fibre_kummer_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
ROSTER=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
BASE=OUT/'construction.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
parent,roster,base=map(read,[PARENT,ROSTER,BASE])
for name,digest in base['inputs'].items():assert sha(ROOT/name)==digest
R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
ai=list(map(dec,parent['a_invariants']))
assert ai[:3]==[1,1,1]
a,b,c=R(5),R(16*ai[3]+8),R(64*ai[4]+16)
D=R(list(map(QQ,base['discriminant_coefficients'])))
assert D==a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
c4=16*(a*a-3*b)
rows=[]
for old in roster['cases']:
    t0=QQ(old['parameter']);value=D(t0);cc=c4(t0)
    coeff=[c(t0),b(t0),a(t0),QQ(1)]
    assert value and cc
    support_factors=[2,3,int(value.denominator()),abs(int(cc.numerator())),
                     int(cc.denominator()),*[int(v.denominator()) for v in coeff]]
    support=math.prod(support_factors);rem=abs(int(value.numerator()));strips=[]
    while True:
        g=math.gcd(rem,support)
        if g==1:break
        strips.append(str(g));rem//=g
    rr=math.isqrt(rem);finite=(rr*rr!=rem)
    rows.append({'label':old['label'],'parameter':str(t0),
        'discriminant':str(value),'c4':str(cc),'cubic_coefficients':list(map(str,coeff)),
        'real_obstruction':bool(value>0),
        'support_factors':list(map(str,support_factors)),
        'gcd_strips':strips,'coprime_remainder':str(rem),'floor_sqrt':str(rr),
        'odd_nodal_prime_obstruction':finite,
        'status':'NOT_SELMER' if value>0 or finite else 'UNKNOWN_NO_OBSTRUCTION_IN_THIS_GATE'})
    print(old['label'],'real',value>0,'odd_nodal_prime',finite,flush=True)
# The smooth infinity cubic, using only the old fixed reduction primes.
finf=PolynomialRing(QQ,'X')([c[12],b[8],0,1]);infinity=[]
for p in [149,151]:
    T=PolynomialRing(GF(p),'X');h=T(finf)
    assert h.degree()==3 and h.discriminant()
    values=[int(h(i)) for i in range(p)]
    infinity.append({'prime':p,'coefficients':[int(v) for v in h],
                     'values':values,'irreducible':all(values)})
print('INFINITY_IRREDUCIBILITY',[(r['prime'],r['irreducible']) for r in infinity],flush=True)
report={'status':'PASS_FIXED_CONTROL_ALGEBRA',
    'classification':'equation-only local obstruction; no solubility claim from gate survival',
    'rows':rows,'infinity_cubic_coefficients':list(map(str,finf.list())),
    'infinity_modular_checks':infinity,
    'limits':{'seconds':25,'point_searches':0,'new_parameters':0,'integer_factorizations':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,ROSTER,BASE,Path(__file__)]}}
dest=OUT/'controls.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
