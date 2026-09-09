#!/usr/bin/env sage-python
"""Exact conic seed sequence indistinguishable on a frozen local footprint.

Reuse one dependent conic address and the old positive u=0 certificate.
No new rational specialization is evaluated; no search is launched.
"""
import hashlib, json, signal
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, gcd, lcm, prod, crt
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OLD=ART/'det1092_conic_seed_progression_v1'
OUT=ART/'det1092_local_shadow_seeds_v2'
FRAME=ART/'det1092_split_descent_v1/dependent-conic/generic-frame.json'
CHART=ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
BRANCHES=[ART/('det1092_split_descent_v1/dependent-conic/branch-%d.json'%i) for i in range(2)]
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    OUT.mkdir(exist_ok=True)
    p=OUT/name
    if p.exists(): assert read(p)==data
    else:
        with p.open('x') as stream: json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
inputs=[OLD/'input.json',OLD/'finite-code.json',OLD/'replay.json',FRAME,CHART,*BRANCHES]
save('protocol.json',{
    'classification':'new constructive counterexample to sufficiency of one fixed local footprint',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
    'script_sha256':sha(Path(__file__)),
    'limits':{'seconds':25,'new_evaluated_rational_parameters':0,'point_searches':0,
              'new_prime_exposures':0,'old_proof_prime_range':[181,1009]},
    'rule':'Use all25 primes in the dependent-conic generic frame. Choose the smaller rational preimage in the old conic slope coordinate. Reuse cached positive-u0 rows above179, skipping only exact map/denominator nonunits; retain rank-increasing prime blocks and the first odd-order prime. If these do not prove18, stop unresolved. CRT gives an explicit infinite rational sequence converging in R to the dependent address, identical on the fixed local footprint, and rank18 at disjoint proof primes.',
    'boundary':'The conic and both anchors are old permitted evidence. This does not select a302 member or give a general rank predictor.'})

inp=read(OLD/'input.json'); finite=read(OLD/'finite-code.json')
old_replay=read(OLD/'replay.json'); frame=read(FRAME); chart=read(CHART)
assert old_replay['status']=='PASS_INDEPENDENT_ORACLE_FREE_RANK18_PROGRESSION'
for packet in [inp,finite,old_replay]:
    for name,digest in packet['inputs'].items(): assert sha(ROOT/name)==digest
assert frame['status']=='PASS_GENERIC_ONLY_INJECTIVE_FOOTPRINT'
S=sorted(r['prime'] for r in frame['records']); assert len(S)==25 and max(S)==149
R=PolynomialRing(QQ,'u'); F=R.fraction_field(); u=R.gen()
def dec(d): return F(R(d['numerator']))/R(d['denominator'])
T=dec(inp['base_map']); ai=list(map(dec,inp['weierstrass_functions']))
points=[[dec(d) for d in P] for P in inp['point_functions']]
s=QQ(-528)/3635; aa,bb,cc,dd=map(QQ,chart['parameter_matrix'])
tau=(aa*s+bb)/(cc*s+dd)
roots=(T.numerator()-tau*T.denominator()).roots(QQ,multiplicities=False)
assert len(roots)==2
r=min(roots); assert T(r)==tau
a,b=r.numerator(),r.denominator()
h=cc*s.numerator()+dd*s.denominator()
kappa=h*h/QQ(chart['weierstrass_u'])
Ef=EllipticCurve(F,ai)
AB=[kappa**4*(-Ef.c4()/48),kappa**6*(-Ef.c6()/864)]
pts=[[kappa**2*(x+Ef.b2()/12),kappa**3*(y+(ai[0]*x+ai[2])/2)] for x,y in points]
Ed=EllipticCurve(QQ,[v(r) for v in AB])
assert list(map(str,Ed.a_invariants()))==frame['curve']
at_r=[Ed([x(r),y(r)]) for x,y in pts]
assert [list(map(str,P.xy())) for P in at_r[:17]]==frame['basis']
for x,y in pts: assert y*y==x*x*x+AB[0]*x+AB[1]
matched=[]
for i,path in enumerate(BRANCHES):
    branch=read(path)
    if list(map(str,at_r[-1].xy()))==branch['original']:
        assert branch['status']=='INHERITED_RATIONAL_SPAN'
        m=ZZ(branch['relation_multiplier']); word=list(map(ZZ,branch['relation_word']))
        assert m*at_r[-1]==sum((w*P for w,P in zip(word,at_r[:17])),Ed(0))
        matched.append({'branch':i,'multiplier':int(m),'word':list(map(int,word))})
assert len(matched)==1

def primitive(fs):
    den=lcm([v.denominator() for v in fs]); ps=[R(v*den) for v in fs]
    ps=[R(v/gcd(ps)) for v in ps]
    scale=lcm([c.denominator() for p in ps for c in p]);ps=[R(p*scale) for p in ps]
    content=gcd([ZZ(c) for p in ps for c in p]);return [R(p/content) for p in ps]
coeff_pairs=[primitive([v,F(1)]) for v in AB]
triples=[primitive([x,y,F(1)]) for x,y in pts]
t_pair=primitive([T,F(1)])
def hev(ps,a,b):
    degree=max(p.degree() for p in ps)
    return [ZZ(sum(c*a**j*b**(degree-j) for j,c in enumerate(p))) for p in ps]
coeff_at_r=[hev(pair,a,b) for pair in coeff_pairs]
point_at_r=[hev(P,a,b) for P in triples]
parameter_at_r=hev(t_pair,a,b)
control=[]
for p in S:
    den_orders=[int(v[1].valuation(p)) for v in coeff_at_r]
    point_orders=[int(min(x.valuation(p) for x in P if x)) for P in point_at_r]
    parameter_order=int(min(x.valuation(p) for x in parameter_at_r if x))
    e=1+max([0,*den_orders,*point_orders,parameter_order])
    assert Ed.discriminant().valuation(p)==0
    control.append({'p':p,'exponent':e,'coefficient_denominator_orders':den_orders,
                    'point_common_orders':point_orders,'parameter_common_order':parameter_order})
L=prod(ZZ(v['p'])**v['exponent'] for v in control)

rows=[]; selected=[]; exposures=[]; odd=None; rank=0
for row in finite['exposures']:
    p=row['p']
    if p<=179: continue
    reason=None
    if row['status']!='PASS_CUBIC_KUMMER_REDUCTION': reason=row['status']
    elif a%p==0 or b%p==0: reason='SKIP_ANCHOR_NONUNIT'
    elif any(pair[1](0)%p==0 for pair in coeff_pairs+[t_pair]): reason='SKIP_MAP_DENOMINATOR'
    elif any(all(v(0)%p==0 for v in P) for P in triples): reason='SKIP_POINT_COMMON_ZERO'
    elif kappa.numerator()%p==0 or kappa.denominator()%p==0: reason='SKIP_CONSTANT_MODEL_NONUNIT'
    if reason:
        exposures.append({'p':p,'status':reason});continue
    candidate=matrix(GF(2),len(rows)+len(row['rows']),18,[x for v in rows+row['rows'] for x in v])
    gain=candidate.rank()>rank
    no2=not row['rational_roots'] and odd is None
    if gain or no2:
        selected.append(p);rows+=row['rows'];rank=int(candidate.rank())
        if no2: odd=p
    exposures.append({'p':p,'status':'REUSED_SAVED_PROOF_ROWS','selected':gain or no2,'rank':rank})
    if rank==18 and odd is not None: break
save('proof-selection.json',{'rank':rank,'no_two_torsion_prime':odd,'selected':selected,'exposures':exposures})
assert rank==18 and odd is not None, 'UNRESOLVED_OLD_DISJOINT_PROOF_POOL'
P=prod(ZZ(p) for p in selected)
assert gcd(L,P)==gcd(a*b,P)==1
target=ZZ((-L*b*ZZ(a).inverse_mod(P))%P)
q0=ZZ(crt(ZZ(1),target,L,P)); period=L*P
assert 0<q0<period and q0%L==1 and (a*q0+b*L)%P==0 and gcd(q0,P)==1
report={
    'status':'CONSTRUCTED_UNIFORM_RANK18_LOCAL_SHADOW_SEQUENCE',
    'classification':'new exact constructive deduction, independent replay required',
    'dependent_reduced_parameter':str(s),'dependent_original_parameter':str(tau),
    'dependent_conic_slope':str(r),'other_slope':str(max(roots)),
    'constant_short_model_scale':str(kappa),'dependence':matched[0],
    'control_primes':S,'control_stability':control,
    'proof_primes':selected,'no_two_torsion_prime':odd,
    'L':str(L),'P':str(P),'q0':str(q0),'period':str(period),
    'sequence':'u_n=r+L/(q0+L*P*n), n any integer; elliptic parameter t=T(u_n)',
    'real_limit':'u_n tends to r as n tends to positive or negative infinity',
    'map_coefficients':{
        'short_coefficient_pairs':[[list(map(str,p.list())) for p in pair] for pair in coeff_pairs],
        'projective_points':[[list(map(str,p.list())) for p in triple] for triple in triples],
        'base_parameter_pair':[list(map(str,p.list())) for p in t_pair]},
    'control_anchor_equation':frame['curve'],
    'boundary':'Every term has rank>=18 but the limiting marked branch is inherited. Only the frozen25-prime panel is indistinguishable; disjoint proof primes distinguish them. No claim of rank17 for the limiting elliptic fibre, and no302 coverage.',
    'inputs':read(OUT/'protocol.json')['inputs'],'protocol_sha256':sha(OUT/'protocol.json')}
save('construction.json',report)
print(report['status'],'controls',len(S),'proof primes',selected,'L/P digits',len(str(L)),len(str(P)),flush=True)
