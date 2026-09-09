#!/usr/bin/env sage-python
"""Independent algebraic replay of the bad-fibre-only Kummer obstruction.

No producer import, exceptional point, new parameter, factoring of integers,
point search, or generic Selmer/group computation. The valuation proof and
the prior all-unramified cohomological equality remain written proofs.
"""
import hashlib,json,signal,math
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,prod,power_mod
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_bad_fibre_kummer_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(d):
    for name,digest in d['inputs'].items():assert sha(ROOT/name)==digest
protocol=read(OUT/'protocol.json');provenance(protocol)
base=read(OUT/'construction.json');controls=read(OUT/'controls.json')
provenance(base);provenance(controls)
unram=read(ART/'det1092_unramified_kummer_v1/replay.json')
assert unram['status']=='PASS_HYPOTHESES_AND_RESIDUE_REGRESSION';provenance(unram)
prior=read(ART/'det1092_polynomial_lift_gate_v1/replay.json')
assert prior['status']=='PASS_INDEPENDENT_POLYNOMIAL_LIFT_OBSTRUCTION';provenance(prior)
old=read(ART/'det1092_polynomial_lift_gate_v1/construction.json')
parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
roster=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')
R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
ai=list(map(dec,parent['a_invariants']));assert ai[:3]==[1,1,1]
a,b,c=R(5),R(16*ai[3]+8),R(64*ai[4]+16)
S=PolynomialRing(K,'z');z=S.gen();f=z**3+a*z*z+b*z+c
D=R(f.discriminant());assert D.degree()==24 and D.gcd(D.derivative())==1
assert list(map(str,D.list()))==base['discriminant_coefficients']
possible=set(range(25))
for saved,p in zip(old['primes'],[149,151]):
    assert saved['prime']==p
    T=PolynomialRing(GF(p),'t');t=T.gen();dp=T(D.monic())
    factors=[T(h) for h in saved['factors']]
    assert prod(factors)==dp and dp.gcd(dp.derivative())==1
    for h in factors:
        n=h.degree();assert h.is_monic() and power_mod(t,p**n,h)==t%h
        for k in range(1,n):
            if n%k==0:assert (power_mod(t,p**k,h)-t).gcd(h)==1
    sums={0}
    for h in factors:sums|={v+h.degree() for v in sums}
    possible&=sums
assert possible=={0,24}
roots={row['label']:R(list(map(QQ,row['root_remainder']))) for row in old['templates']}
r,s=roots['double'],roots['simple']
assert (r-s).gcd(D)==1 and all(R(v)%D==0 for v in (f-(z-r)**2*(z-s)).list())
beta=S([dec(row) for row in base['beta_coefficients']]);assert beta==-D*f.derivative()
# Multiplication determinant rather than the constructor's resultant.
mult=matrix(K,3,3,lambda i,j:((beta*z**j)%f)[i])
assert mult.det()==K(D)**4
finf=PolynomialRing(QQ,'X')([c[12],b[8],0,1])
assert finf.discriminant()==D.leading_coefficient()!=0
assert base['valuations_at_bad_points']==[3,1]
assert base['valuations_at_infinity']==[-32,-32,-32]

# Independent normalized negative-remainder recurrence, using positive
# leading coefficient rather than the producer's primitive integral content.
chain=[D/abs(D.leading_coefficient()),D.derivative()/abs(D.derivative().leading_coefficient())]
while chain[-1].degree()>0:
    h=-(chain[-2]%chain[-1]);assert h
    chain.append(h/abs(h.leading_coefficient()))
neg=[int(h.leading_coefficient().sign())*(-1)**h.degree() for h in chain]
pos=[int(h.leading_coefficient().sign()) for h in chain]
assert neg==base['Sturm_signs_minus_infinity'] and pos==base['Sturm_signs_plus_infinity']
def changes(seq):return sum(x!=y for x,y in zip(seq,seq[1:]))
assert changes(neg)-changes(pos)==base['real_discriminant_roots']==18

rows=[]
assert len(controls['rows'])==len(roster['cases'])==9
for saved,source in zip(controls['rows'],roster['cases']):
    t0=QQ(source['parameter']);value=D(t0);c4=16*(a(t0)**2-3*b(t0))
    coeff=[c(t0),b(t0),a(t0),QQ(1)]
    assert saved['label']==source['label'] and QQ(saved['parameter'])==t0
    assert value==QQ(saved['discriminant']) and c4==QQ(saved['c4'])
    assert coeff==list(map(QQ,saved['cubic_coefficients']))
    support_factors=[ZZ(2),ZZ(3),value.denominator(),abs(c4.numerator()),c4.denominator(),
                     *[v.denominator() for v in coeff]]
    assert support_factors==list(map(ZZ,saved['support_factors']))
    support=ZZ(prod(support_factors));rem=abs(value.numerator())
    for val in saved['gcd_strips']:
        divisor=ZZ(val);assert rem.gcd(support)==divisor and divisor>1
        rem//=divisor
    assert rem==ZZ(saved['coprime_remainder'])
    gg,u,v=rem.xgcd(support);assert gg==1 and u*rem+v*support==1
    rr=ZZ(saved['floor_sqrt']);assert rr*rr<rem<(rr+1)**2
    assert saved['odd_nodal_prime_obstruction'] is True and saved['status']=='NOT_SELMER'
    assert saved['real_obstruction']==bool(value>0)
    rows.append({'label':source['label'],'real_obstruction':bool(value>0),
        'odd_nodal_prime_obstruction':True,'coprime_remainder':str(rem),
        'floor_sqrt':str(rr),'Bezout_with_support':[str(u),str(v)],
        'beta_times_every_generic_class':'NOT_SELMER',
        'any_rational_Kummer_specialization_from_A_bad':'inherited mod2 image only'})
assert sum(r['real_obstruction'] for r in rows)==7
assert [r['label'] for r in rows if not r['real_obstruction']]==['scale-0131232','scale-0487239']
for saved,p in zip(controls['infinity_modular_checks'],[149,151]):
    T=PolynomialRing(GF(p),'X');h=T(finf)
    values=[int(h(i)) for i in range(p)]
    assert saved['values']==values and saved['irreducible']==all(values)==False

report={'status':'PASS_INDEPENDENT_BAD_FIBRE_KUMMER_OBSTRUCTION',
    'classification':'new deduction and independently replayed equation-only local obstructions',
    'discriminant_irreducible':True,'bad_only_coset_representatives':['1','-Delta*f_prime(theta)'],
    'rows':rows,'real_discriminant_roots':18,
    'cohomology_dependency':'All-unramified norm kernel equals delta(M17), proved in the prior canonical note.',
    'local_proof_gate':'At an odd nodal prime with odd discriminant valuation n, every rational Kummer image has even valuations; beta has (3n,n). Proven by Hensel simple-root factorization and the parity of min(2m,n).',
    'boundary':'No specialized rank upper bound, no exclusion of nongeneric doubled points, and no solubility claim from the absence of either obstruction. Smooth-fibre ramification and base changes remain open.',
    'inputs':protocol['inputs'],'protocol_sha256':sha(OUT/'protocol.json'),'checker_sha256':sha(Path(__file__))}
dest=OUT/'independent-replay.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],'real',sum(r['real_obstruction'] for r in rows),'odd_nodal',len(rows),flush=True)
