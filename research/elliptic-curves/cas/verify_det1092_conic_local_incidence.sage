#!/usr/bin/env sage-python
"""Nine fixed conic non-incidences with explicit odd-prime open balls.

No factorization, new parameter, or point search. Prime witness pool<=179.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,prime_range,ceil
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_local_shadow_seeds_v2'
COVER=ART/'det1092_orbit8044_rank18_base_change_v2.json'
ROSTER=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=PolynomialRing(QQ,'t');q=R(read(COVER)['curve_over_Q']['q_coefficients'])
assert q.degree()==2 and q.discriminant()!=0
cases=read(ROSTER)['cases'];assert len(cases)==9
rows=[]
for case in cases:
    tau=QQ(case['parameter']);value=q(tau);assert value!=0
    witnesses=[];exposures=[]
    for p in prime_range(3,180):
        v=ZZ(value.valuation(p));unit=value/p**v
        residue=ZZ(unit.numerator()*unit.denominator().inverse_mod(p)%p)
        square=int(residue) in {int(y*y%p) for y in range(int(p))}
        assert square==(pow(int(residue),(int(p)-1)//2,int(p))==1)
        obstruction=v%2!=0 or not square
        exposures.append({'p':int(p),'valuation':int(v),'unit_residue':int(residue),
                          'unit_is_square_mod_p':square,'obstructed':obstruction})
        if obstruction:
            derivative=q.derivative()(tau)
            N=max([ZZ(1),ZZ(ceil((v+1-q[2].valuation(p))/2))]+
                  ([ZZ(v+1-derivative.valuation(p))] if derivative else []))
            assert q[2].valuation(p)+2*N>v
            assert not derivative or derivative.valuation(p)+N>v
            witnesses.append({'p':int(p),'valuation':int(v),'unit_residue':int(residue),
                              'unit_is_square_mod_p':square,'parameter_ball_exponent':int(N)})
            break
    assert len(witnesses)==1,'UNRESOLVED_FIXED_PRIME_POOL'
    rows.append({'label':case['label'],'parameter':str(tau),'q_value':str(value),
                 'witness':witnesses[0],'exposures':exposures,
                 'conclusion':'No rational splitting of this conic anywhere in the indicated p-adic parameter ball; not an exclusion of other seeds.'})
result={'status':'PASS_NINE_EXPLICIT_LOCAL_CONIC_NONINCIDENCE_BALLS',
        'classification':'verified application of the elementary odd-prime square criterion',
        'rule':'First local obstruction at an odd prime<=179, on the unchanged nine-address roster.',
        'limits':{'seconds':25,'cases':9,'prime_bound':179,'new_parameters':0,'factorizations':0,'point_searches':0},
        'cases':rows,
        'uniformity':'q(tau+delta)-q(tau)=q_prime(tau)*delta+q2*delta^2. The certified delta valuation makes the difference higher than v_p(q(tau)), preserving its obstructed leading squareclass.',
        'boundary':'These balls exclude only orbit8044 splitting. Curve302 already has an independent seed from another carrier. No fibre rank upper bounds follow.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [COVER,ROSTER,Path(__file__)]}}
dest=OUT/'local-nonsplitting.json'
if dest.exists():assert read(dest)==result
else:
    with dest.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(result['status'],[(r['label'],r['witness']['p'],r['witness']['parameter_ball_exponent']) for r in rows],flush=True)
