#!/usr/bin/env sage-python
"""Exact splitting only on the unchanged nine original-parameter controls."""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,prime_range
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rational_bisection_index_v1'
ROSTER=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=PolynomialRing(QQ,'t');cases=read(ROSTER)['cases'];assert len(cases)==9
rows=[]
for orbit in [8044,47755,103186]:
    data=read(OUT/('orbit-%d.json'%orbit));q=R(data['q'])
    outcome=[]
    for case in cases:
        tau=QQ(case['parameter']);value=QQ(q(tau))
        split=value.is_square();witness=None
        if not split:
            for p in prime_range(3,180):
                v=value.valuation(p);unit=value/p**v
                c=ZZ(unit.numerator()*unit.denominator().inverse_mod(p)%p)
                if v%2 or pow(int(c),(int(p)-1)//2,int(p))==p-1:
                    witness={'p':int(p),'valuation':int(v),'unit_residue':int(c)};break
        outcome.append({'label':case['label'],'parameter':str(tau),'q_value':str(value),
                        'split':bool(split),'local_nonsquare_witness':witness,
                        'rational_ordinate':str(value.sqrt()) if split else None,
                        'independence':'NOT_TESTED' if split else 'NO_POINT_FROM_THIS_COVER'})
    rows.append({'orbit':orbit,'cases':outcome})
report={'status':'PASS_THREE_FIXED_BISECTION_SPLITTING_AUDIT','classification':'retrospective fixed-control evaluation after generic construction',
    'limits':{'seconds':25,'curves':3,'old_parameters':9,'new_parameters':0,'point_searches':0,'local_prime_bound':179},
    'rows':rows,'boundary':'Split branches require a separate independence certificate. Nonsplitting excludes only the named orbit and all its generic translations, not other seeds.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ROSTER,Path(__file__),*[OUT/('orbit-%d.json'%i) for i in [8044,47755,103186]]]}}
dest=OUT/'controls.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],[(r['orbit'],[v['label'] for v in r['cases'] if v['split']]) for r in rows],flush=True)
