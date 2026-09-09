#!/usr/bin/env sage-python
"""Exactly the five old UNKNOWN cover fibres; one factorization per process.

No new cover, parameter, prime sweep, elliptic point search, or production
change. Each --case runs under timeout25s. Protocol freezes the entire
five-case roster before any factorization.
"""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
SOURCE=ART/'det1092_pencil_multiples_v2';OUT=ART/'det1092_pencil_exact_incidence_v1'
CONTROL=SOURCE/'controls.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert json.loads(p.read_text())==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,required=True);args=parser.parse_args()
data=json.loads(CONTROL.read_text());cases=[{k:r[k] for k in ['n','label','parameter']} for r in data['cases'] if r['status']=='UNKNOWN']
assert len(cases)==5 and 0<=args.case<5
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'exact follow-up to five frozen UNKNOWN cover-incidence cases',
    'cases':cases,'limits':{'seconds_per_case':25,'cases':5,'factorizations_per_case':1,
                          'new_parameters':0,'new_covers':0,'point_searches':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CONTROL,SOURCE/'multiple2.json',SOURCE/'multiple3.json',Path(__file__)]},
    'selection':'All and only the five UNKNOWN entries in the old control certificate, in its saved order. No residual point input.'})
c=cases[args.case];row=json.loads((SOURCE/('multiple%d.json'%c['n'])).read_text())
R=PolynomialRing(QQ,'z');z=R.gen();N=R(row['t_of_z']['numerator']);D=R(row['t_of_z']['denominator'])
assert N.gcd(D).degree()==0
tau=QQ(c['parameter']);f=N-tau*D
f*=lcm([q.denominator() for q in f.list()]);f/=gcd([ZZ(q) for q in f.list()]);f=R(f)
if f.leading_coefficient()<0:f=-f
save('input-%02d.json'%args.case,{'case':c,'primitive_incidence_polynomial':list(map(str,f.list())),
    'map_degree':row['parameter_map_degree'],'infinity_preimage':f.degree()<row['parameter_map_degree']})
fac=f.factor();assert fac.value()==f
roots=[-g[0]/g[1] for g,e in fac if g.degree()==1]
for root in roots:assert f(root)==0 and D(root)!=0 and N(root)/D(root)==tau
out={'classification':'exact cover-incidence factorization, not a rank certificate',
    'case':c,'status':'RATIONAL_PREIMAGE' if roots or f.degree()<row['parameter_map_degree'] else 'NO_RATIONAL_PREIMAGE',
    'unit':str(fac.unit()),'factors':[{'coefficients':list(map(str,g.list())),'degree':int(g.degree()),'exponent':int(e)} for g,e in fac],
    'rational_roots':list(map(str,roots)),'infinity_preimage':f.degree()<row['parameter_map_degree'],
    'boundary':'A rational preimage still needs a chart-limit map and independence check; no rational preimage excludes only this fixed cover.',
    'checker_sha256':sha(Path(__file__))}
save('case-%02d.json'%args.case,out)
print(args.case,c['label'],c['n'],out['status'],'factors',[(g.degree(),e) for g,e in fac],
      'roots',list(map(str,roots)),flush=True)
