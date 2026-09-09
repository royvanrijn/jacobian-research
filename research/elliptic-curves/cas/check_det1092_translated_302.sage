#!/usr/bin/env sage-python
"""Exact302 incidence of the completed six generic translation maps, <=25s."""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_translated_sections_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert json.loads(path.read_text())==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,required=True);a=parser.parse_args();assert 0<=a.case<6
maps=[OUT/('map-%02d.json'%i) for i in range(6)];assert all(p.exists() for p in maps)
save('302-incidence-protocol.json',{'classification':'retrospective positive-control evaluation, after generic maps freeze',
    'parameter':'0','cases':list(range(6)),
    'limits':{'seconds_per_case':25,'cases':6,'rational_factorizations_per_case':1,
              'point_searches':0,'exceptional_point_inputs':0,'map_retuning':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',*maps,Path(__file__)]}})
d=json.loads(maps[a.case].read_text());R=PolynomialRing(QQ,'u')
N=R(d['t_of_u']['numerator']);D=R(d['t_of_u']['denominator']);assert N.gcd(D).degree()==0
f=N*lcm([c.denominator() for c in N.list()]);f/=gcd([ZZ(c) for c in f.list()]);f=R(f)
if f.leading_coefficient()<0:f=-f
fac=f.factor();assert fac.value()==f
roots=[str(-g[0]/g[1]) for g,e in fac if g.degree()==1]
infinity=f.degree()<d['degree']
record={'status':'RATIONAL_PREIMAGE_PENDING_POINT_MAP' if roots or infinity else 'NO_RATIONAL_PREIMAGE',
    'classification':'exact302 incidence; not an elliptic rank bound',
    'case':d['case'],'primitive_polynomial':list(map(str,f.list())),
    'degree':d['degree'],'factor_degrees':[[int(g.degree()),int(e)] for g,e in fac],
    'rational_roots':roots,'infinity_preimage':infinity,'checker_sha256':sha(Path(__file__)),
    'boundary':'Only these six fixed rational curves are tested. No original point coordinates are read.'}
save('302-incidence-%02d.json'%a.case,record)
print(a.case,d['case'],record['status'],'roots',roots,'degrees',record['factor_degrees'],flush=True)
