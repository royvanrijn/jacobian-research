#!/usr/bin/env python3
"""Independent, Sage-free point/rank and snapshot-novelty certificates.

No search score, numerical height, MWState assertion, or descent result is
trusted. Good finite quotients and the no-rational-2-torsion witness are
recomputed using the older standalone implementation.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'elliptic-curves'),str(ROOT/'elliptic-curves/cas')]
sys.set_int_max_str_digits(0)
from elliptic_candidate_record import is_on_weierstrass_curve, weierstrass_invariants
from mod2_reduction_independence import (find_mod2_reduction_certificate,combined_mod2_rank,
    find_two_torsion_certificate_prime,short_curve_has_no_rational_2_torsion_modular_certificate,
    mod2_reduction_signature)
from ecsearch.q12o5867_specialization import load_q12o5867_data,evaluate_projective_specialization

MODEL=ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
SECTIONS=ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json'
DATABASE=ROOT/'artifacts/local/elliptic-curves/breakthrough-audit-20260905/icarm-database-20260905.json'


def read(path):return json.loads(path.read_bytes())
def hashed(path):return sha256(path.read_bytes()).hexdigest()
def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n');temp.replace(path)


def square_root(q):
    if q<0:return None
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    return F(a,b) if a*a==q.numerator and b*b==q.denominator else None


def isomorphic(left,right):
    """Exact characteristic-zero Q-isomorphism test for nonzero c4,c6.

    All candidates in this family trial are checked to be in this branch.
    Special j=0,1728 database rows then cannot match their j-invariants.
    """
    a,b=weierstrass_invariants(left),weierstrass_invariants(right)
    if not a['c4'] or not a['c6'] or not a['discriminant']:raise ValueError('candidate outside nonexceptional-j protocol')
    if a['c4']**3*b['discriminant']!=b['c4']**3*a['discriminant']:return False
    if not b['c4'] or not b['c6'] or not b['discriminant']:return False
    u=square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
    return u is not None and a['c4']==u**4*b['c4'] and a['c6']==u**6*b['c6']


def checked_rank(model,points,primes=None,torsion_prime=None):
    if not points or any(not is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('point membership failed')
    if primes is None:sigs=find_mod2_reduction_certificate(model,points,prime_bound=1000)
    else:sigs=tuple(mod2_reduction_signature(model,points,p) for p in primes)
    if combined_mod2_rank(sigs,len(points))!=len(points):raise ArithmeticError('independent finite quotient columns not established')
    torsion_prime=torsion_prime or find_two_torsion_certificate_prime(model,prime_bound=200)
    if not short_curve_has_no_rational_2_torsion_modular_certificate(model,torsion_prime):raise ArithmeticError('2-torsion witness failed')
    return {'rank_lower_bound':len(points),'no_rational_2_torsion_prime':torsion_prime,
        'signatures':[asdict(s) for s in sigs],
        'argument':'Every integral relation is divisible by 2; E(Q)[2]=0 permits infinite descent. Hence all listed points are independent.'}


def family_check(parameter,model,points):
    f=load_q12o5867_data(MODEL,SECTIONS);t=F(parameter)
    s=evaluate_projective_specialization(f,t.numerator,t.denominator)
    if not isomorphic(s.model,model):raise ArithmeticError('curve is not the claimed family specialization')
    a,b=weierstrass_invariants(s.model),weierstrass_invariants(model)
    u=square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
    if u is None:raise ArithmeticError('missing rational scale')
    for sign in (1,-1):
        scale=sign*u
        if tuple((x/scale**2,y/scale**3) for x,y in s.points)==tuple(points[:17]):return str(scale)
    raise ArithmeticError('first seventeen points do not equal the transported sections')


def build(directory,output):
    directory=directory.resolve()
    if output.exists():raise FileExistsError('use a new immutable certificate path')
    db=read(DATABASE)
    projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in db['curves']]
    protocol=read(directory/'protocol.json');pop=read(directory/'population.json')
    candidates=[]
    for folder in sorted(directory.glob('candidate-*')):
        if not folder.is_dir() or not (folder/'result.json').exists():continue
        result=read(folder/'result.json')
        witnesses=[]
        if 'final_state' in result:
            witnesses.append((len(result['final_state']['state']['reductions']['points']),
                result['curve'],result['final_state']['state']['reductions']['points'],folder/'result.json'))
        # A timed-out worker may retain a stronger basis in an immutable chart
        # than its last top-level checkpoint. Never infer rank from its length.
        for path in sorted((folder/'charts').rglob('*.json')):
            if 'previous' in path.parts:continue
            row=read(path)
            if 'record' not in row:continue
            inp=row['record']['input'];pts=[[p['x'],p['y']] for p in inp['subgroup']]
            witnesses.append((len(pts),inp['curve'],pts,path))
        if not witnesses:continue
        count,rawmodel,rawpoints,source=max(witnesses,key=lambda row:row[0])
        if count<22:continue
        model=tuple(map(F,rawmodel));points=tuple(tuple(map(F,p)) for p in rawpoints)
        scale=family_check(result['parameter'],model,points)
        cert=checked_rank(model,points)
        matches=[r['id'] for r in projection if isomorphic(model,r['ainvs'])]
        row={'id':folder.name,'parameter':result['parameter'],'curve':list(map(str,model)),
            'points':[list(map(str,p)) for p in points],'rank_certificate':cert,
            'family_to_curve_scale_u':scale,'icarm_snapshot_isomorphism_matches':matches,
            'discovery_witness':{'path':str(source.relative_to(ROOT)),'sha256':hashed(source)},
            'search_completion_at_export':result['status'],'completed_chart_records_at_export':len(result['charts'])}
        candidates.append(row);print('CERTIFIED',folder.name,result['parameter'],count,'snapshot_matches',matches,flush=True)
    doc={'schema':'elliptic-curves.compact-r17-new-rank-certificates.v1','status':'PASS_EXACT_RANK_LOWER_BOUNDS',
        'curves':candidates,'protocol_sha256':hashed(directory/'protocol.json'),'population_sha256':hashed(directory/'population.json'),
        'source_hashes':{str(p.relative_to(ROOT)):hashed(p) for p in
            (Path(__file__).resolve(),MODEL,SECTIONS,ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')},
        'novelty_snapshot':{'url':'https://elliptic-rank.icarm.cloud/database.json','date':'2026-09-05',
            'raw_sha256':hashed(DATABASE),'curve_count':len(projection),'equation_projection':projection,
            'acknowledgement':'ICARM and NSF Grant DMS 2425401'},
        'claim_boundary':'Rank lower bounds only. No exact rank or conductor claim. No match in this pinned database is not proof that a curve has never appeared anywhere.'}
    write(output,doc)


def check(path):
    doc=read(path)
    for name,h in doc['source_hashes'].items():
        if hashed(ROOT/name)!=h:raise ArithmeticError('certificate dependency changed: '+name)
    seen=[]
    for row in doc['curves']:
        model=tuple(map(F,row['curve']));points=tuple(tuple(map(F,p)) for p in row['points'])
        old=row['rank_certificate']
        actual=checked_rank(model,points,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
        if json.dumps(actual,sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('finite certificate mismatch')
        if family_check(row['parameter'],model,points)!=row['family_to_curve_scale_u']:raise ArithmeticError('section transport changed')
        matches=[r['id'] for r in doc['novelty_snapshot']['equation_projection'] if isomorphic(model,r['ainvs'])]
        if matches!=row['icarm_snapshot_isomorphism_matches']:raise ArithmeticError('snapshot comparison changed')
        if any(isomorphic(model,m) for m in seen):raise ArithmeticError('duplicate candidate curve')
        seen.append(model);print('REPLAYED',row['parameter'],'rank >=',len(points),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True);g.add_argument('--directory',type=Path);g.add_argument('--check',type=Path)
    p.add_argument('--output',type=Path);a=p.parse_args()
    if a.check:check(a.check)
    else:build(a.directory,a.output)


if __name__=='__main__':main()
