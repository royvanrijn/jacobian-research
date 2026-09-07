#!/usr/bin/env sage-python
"""Open masked oracles only after search, certify exact recovered directions.

Finite reductions propose bounded integer/rational group words; exact group
law verifies them. Failure remains UNKNOWN. No numerical heights or search.
"""
import gzip,json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys
from sage.all import EllipticCurve,QQ,ZZ,GF,matrix,vector,lcm,prime_range

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from mod_l_reduction_independence import find_mod_l_reduction_certificate,mod_l_reduction_signature,combined_mod_l_rank
from research_runtime.store import checkpoint
from search_observability import point_visibility

ART=ROOT/'artifacts/generated-results/elliptic-curves'
SEARCH=ART/'ordinary_masked_controls_v1.json.gz'
POP=ART/'fibre_height_population_v1.json.gz'


def read(p):return json.loads(gzip.decompress(p.read_bytes()))

def run():
    from dataclasses import asdict
    search=read(SEARCH);population=read(POP)
    rows={r['id']:r for f in population['families'].values() for r in f['rows']}
    output=dict(schema='elliptic-curves.ordinary-masked-relations.v1',results={},
        inputs={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in (SEARCH,POP,Path(__file__))},
        endpoint='WITHHELD_KNOWN_DIRECTIONS_NOT_NEW_RANK',moduli=[3,5,7,11],prime_bound=1000)
    for identifier,result in search['results'].items():
        row=rows[identifier];n=len(row['subgroup']);model=tuple(map(Fraction,row['search_model']))
        E=EllipticCurve(QQ,row['search_model']);points=[E(QQ(p['x']),QQ(p['y'])) for p in row['subgroup']]
        records={tuple(p.values()):p for ch in result['charts'] for p in ch['finite_curve_points']}
        discovered=[E(QQ(p['x']),QQ(p['y'])) for p in records.values()]
        allpoints=tuple(tuple(Fraction(str(c)) for c in p.xy()) for p in points+discovered)
        words={};attempts=[]
        for modulus in output['moduli']:
            if len(words)==len(discovered):break
            signatures=[]; proposal_rows=[]
            for prime in prime_range(3,1001):
                try:s=mod_l_reduction_signature(model,allpoints,int(prime),modulus)
                except ValueError:continue
                if not s.rows:continue
                signatures.append(s);proposal_rows.extend(s.rows)
                if matrix(GF(modulus),proposal_rows)[:,:n].rank()==n:break
            if not proposal_rows or matrix(GF(modulus),proposal_rows)[:,:n].rank()!=n:
                attempts.append(dict(modulus=modulus,status='INSUFFICIENT_FINITE_RANK'));continue
            m=matrix(GF(modulus),[r for s in signatures for r in s.rows]);base=m[:,:n]
            attempts.append(dict(modulus=modulus,status='FULL_INPUT_RANK',signatures=[asdict(s) for s in signatures]))
            for i,point in enumerate(discovered):
                if str(i) in words:continue
                try:coefficients=base.solve_right(m.column(n+i))
                except ValueError:continue
                candidates=[]
                candidates.append([QQ(int(c) if int(c)<=modulus//2 else int(c)-modulus) for c in coefficients])
                # Scalar denominators 2,3,4 cover a modest saturation relation;
                # verified exact equality, not a rank inference from congruences.
                for d in (2,3,4):
                    if d%modulus:
                        integers=[int(d*c) if int(d*c)<=modulus//2 else int(d*c)-modulus for c in coefficients]
                        candidates.append([QQ(c)/d for c in integers])
                for word in candidates:
                    denominator=lcm(c.denominator() for c in word)
                    rhs=sum((ZZ(denominator*c)*p for c,p in zip(word,points) if c),E(0))
                    if denominator*point==rhs:
                        words[str(i)]=dict(point=list(map(str,point.xy())),coefficients=list(map(str,word)),
                            relation_denominator=int(denominator),withheld_coefficient=str(word[0]),exact_group_relation=True)
                        break
        # A returned partner in a pointed chart is often Q-R. Use only
        # retained centre words and already exactly reconstructed points.
        # Check the group law explicitly, rather than infer a word from x.
        for i,point in enumerate(discovered):
            if str(i) in words:continue
            for ch in result['charts']:
                centre=ch['input']['centre']['coefficients']
                centre=vector(QQ,[0]+centre)
                for known in list(words.values()):
                    for sign in (1,-1):
                        word=centre+sign*vector(QQ,known['coefficients'])
                        d=lcm(c.denominator() for c in word)
                        if d*point==sum((ZZ(d*c)*p for c,p in zip(word,points) if c),E(0)):
                            words[str(i)]=dict(point=list(map(str,point.xy())),coefficients=list(map(str,word)),
                                relation_denominator=int(d),withheld_coefficient=str(word[0]),exact_group_relation=True,
                                proposal='retained centre plus signed exactly reconstructed point')
                            break
                    if str(i) in words:break
                if str(i) in words:break
        unknown=[i for i in range(len(discovered)) if str(i) not in words]
        recovery=any(QQ(r['withheld_coefficient']) for r in words.values())
        coverage=[]
        for j,ch in enumerate(result['charts']):
            for p in ch['finite_curve_points']:
                v=point_visibility(ch,p)
                if v['status'] not in ('VISIBLE_AND_RECORDED','KNOWN_POINTED_ENDPOINT'):
                    raise ArithmeticError('discovered point failed coverage replay')
                coverage.append(v['status'])
        output['results'][identifier]=dict(family=result['family'],finite_point_count=len(discovered),
            exact_relations=words,finite_reduction_proposals=attempts,unknown_point_indices=unknown,
            withheld_direction_recovered=bool(recovery),
            withheld_direction_status='RECOVERED' if recovery else 'UNKNOWN' if unknown else 'BOUNDED_MISS',
            coverage_replays=len(coverage))
        print('RELATIONS',identifier,'points',len(discovered),'recovered',bool(recovery),'unknown',len(unknown),flush=True)
    output['status']='COMPLETE_EXACT_RELATIONS' if not any(r['unknown_point_indices'] for r in output['results'].values()) else 'UNKNOWN_RELATIONS_REMAIN'
    checkpoint(ART/'ordinary_masked_relations_v1.json',output)


if __name__=='__main__':run()
