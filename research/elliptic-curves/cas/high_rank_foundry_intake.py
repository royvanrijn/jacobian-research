"""Pinned equation intake, exact Q-deduplication, and inexhaustible height shells."""
from fractions import Fraction as F
import json
import math
from pathlib import Path

import certify_compact_r17_candidates as cert
import compact_atlas_specialization as atlas
from high_rank_foundry_policy import FAMILIES, hashed


def jkey(model):
    inv = cert.weierstrass_invariants(tuple(map(F,model)))
    if not inv['discriminant']:
        raise ValueError('singular equation')
    return str(inv['c4']**3/inv['discriminant'])


def supported(model):
    inv = cert.weierstrass_invariants(tuple(map(F,model)))
    return bool(inv['discriminant'] and inv['c4'] and inv['c6'])


def matches(model, others):
    return any(cert.isomorphic(tuple(map(F,model)), tuple(map(F,q))) for q in others)


def shell(height):
    return [(p,q) for q in range(1,height+1) for p in
            (range(-height,height+1) if q == height else (-height,height))
            if math.gcd(p,q) == 1]


def streams(pool):
    result = {}
    for family in FAMILIES:
        rows = [r for r in pool if r['family']==family]
        ordered = sorted(rows,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
        result[family] = {
            'low': [r for r in ordered if max(abs(r['numerator']),r['denominator'])<=1024],
            'high': [r for r in ordered if max(abs(r['numerator']),r['denominator'])>1024],
            'all': sorted(rows,key=lambda r:hashed('score-independent-v1',family,r['parameter'])),
        }
    return result


class Intake:
    def __init__(self, snapshot, cursor=None):
        self.snapshot = snapshot
        self.streams = streams(snapshot['pool'])
        self.families = {f['family']:f for f in json.loads(atlas.ATLAS.read_text())['families']}
        self.cursor = cursor or {'ordinal':0,'offsets':{},'shells':{}}
        self.excluded = {}
        for q in snapshot['excluded_models']:
            try:
                self.excluded.setdefault(jkey(q),[]).append(q)
            except ValueError:
                pass
        self.addresses = {tuple(a) for a in snapshot['reserved_addresses']}

    def next(self, seen):
        """Return a fresh row and an explicit rejection ledger; seen maps j to models."""
        rejected = []
        while True:
            n = self.cursor['ordinal']
            family = FAMILIES[n%6]
            lane = ('low','low','all','low','high','raw','all','raw')[(n//6)%8]
            if lane != 'raw':
                key = family+'/'+lane
                offset = self.cursor['offsets'].get(key,0)
                available = self.streams[family][lane]
                self.cursor['offsets'][key] = offset+1
                if offset < len(available):
                    row = dict(available[offset])
                else:
                    lane = 'raw'
            if lane == 'raw':
                h,i = self.cursor['shells'].get(family,[1,0])
                pairs = shell(h)
                if i == len(pairs):
                    h,i = h+1,0; pairs = shell(h)
                p,q = pairs[i]
                self.cursor['shells'][family] = [h,i+1]
                row = {'family':family,'parameter':str(F(p,q)),'numerator':p,'denominator':q,
                       'combined_selection_units':None,'combined_good':None}
                try:
                    model,_ = atlas.specialize(self.families[family],row['parameter'])
                    row['model'] = list(map(str,model))
                except (ArithmeticError,ZeroDivisionError) as e:
                    rejected.append({'family':family,'parameter':row['parameter'],'reason':'UNSUPPORTED_SPECIALIZATION','detail':str(e)})
                    continue
            addr = family,str(F(row['parameter']))
            if addr in self.addresses:
                rejected.append({'family':family,'parameter':addr[1],'reason':'RESERVED_ADDRESS'})
                continue
            if not supported(row['model']):
                rejected.append({'family':family,'parameter':addr[1],'reason':'UNSUPPORTED_EQUATION_NOT_RANK_EXCLUSION'})
                continue
            jk = jkey(row['model'])
            if matches(row['model'], self.excluded.get(jk,[])+seen.get(jk,[])):
                rejected.append({'family':family,'parameter':addr[1],'reason':'EXACT_Q_ISOMORPH'})
                continue
            self.cursor['ordinal'] += 1
            row.update(id='f-'+hashed(family,addr[1])[:20],parameter=addr[1],lane=lane,
                       selection_ordinal=n,parameter_height=max(abs(row['numerator']),row['denominator']),j_key=jk)
            return row,rejected


def novelty(model, catalogue):
    jk = jkey(model)
    found = [r['id'] for r in catalogue['curves'] if jkey(r['ainvs']) == jk and matches(model,[r['ainvs']])]
    return {'status':'MATCHES_PINNED_CATALOGUE' if found else 'ABSENT_FROM_PINNED_CATALOGUE',
            'matches':found,'catalogue_count':catalogue['count'],
            'highest_reported_lower_bound':max((r['rank_lower_bound'] for r in catalogue['curves'] if r['id'] in found),default=None),
            'boundary':'Exact Q-isomorphism comparison to a frozen snapshot; worldwide novelty remains UNKNOWN.'}


def conductor_gate(model, rank, catalogue):
    thresholds = [int(r['conductor']) for r in catalogue['curves']
                  if r['rank_lower_bound']>=rank and r.get('conductor')]
    if not thresholds:
        return {'eligible':False,'reason':'NO_REPORTED_BENCHMARK','threshold':None}
    threshold = min(thresholds)
    # A small raw discriminant can flag inexpensive fibres without factoring.
    # This is only a positive scheduling gate. Large delta is not a lower bound on N.
    q=tuple(map(F,model));scale=1
    if q[:3]==(0,0,0) and q[3].denominator==q[4].denominator==1 and q[3] and q[4]:
        from mod2_reduction_independence import _primes_up_to
        a,b=abs(q[3].numerator),abs(q[4].numerator)
        for p in _primes_up_to(1000):
            while a % (p**4)==0 and b % (p**6)==0:
                a//=p**4;b//=p**6;scale*=p
    delta = abs(cert.weierstrass_invariants(q)['discriminant'])/scale**12
    eligible = delta <= 100*threshold
    return {'eligible':eligible,'reason':'SMALL_RAW_DISCRIMINANT' if eligible else 'LOW_CONDUCTOR_PRIORITY',
            'threshold':str(threshold),'discriminant_abs':str(delta),'proven_short_model_scale':str(scale),
            'missing_catalogue_conductors':sum(r['rank_lower_bound']>=rank and not r.get('conductor') for r in catalogue['curves'])}
