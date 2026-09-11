#!/usr/bin/env sage-python
"""Exact event audit and support-cached all-subset accessibility experiment.

No point search. Fixed raw slope coordinate, integer height comparisons, and
unchanged held-out points throughout. Stages checkpoint deterministic files;
--check replays events, rational chart costs, and subset minima independently.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gc
import gzip
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd, log2, comb
from pathlib import Path
from statistics import median

from sage.all import QQ, ZZ, EllipticCurve, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
SOURCE = ART / 'rank_accessibility_atlas_curve302_11952_v1.json.gz'
INDEX = ART / 'rank_accessibility_atlas_curve302_11952_v1.summary.json'
OUT = ART / 'rank_accessibility_subsets_v1'
POLICIES = [('original', 0, 0), ('generic_plus', 1, 0),
            ('generic_minus', -1, 0), ('exceptional_plus', 0, 1),
            ('exceptional_minus', 0, -1), ('both_plus', 1, 1)]


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def save(path, value):
    raw = canonical(value)
    if path.suffix == '.gz':
        raw = gzip.compress(raw, compresslevel=6, mtime=0)
    if path.exists():
        require(path.read_bytes() == raw, f'Existing checkpoint differs: {path}')
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + '.tmp')
        temporary.write_bytes(raw)
        temporary.replace(path)


def word(values):
    return tuple(map(int, values))


def axis(n, i):
    return tuple(int(j == i) for j in range(n))


def combine(left, right, sign=1):
    return tuple(a + sign * b for a, b in zip(left, right))


def negate(values):
    return tuple(-x for x in values)


def coords(point):
    return [str(x) for x in point.xy()]


def label(w):
    pieces = []
    for i, c in enumerate(w):
        if c:
            name = f'M{i+1}' if i < 17 else f'E{i-16}'
            pieces.append(('+' if c > 0 else '-') + (str(abs(c)) + '*' if abs(c) != 1 else '') + name)
    return ''.join(pieces).lstrip('+') or '0'


def height(row):
    p = row['parameter']
    require(p != 'infinity', 'Held-out point became an anchor')
    return max(abs(int(p['numerator'])), int(p['denominator']))


def logratio(value):
    return log2(value.numerator) - log2(value.denominator)


def protocol():
    return {
        'schema': 'rank-accessibility-subsets.protocol.v1',
        'inputs': {str(p.relative_to(ROOT)): digest(p) for p in (SOURCE, INDEX, Path(__file__))},
        'metric': 'Raw t=(y(P)+y(Q))/(x(P)-x(Q)); compare reduced integer H=max(abs(n),d) exactly; logs only for display.',
        'normalization': 'Fixed short model and primitive integral common N/F coefficient pair per anchor; no surrounding-subgroup input. Gauss quartics in old atlas do not define its kappa.',
        'policies': POLICIES,
        'rebasing': 'In each selected block, disjoint consecutive pairs (a,b) become (a+s*b,b), s=+1 or -1. Generic indices 0..16; exceptional indices in original order with the fixed hidden target omitted. Both matrices have determinant one and never mix M17 and exceptional blocks.',
        'heldout': 'Original target points remain fixed. No rebased supplied generator has a component in that target direction. Basis choice uses indices only, never coordinates or costs.',
        'limits': {'fibres': 2, 'policies_per_fibre': 6, 'states_per_target_max': 8192,
                   'anchors_per_fixed_full_original_basis_max': 1922,
                   'anchors_per_target_native_bank_max': 1800, 'point_searches': 0},
        'units': 'Targets/fibres and explicitly counted subset edges; descriptive, no independent-sample inference.',
    }


def event_chart(curve, basis, inverse, row):
    E = EllipticCurve(QQ, curve['curve'])
    key = ','.join(map(str, row['anchor_word_in_D']))
    chart = curve['chart_atlas'][key]
    q = E(list(map(QQ, chart['anchor'])))
    n, d = int(row['parameter']['numerator']), int(row['parameter']['denominator'])
    H = max(abs(n), d)
    def evaluate(cs):
        return sum(int(c) * n**i * d**(4-i) for i, c in enumerate(cs))
    pair = chart['primitive_integral_x_map']
    N, F = evaluate(pair['N_ascending']), evaluate(pair['F_ascending'])
    require((str(N), str(F)) == (row['raw_x_map']['N_homogeneous'], row['raw_x_map']['F_homogeneous']), 'Stored N/F does not replay')
    require(F != 0, 'Exceptional rational-map point')
    g, M = gcd(abs(N), abs(F)), max(abs(N), abs(F))
    require(str(g) == row['finite_cancellation_gcd'], 'Stored cancellation differs')
    bw = word(inverse * vector(ZZ, row['anchor_word_in_D']))
    require(sum((int(c) * p for c, p in zip(bw, basis) if c), E(0)) == q, 'Anchor group word mismatch')
    return {
        'word_in_D': row['anchor_word_in_D'], 'word_in_M17_E': list(bw), 'label': label(bw),
        'point': coords(q), 'chart': chart, 'parameter': {'n': str(n), 'd': str(d), 'H': str(H)},
        'N': str(N), 'F': str(F), 'gcd': str(g), 'max_NF': str(M),
        'reduced_x_height_integer': str(M // g),
        'components_bits': {'k': log2(H), 'quarter_hx': log2(M // g)/4,
                            'a': log2(H) - log2(M)/4, 'c': log2(g)/4},
    }


def prepare():
    save(OUT / 'protocol.json', protocol())
    if (OUT / 'inputs.json.gz').exists() and (OUT / 'events.json').exists():
        return
    raw = gzip.decompress(SOURCE.read_bytes())
    require(sha256(raw).hexdigest() == load(INDEX)['payload_sha256'], 'Original atlas payload hash mismatch')
    atlas = json.loads(raw)
    del raw
    inputs, events = [], []
    for cid, curve in zip(('302', '11952'), atlas['curves']):
        bwords = [x['word_in_D'] for x in curve['M17_specialization'] + curve['exceptional_directions']]
        C = matrix(ZZ, bwords).transpose()
        require(abs(C.det()) == 1, 'Displayed basis is not unimodular')
        inverse = C.inverse().change_ring(ZZ)
        E = EllipticCurve(QQ, curve['curve'])
        basis = [E(list(map(QQ, x['point']))) for x in curve['M17_specialization'] + curve['exceptional_directions']]
        panels = {(p['stage'], p['target_exceptional_index_zero_based']): p for p in curve['nested_target_blind_panels']}
        # Record each old pair only once; all duplicates must agree exactly.
        costs, nested = {}, []
        inverse_cache = {}
        for (stage, j), panel in sorted(panels.items()):
            for bank in ('native_bank', 'generic_only_control'):
                entries = panel[bank]['measurement']['all_anchor_measurements']
                winner = panel[bank]['measurement']['minima']['kappa_parameter_height_bits']
                require(height(winner) == min(map(height, entries)), 'Float-ranked old winner is not an exact height minimum')
                nested.append({'stage': stage, 'target': j, 'bank': bank, 'H': str(height(winner))})
                for row in entries:
                    dw = word(row['anchor_word_in_D'])
                    if dw not in inverse_cache:
                        inverse_cache[dw] = word(inverse * vector(ZZ, dw))
                    bw = inverse_cache[dw]
                    k = (j, bw)
                    value = (row['parameter']['numerator'], row['parameter']['denominator'])
                    if k in costs:
                        require(costs[k] == value, 'Old chart cost depended on surrounding subgroup')
                    costs[k] = value
            if stage == 0:
                continue
            before = panels[(stage-1, j)]['native_bank']['measurement']['minima']['kappa_parameter_height_bits']
            after = panel['native_bank']['measurement']['minima']['kappa_parameter_height_bits']
            if height(before) <= height(after):
                continue
            minus, plus = [event_chart(curve, basis, inverse, row) for row in (before, after)]
            P, R = basis[17+j], basis[17+stage-1]
            for chart in (minus, plus):
                Q = E(list(map(QQ, chart['point'])))
                t = QQ(int(chart['parameter']['n'])) / int(chart['parameter']['d'])
                require(t == (P[1]+Q[1])/(P[0]-Q[0]), 'Slope differs')
                residual = 2*P-Q
                require(QQ(chart['N'])/QQ(chart['F']) == residual[0], 'x(2P-Q) does not equal recorded N/F')
                chart['residual_point'] = coords(residual)
                chart['residual_word_in_M17_E'] = list(combine(tuple(2*v for v in axis(len(basis),17+j)), word(chart['word_in_M17_E']), -1))
            newword = plus['word_in_M17_E']
            require(abs(newword[17+stage-1]) == 1 and sum(abs(v) for v in newword) <= 2, 'Winner is not a new unit/pair with added R')
            require(all(not v or i < 17+stage for i,v in enumerate(newword)), 'New winner uses an unsupplied point')
            hratio = Fraction(int(minus['parameter']['H']), int(plus['parameter']['H']))
            xr = Fraction(int(minus['reduced_x_height_integer']), int(plus['reduced_x_height_integer']))
            cr = Fraction(int(minus['gcd']), int(plus['gcd']))
            ar = hratio**4 * Fraction(int(plus['max_NF']), int(minus['max_NF']))
            require(hratio**4 == xr*ar*cr, 'Exact additive identity fails')
            parts = {k: minus['components_bits'][k]-plus['components_bits'][k] for k in ('quarter_hx','a','c')}
            require(abs(sum(parts.values()) - logratio(hratio)) < 1e-9, 'Displayed event decomposition fails')
            events.append({'curve': cid, 'target': j, 'added': stage-1,
                           'target_point': coords(P), 'target_word_in_D': bwords[17+j],
                           'target_word_in_M17_E': list(axis(len(basis),17+j)),
                           'added_point': coords(R), 'added_word_in_D': bwords[17+stage-1],
                           'added_word_in_M17_E': list(axis(len(basis),17+stage-1)),
                           'before': minus, 'after': plus, 'drop_bits': logratio(hratio),
                           'height_ratio': str(hratio), 'height_ratio_approx': float(hratio),
                           'drop_components_bits': parts,
                           'exact_fourth_power_ratio_factors': {'residual_x':str(xr), 'archimedean':str(ar), 'finite':str(cr)},
                           'exact_identity_checked': True})
        inputs.append({'id':cid, 'model':curve['curve'], 'basis_points':list(map(coords,basis)),
                       'basis_words_in_D':bwords, 'nested_reference':nested,
                       'baseline_cost_cache':[{'target':j, 'word':list(w), 'n':n, 'd':d} for (j,w),(n,d) in sorted(costs.items())]})
    require(len(events) == 4 and all(e['curve']=='302' for e in events), 'Four-event baseline changed')
    save(OUT / 'events.json', {'schema':'rank-accessibility.events.v1', 'protocol_sha256':digest(OUT/'protocol.json'), 'events':events})
    save(OUT / 'inputs.json.gz', {'schema':'rank-accessibility.inputs.v1', 'protocol_sha256':digest(OUT/'protocol.json'), 'curves':inputs})
    print('PREPARE PASS: four exact event identities; old duplicate costs agree.', flush=True)


def generators(q, target, gs, es):
    n = 17+q
    result = [axis(n,i) for i in range(n)]
    for indices, s in ((list(range(17)),gs), ([17+i for i in range(q) if i!=target], es)):
        if s:
            for a,b in zip(indices[::2], indices[1::2]):
                result[a] = combine(result[a], result[b], s)
    C = matrix(ZZ,result).transpose()
    require(C.det()==1, 'Nonunimodular rebase')
    require(all(not result[i][17+target] for i in range(n) if i!=17+target), 'Hidden target leaked into supplied basis')
    return result


def native_universe(gens):
    answer = []
    for i,w in enumerate(gens):
        support = 1 << (i-17) if i>=17 else 0
        answer.extend([(w,support), (negate(w),support)])
    for i,j in combinations(range(len(gens)),2):
        support = (1 << (i-17) if i>=17 else 0) | (1 << (j-17) if j>=17 else 0)
        for s,t in product((1,-1),repeat=2):
            answer.append((tuple(s*a+t*b for a,b in zip(gens[i],gens[j])),support))
    require(len(answer)==2*len(gens)**2 and len({w for w,_ in answer})==len(answer), 'Full universe cardinality mismatch')
    return answer


def generic_stream(gens, count):
    ans=[]
    for weight in (1,2,3):
        for support in combinations(range(17),weight):
            for signs in product((-1,1),repeat=weight-1):
                w=gens[support[0]]
                for i,s in zip(support[1:],signs):
                    w=combine(w,gens[i],s)
                ans.extend((w,negate(w)))
                if len(ans)==count:
                    return ans
    raise ArithmeticError('Generic control bound exceeds frozen triple shell')


class Costs:
    def __init__(self, data):
        self.E=EllipticCurve(QQ,data['model'])
        self.basis=[self.E(list(map(QQ,p))) for p in data['basis_points']]
        self.points={axis(len(self.basis),i):p for i,p in enumerate(self.basis)}
        self.values={(r['target'],word(r['word'])):(int(r['n']),int(r['d'])) for r in data['baseline_cost_cache']}
        self.fresh=0

    def point(self,w):
        if w not in self.points:
            nw=negate(w)
            if nw in self.points:
                self.points[w]=-self.points[nw]
            else:
                self.points[w]=sum((int(c)*p for c,p in zip(w,self.basis) if c),self.E(0))
        return self.points[w]

    def cost(self,j,w):
        require(w[17+j]==0, 'Target-coordinate leak')
        key=(j,w)
        if key not in self.values:
            P,Q=self.basis[17+j],self.point(w)
            require(P[0]!=Q[0], 'Hidden target has exceptional slope')
            t=(P[1]+Q[1])/(P[0]-Q[0])
            self.values[key]=(int(t.numerator()),int(t.denominator()))
            self.fresh+=1
        n,d=self.values[key]
        return max(abs(n),d)


def masks_without(q,j):
    low=(1<<j)-1
    return [(m&low)|((m&~low)<<1) for m in range(1<<(q-1))]


def compact_mask(mask,j):
    return (mask & ((1<<j)-1)) | ((mask>>(j+1))<<j)


def ratio_spectrum(counter):
    """Full discrete survival curve (including zero/negative effects)."""
    remain=sum(counter.values())
    out=[]
    for ratio,count in sorted(counter.items()):
        out.append({'ratio':str(ratio),'bits':logratio(ratio),'multiplicity':count,'A_ge':remain})
        remain-=count
    return out


def spectrum_stats(counter):
    n=sum(counter.values())
    s=sorted(counter.items())
    def at(index):
        c=0
        for v,k in s:
            c+=k
            if c>index:return logratio(v)
    return {'count':n,'positive_count':sum(k for v,k in s if v>1),
            'negative_count':sum(k for v,k in s if v<1),
            'median_bits':(at((n-1)//2)+at(n//2))/2 if n else None,
            'max_bits':logratio(s[-1][0]) if n else None,
            'positive_AUC_bits':sum(logratio(v)*k for v,k in s if v>1),
            'mean_positive_part_bits':sum(logratio(v)*k for v,k in s if v>1)/n if n else None}


def landscape(data,costs,policy):
    name,gs,es=policy
    q=len(costs.basis)-17
    all_edges, all_g=Counter(),Counter()
    records,targets=[],[]
    n_native=2*(16+q)**2
    for j in range(q):
        gens=generators(q,j,gs,es)
        universe=[(w,m) for w,m in native_universe(gens) if not m&(1<<j)]
        require(len(universe)==n_native,'Eligible anchor count')
        controls=generic_stream(gens,n_native)
        # Costs sorted once; minima propagate indices, never float logarithms.
        words=sorted({w for w,_ in universe}|set(controls))
        ids={w:i for i,w in enumerate(words)}
        heights=[costs.cost(j,w) for w in words]
        masks=masks_without(q,j)
        best=[None]*len(masks)
        def better(a,b):
            if a is None:return b
            if b is None:return a
            return a if (heights[a],a)<=(heights[b],b) else b
        for w,support in universe:
            s=compact_mask(support,j)
            best[s]=better(best[s],ids[w])
        for bit in range(q-1):
            for mask in range(len(best)):
                if mask&(1<<bit):
                    best[mask]=better(best[mask],best[mask^(1<<bit)])
        require(all(i is not None for i in best),'Missing empty-support bank')
        control_best=[]
        current=None
        for k in range(q):
            start=0 if k==0 else 2*(16+k)**2
            end=2*(17+k)**2
            for w in controls[start:end]:current=better(current,ids[w])
            control_best.append(current)
        target_edges, target_g=Counter(),Counter()
        effects=Counter()
        edge_examples={}
        for mask,T in enumerate(masks):
            h=heights[best[mask]]
            k=T.bit_count()
            target_g[Fraction(heights[control_best[k]],h)]+=1
            for r in range(q):
                if r==j or T&(1<<r):continue
                after=compact_mask(T|(1<<r),j)
                ratio=Fraction(h,heights[best[after]])
                require(ratio>=1,'Subset minimum increased')
                target_edges[ratio]+=1
                if ratio>1:
                    effects[r]+=1
                    edge_examples.setdefault((r,ratio), {'T':T,'added':r,'before_anchor':best[mask],'after_anchor':best[after],'ratio':str(ratio)})
        # An alternate support scan verifies selected levels independently of DP.
        test_masks=sorted(set([0,len(masks)-1]+[1<<i for i in range(q-1)]+[i for i in range(len(masks)) if (i*104729)%257==0]))
        for m in test_masks:
            T=masks[m]
            direct=min((ids[w] for w,s in universe if s&T==s),key=lambda i:(heights[i],i))
            require(best[m]==direct,'Direct support scan disagrees with subset DP')
        if name=='original':
            for ref in data['nested_reference']:
                if ref['target']!=j:continue
                k=ref['stage']
                idx=best[compact_mask((1<<k)-1,j)] if ref['bank']=='native_bank' else control_best[k]
                require(str(heights[idx])==ref['H'],'Original acquisition path does not replay')
        coefficients=[{'word_in_original_M17_E':list(w), 'n':str(costs.values[j,w][0]), 'd':str(costs.values[j,w][1])} for w in words]
        per_k=[]
        for k in range(q):
            ratios=Counter(Fraction(heights[control_best[k]],heights[best[m]]) for m,T in enumerate(masks) if T.bit_count()==k)
            per_k.append({'supplied_exceptional_count':k,'G':spectrum_stats(ratios)})
        summary={'target':j,'initial_kappa_bits':log2(heights[best[0]]),
                 'full_hidden_kappa_bits':log2(heights[best[-1]]),
                 'best_total_drop_bits':logratio(Fraction(heights[best[0]],heights[best[-1]])),
                 'G_full_hidden_bits':logratio(Fraction(heights[control_best[-1]],heights[best[-1]])),
                 'G':spectrum_stats(target_g),'edges':spectrum_stats(target_edges),
                 'helpful_added_directions':dict(sorted(effects.items())), 'by_subset_size':per_k}
        targets.append(summary)
        records.append({'target':j,'generators_in_original_M17_E':[list(w) for w in gens],
                        'measurements':coefficients,'native_anchor_supports':[[ids[w],s] for w,s in universe],
                        'generic_control_anchor_ids':[ids[w] for w in controls],
                        'subset_masks':masks, 'winning_anchor_ids':best,'generic_control_winners_by_size':control_best,
                        'positive_edge_examples':list(edge_examples.values())})
        all_edges.update(target_edges)
        all_g.update(target_g)
    require(sum(all_g.values())==q*2**(q-1),'Cell count mismatch')
    require(sum(all_edges.values())==q*(q-1)*2**(q-2),'Edge count mismatch')
    summary={'curve':data['id'],'policy':name,'known_exceptional_dimensions':q,
             'cells':sum(all_g.values()),'edges':spectrum_stats(all_edges),'G':spectrum_stats(all_g),
             'targets_with_any_native_gain':sum(t['best_total_drop_bits']>0 for t in targets),
             'targets_with_any_positive_G':sum(t['G']['positive_count']>0 for t in targets),
             'target_mean_G_median_bits':median(t['G']['median_bits'] for t in targets),
             'targets':targets,'edge_survival':ratio_spectrum(all_edges),'G_survival':ratio_spectrum(all_g)}
    return {'schema':'rank-accessibility.subsets.v1','protocol_sha256':digest(OUT/'protocol.json'),
            'input_sha256':digest(OUT/'inputs.json.gz'),'summary':summary,'target_landscapes':records}


def run_landscapes(original_only=False):
    require(load(OUT/'protocol.json')==json.loads(canonical(protocol())), 'Frozen protocol or code changed')
    for data in load(OUT/'inputs.json.gz')['curves']:
        costs=Costs(data)
        for policy in POLICIES[:1] if original_only else POLICIES:
            path=OUT/f"{data['id']}-{policy[0]}.json.gz"
            if path.exists():continue
            payload=landscape(data,costs,policy)
            save(path,payload)
            s=payload['summary']
            print(f"{data['id']} {policy[0]} cells={s['cells']} gain-targets={s['targets_with_any_native_gain']} G-targets={s['targets_with_any_positive_G']} max-edge={s['edges']['max_bits']:.6f} fresh-costs={costs.fresh}",flush=True)
            del payload
            gc.collect()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('stage',choices=['prepare','original','all'])
    args=ap.parse_args()
    if args.stage=='prepare':prepare()
    else:run_landscapes(args.stage=='original')


if __name__=='__main__':main()
