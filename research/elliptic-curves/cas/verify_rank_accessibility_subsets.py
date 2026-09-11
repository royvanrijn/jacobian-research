#!/usr/bin/env python3
"""Independent Fraction group law, slope replay, and sparse support envelopes.

This verifier does not import the Sage producer or use its subset DP. Every
stored rational slope is checked, and every cell is obtained by scanning a
minimal list of support offers. No new points are sought.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import gzip
from itertools import combinations, product
import json
from math import gcd, log2
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elliptic-curves/rank_accessibility_subsets_v1'


def require(value,message):
    if not value:raise ArithmeticError(message)


def load(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)


def digest(path):return sha256(path.read_bytes()).hexdigest()


def add(P,R,A):
    if P is None:return R
    if R is None:return P
    x,y=P;u,v=R
    if x==u and y==-v:return None
    slope=(3*x*x+A)/(2*y) if P==R else (v-y)/(u-x)
    z=slope*slope-x-u
    return z,-y+slope*(x-z)


def scale(P,n,A):
    if n<0:return scale((P[0],-P[1]),-n,A)
    result=None
    while n:
        if n&1:result=add(result,P,A)
        P=add(P,P,A);n>>=1
    return result


class Points:
    def __init__(self,data):
        self.A,self.B=map(Q,data['model'][3:])
        self.basis=[tuple(map(Q,p)) for p in data['basis_points']]
        for x,y in self.basis:require(y*y==x*x*x+self.A*x+self.B,'Off-curve basis point')
        self.cache={}
    def get(self,w):
        if w not in self.cache:
            nw=tuple(-v for v in w)
            if nw in self.cache:
                x,y=self.cache[nw];self.cache[w]=(x,-y)
            else:
                P=None
                for c,R in zip(w,self.basis):
                    if c:P=add(P,scale(R,c,self.A),self.A)
                require(P is not None,'Zero chart anchor')
                self.cache[w]=P
        return self.cache[w]


def word_label(w):
    ans=[]
    for i,c in enumerate(w):
        if c:ans.append(('+' if c>0 else '-')+(str(abs(c))+'*' if abs(c)!=1 else '')+(f'M{i+1}' if i<17 else f'E{i-16}'))
    return ''.join(ans).lstrip('+')


def verify_events(inputs):
    events=load(OUT/'events.json')['events']
    for event in events:
        data=inputs[event['curve']];points=Points(data)
        P=tuple(map(Q,event['target_point']))
        require(P==points.get(tuple(event['target_word_in_M17_E'])),'Event target word')
        require(tuple(map(Q,event['added_point']))==points.get(tuple(event['added_word_in_M17_E'])),'Added point word')
        for side in ('before','after'):
            row=event[side];anchor=points.get(tuple(row['word_in_M17_E']))
            require(anchor==tuple(map(Q,row['point'])),'Event anchor point')
            residual=add(scale(P,2,points.A),(anchor[0],-anchor[1]),points.A)
            require(residual==tuple(map(Q,row['residual_point'])),'Event residual point')
            n,d=int(row['parameter']['n']),int(row['parameter']['d'])
            t=Q(n,d)
            require(t==(P[1]+anchor[1])/(P[0]-anchor[0]),'Event slope')
            a,b=anchor
            rawF=[-3*a*a-4*points.A,-8*b,-6*a,0,1]
            rawN=[a**3+4*points.B,4*a*b,6*a*a+4*points.A,4*b,a]
            stored=row['chart']['primitive_integral_x_map']
            nf=list(map(int,stored['N_ascending']))+list(map(int,stored['F_ascending']))
            raw=rawN+rawF
            factor=Q(nf[-1])/Q(raw[-1])
            require(all(Q(i)==factor*r for i,r in zip(nf,raw)),'Primitive N/F polynomial mismatch')
            require(gcd(*nf)==1,'Nonprimitive N/F pair')
            N=sum(nf[i]*n**i*d**(4-i) for i in range(5))
            F=sum(nf[5+i]*n**i*d**(4-i) for i in range(5))
            require((N,F)==(int(row['N']),int(row['F'])),'Evaluated N/F differs')
            require(F!=0 and Q(N,F)==residual[0],'Group residual N/F identity')
            require(gcd(abs(N),abs(F))==int(row['gcd']),'Event gcd')
            require(int(row['reduced_x_height_integer'])==max(abs(residual[0].numerator),residual[0].denominator),'Residual rational height')
        factors=event['exact_fourth_power_ratio_factors']
        require(Q(event['height_ratio'])**4==Q(factors['residual_x'])*Q(factors['archimedean'])*Q(factors['finite']),'Event multiplicative identity')
    return len(events)


def expected_generators(n,j,gs,es):
    # Independently reconstruct the frozen disjoint shear matrices.
    gens=[tuple(int(i==k) for k in range(n)) for i in range(n)]
    for indices,s in ((list(range(17)),gs),([i for i in range(17,n) if i!=17+j],es)):
        if s:
            for k in range(0,len(indices)-1,2):
                a,b=indices[k:k+2]
                gens[a]=tuple(int(i==a)+s*int(i==b) for i in range(n))
    return gens


def verify_landscape(path,data,points,gs,es,checked_costs):
    payload=load(path);summ=payload['summary'];q=summ['known_exceptional_dimensions']
    require(payload['protocol_sha256']==digest(OUT/'protocol.json'),'Policy checkpoint hash')
    require(payload['input_sha256']==digest(OUT/'inputs.json.gz'),'Input checkpoint hash')
    cells=0;slopes=0;envelopes=[]
    edge_counter=Counter();g_counter=Counter()
    for record in payload['target_landscapes']:
        j=record['target'];P=points.basis[17+j]
        gens=expected_generators(17+q,j,gs,es)
        require(list(map(list,gens))==record['generators_in_original_M17_E'],'Frozen basis transformation differs')
        ids={tuple(r['word_in_original_M17_E']):i for i,r in enumerate(record['measurements'])}
        hs=[]
        for r in record['measurements']:
            w=tuple(r['word_in_original_M17_E']);n,d=int(r['n']),int(r['d'])
            require(d>0 and gcd(n,d)==1 and w[17+j]==0,'Nonprimitive slope or hidden-target contamination')
            key=(j,w,n,d)
            if key not in checked_costs:
                a,b=points.get(w)
                require((P[1]+b)*d==(P[0]-a)*n and P[0]!=a,'Exact Fraction slope replay failed')
                checked_costs.add(key);slopes+=1
            hs.append(max(abs(n),d))
        native={}
        supplied=[i for i in range(17+q) if i!=17+j]
        for size in (1,2):
            for support in combinations(supplied,size):
                for signs in product((-1,1),repeat=size):
                    w=tuple(sum(s*gens[i][k] for s,i in zip(signs,support)) for k in range(17+q))
                    mask=sum(1<<(i-17) for i in support if i>=17)
                    native[ids[w]]=mask
        require(native==dict(record['native_anchor_supports']),'Native anchor universe or exceptional support differs')
        require(len(native)==2*(16+q)**2,'Native universe count')
        # Verify the generic control order, including every truncated shell.
        generic=[]
        for weight in (1,2,3):
            for support in combinations(range(17),weight):
                for suffix in product((-1,1),repeat=weight-1):
                    w=tuple(sum(s*gens[i][k] for s,i in zip((1,)+suffix,support)) for k in range(17+q))
                    generic.extend((w,tuple(-x for x in w)))
                    if len(generic)==2*(16+q)**2:break
                if len(generic)==2*(16+q)**2:break
            if len(generic)==2*(16+q)**2:break
        require([ids[w] for w in generic]==record['generic_control_anchor_ids'],'Generic-only control policy differs')
        control=[]
        for k in range(q):
            selected=record['generic_control_anchor_ids'][:2*(17+k)**2]
            c=min(selected,key=lambda i:(hs[i],i));control.append(c)
        require(control==record['generic_control_winners_by_size'],'Generic prefix minimum')
        by_support={}
        for i,s in native.items():
            if s not in by_support or (hs[i],i)<(hs[by_support[s]],by_support[s]):by_support[s]=i
        # Drop a support only if a proper sub-support is at least as cheap.
        # This envelope is a direct logical predicate, independent of DP order.
        offers={s:i for s,i in by_support.items() if not any(t!=s and t&s==t and hs[k]<=hs[i] for t,k in by_support.items())}
        masks=record['subset_masks']
        require(masks==[m for m in range(1<<q) if not m&(1<<j)],'Held-out subset census')
        winners=record['winning_anchor_ids'];by_mask=dict(zip(masks,winners))
        for T,idx in zip(masks,winners):
            require(idx in native and native[idx]&T==native[idx],'Winning anchor not supplied')
            direct=min(hs[i] for s,i in offers.items() if s&T==s)
            require(hs[idx]==direct,'Full-cell sparse envelope differs from DP')
            g_counter[Q(hs[control[T.bit_count()]],direct)]+=1
            for r in range(q):
                if r==j or T&(1<<r):continue
                edge_counter[Q(direct,hs[by_mask[T|(1<<r)]])]+=1
        h0=hs[by_support[0]]
        envelopes.append({'target':j,'initial_H':str(h0),'initial_kappa_bits':log2(h0),
                          'offers':[{'exceptional_support_mask':s,'required_exceptional_indices':[r for r in range(q) if s&(1<<r)],
                                     'anchor_word_in_original_M17_E':record['measurements'][i]['word_in_original_M17_E'],
                                     'anchor_label':word_label(record['measurements'][i]['word_in_original_M17_E']),
                                     'H':str(hs[i]),'improvement_over_initial_bits':log2(h0)-log2(hs[i]),
                                     'parameter':{'n':record['measurements'][i]['n'],'d':record['measurements'][i]['d']}}
                                    for s,i in sorted(offers.items()) if s]})
        cells+=len(masks)
    for name,counter in (('edge_survival',edge_counter),('G_survival',g_counter)):
        spectrum=summ[name]
        require({Q(r['ratio']):r['multiplicity'] for r in spectrum}==counter,'Survival multiplicities')
        remain=sum(counter.values())
        for row in spectrum:
            require(row['A_ge']==remain,'Survival cumulative count')
            remain-=row['multiplicity']
    require(cells==q*2**(q-1),'Complete held-out cell count')
    return {'file':path.name,'sha256':digest(path),'cells_checked':cells,'new_rational_slopes_checked':slopes,
            'edges_checked':sum(edge_counter.values()),'minimal_support_envelopes':envelopes}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
    protocol=load(OUT/'protocol.json')
    for relative,expected in protocol['inputs'].items():require(digest(ROOT/relative)==expected,'Protocol source hash mismatch: '+relative)
    inputs={r['id']:r for r in load(OUT/'inputs.json.gz')['curves']}
    count=verify_events(inputs);results=[]
    for cid,data in inputs.items():
        points=Points(data);checked=set()
        for name,gs,es in protocol['policies']:
            path=OUT/f'{cid}-{name}.json.gz'
            result=verify_landscape(path,data,points,gs,es,checked);results.append(result)
            print(f"VERIFIED {path.name}: {result['cells_checked']} cells; {result['new_rational_slopes_checked']} newly replayed slopes",flush=True)
    output={'schema':'rank-accessibility.subsets-verification.v1','status':'PASS',
            'verifier_sha256':digest(Path(__file__)),'protocol_sha256':digest(OUT/'protocol.json'),
            'four_event_group_and_exact_height_identities_checked':count,
            'independence':'Python Fraction elliptic group law and slope checks; direct sparse support envelopes for every subset cell; complete edge/control survival counts.',
            'results':results}
    raw=(json.dumps(output,sort_keys=True,indent=2)+'\n').encode();path=OUT/'verified.json'
    if args.check:require(path.read_bytes()==raw,'Verification replay differs')
    else:path.write_bytes(raw)
    print('RANK_ACCESSIBILITY_SUBSETS_VERIFY PASS',flush=True)


if __name__=='__main__':main()
