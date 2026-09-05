#!/usr/bin/env python3
"""Frozen 16-cover / six-fibre retrospective test. Exact arithmetic, no searches."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import retrospective as r

PROTOCOL=Path(__file__).resolve().parent/'EXPERIMENT.json'
INPUT=r.OUT/'rank_jump_cover_experiment_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_cover_experiment_v1.json'


def trim(a):
    a=list(a)
    while len(a)>1 and not a[-1]:a.pop()
    return a

def mul(a,b):
    c=[r.F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

def sub(a,b):
    return trim([(a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])

def sqrtq(x):
    if x<0:return None
    a,b=r.isqrt(x.numerator),r.isqrt(x.denominator)
    return r.F(a,b) if a*a==x.numerator and b*b==x.denominator else None

def polynomial_square_root(poly):
    p=trim(poly)
    if p==[0]:return [r.F(0)]
    if len(p)%2==0:return None
    n=(len(p)-1)//2; lead=sqrtq(p[-1])
    if lead is None:return None
    v=[r.F(0)]*(n+1);v[n]=lead
    for k in range(n-1,-1,-1):
        known=sum(v[i]*v[n+k-i] for i in range(k+1,n+1) if k<n+k-i<=n)
        v[k]=(p[n+k]-known)/(2*lead)
    return v if mul(v,v)==p else None

def evaluate(p,t):
    ans=r.F(0)
    for a in reversed(p):ans=ans*t+a
    return ans


def capture():
    protocol=r.read(PROTOCOL)
    raw=(r.ROOT/protocol['atlas']).read_bytes()
    # The large atlas is ignored; the committed deep-cover certificate pins its bytes.
    binding=json.loads(subprocess.check_output(['git','show',r.BASE+':artifacts/generated-results/elliptic-curves/elkies_2026_deep_cover_exceptional_quotients_v1.json'],cwd=r.ROOT))
    assert binding['generation']['inputs'][protocol['atlas']]==protocol['atlas_sha256']
    assert r.digest(raw)==protocol['atlas_sha256']
    atlas=json.loads(raw)
    covers=sorted(atlas['bisections'],key=lambda b:b['label'])[:16]
    assert len({c['label'] for c in covers})==16
    r.write_new(INPUT,{'schema':'rank-jump.cover-experiment-inputs.v1','protocol_sha256':r.digest(PROTOCOL.read_bytes()),
      'atlas_sha256':r.digest(raw),'covers':[{'label':c['label'],'quadratic_cover':c['quadratic_cover']} for c in covers]})


def build(check=False):
    protocol=r.read(PROTOCOL);inp=r.read(INPUT)
    assert inp['protocol_sha256']==r.digest(PROTOCOL.read_bytes())
    ds=[]
    for c in inp['covers']:
        q=c['quadratic_cover'];a,b,e=[list(map(r.F,q[k])) for k in ('leading_coefficients','linear_coefficients','constant_coefficients')]
        ds.append(sub(mul(b,b),[4*x for x in mul(a,e)]))
    # f/g is a square in Q(t) iff f*g is a square; valid for nonzero polynomials.
    groups=[]
    for i,p in enumerate(ds):
        assert p!=[0]
        g=next((g for g in groups if polynomial_square_root(mul(p,ds[g[0]])) is not None),None)
        if g is None:groups.append([i])
        else:g.append(i)
    rows=[]
    for pair in protocol['pairs']:
        for role in ('high','low'):
            t=r.F(pair[role]);vals=[evaluate(p,t) for p in ds];splits=[];degenerate=[];valued_groups=[]
            for i,x in enumerate(vals):
                if not x:degenerate.append(i);continue
                if sqrtq(x) is not None:splits.append(i)
                g=next((g for g in valued_groups if sqrtq(x*vals[g[0]]) is not None),None)
                if g is None:valued_groups.append([i])
                else:g.append(i)
            rows.append({'parameter':str(t),'role':role,'nonzero_rational_split_indices':splits,
              'degenerate_indices':degenerate,'value_squareclass_groups':valued_groups,
              'discriminant_values':list(map(str,vals)),
              'split_roots':{str(i):str(sqrtq(vals[i])) for i in splits}})
    nontrivial=[g for g in groups if len(g)>=3 and polynomial_square_root(ds[g[0]]) is None]
    positive=[{'parameter':row['parameter'],'group':g} for row in rows if row['role']=='high' for g in nontrivial if set(g)<=set(row['nonzero_rational_split_indices'])]
    out={'schema':'rank-jump.cover-experiment.v1','input_sha256':r.digest(INPUT.read_bytes()),
      'script_sha256':r.digest(Path(__file__).read_bytes()),'protocol_sha256':r.digest(PROTOCOL.read_bytes()),
      'discriminant_polynomials_ascending':[list(map(str,d)) for d in ds],
      'generic_squareclass_groups':groups,'eligible_shared_nontrivial_carriers':nontrivial,'specializations':rows,
      'primary_positive_witnesses':positive,'primary_endpoint':'SIMULTANEOUS_LIFT_GATE_REQUIRES_INDEPENDENCE' if positive else 'FIXED_DICTIONARY_SHARED_CARRIER_HYPOTHESIS_REFUTED',
      'scope':'Only this frozen dictionary. No negative rank conclusion; a positive split need not add a quotient direction.'}
    if check:
        if r.read(OUTPUT)!=out:raise ValueError('cover experiment mismatch')
        print('PASS cover replay')
    else:r.write_new(OUTPUT,out)
    print(out['primary_endpoint']);print('groups',groups)
    for row in rows:print(row['role'],row['parameter'],'split',row['nonzero_rational_split_indices'],'degenerate',row['degenerate_indices'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('capture','build','check'));a=p.parse_args()
    capture() if a.mode=='capture' else build(a.mode=='check')
