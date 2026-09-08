"""Bounded retrospective M18-only features. Labels are a separate, sealed step.

No chart transcripts or M19+ basis are feature inputs. Exactness concerns the
retained integer decision metric, not the numerical canonical heights.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from visibility_lattice_v2 import ExactParity

ROOT = Path(__file__).resolve().parents[3]
LOCAL = ROOT/'research/artifacts/local/elliptic-curves'
OUT = ROOT/'research/artifacts/generated-results/elliptic-curves/m18_landscape_comparison_v1'
SNAP = LOCAL/'m18-landscape-comparison-v1'
THRESHOLDS = ['1/2', '1', '3/2', '2', '5/2', '3', '4']
PATTERNS = ['curve302-seeded-v3-amplifier-v1/*/seed-input.json',
 'curve302-seed-universality-panel-v1/*/seed-input.json',
 'det1092-funnel-conic-split-seeds-v1/amplifiers/*/seed-input.json',
 'det1092-small-conic-seed-v1/amplifiers/*/seed-input.json',
 'det1092-funnel-production-v1/amplifiers/*/seed-input.json',
 'orbit8044-seed-factory-pilot-v1/amplifiers/*/seed-input.json']
POLICY_KEYS = ['anchors_per_shell','canonical_per_shell','exact_cvp_node_limit',
 'height','seconds_per_chart','max_charts','max_epochs','target_rank','metric']

def require(ok, message):
    if not ok: raise ValueError(message)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True)+'\n'
    if p.exists(): require(p.read_text()==text, f'immutable output changed: {p}')
    else:
        with p.open('x') as f: f.write(text)

def bits(x):
    x=Q(x)
    return abs(x.numerator).bit_length()+x.denominator.bit_length()

def quantile(xs, numerator, denominator):
    xs=sorted(xs)
    return xs[max(0, math.ceil(len(xs)*numerator/denominator)-1)]

def stats(name, xs, result):
    require(bool(xs), 'empty feature population')
    for tag,n,d in [('min',0,1),('q10',1,10),('q25',1,4),('median',1,2),('q75',3,4),('q90',9,10),('max',1,1)]:
        result[name+'_'+tag]=float(quantile(xs,n,d))

def entropy(counts):
    total=sum(counts)
    return -sum((n/total)*math.log2(n/total) for n in counts if n) if total else None

def freeze():
    require(not (OUT/'protocol.json').exists(), 'protocol already frozen; use extract')
    cases=[]
    for pattern in PATTERNS:
        for seed in sorted(LOCAL.glob(pattern)):
            folder=seed.parent
            selection=folder/'replay-M17/epoch-00/selection.json'
            if not selection.exists(): continue
            case=folder.parent.name+'__'+folder.name
            if folder.parent.name=='amplifiers': case=folder.parent.parent.name+'__'+folder.name
            dest=SNAP/case
            srcs=[seed,folder/'seed-proof.json',selection]+sorted(selection.parent.glob('anchor-*-full.npz'))
            inputs={}
            for src in srcs:
                target=dest/src.name
                target.parent.mkdir(parents=True,exist_ok=True)
                if not target.exists(): shutil.copyfile(src,target)
                require(sha(src)==sha(target),'snapshot mismatch')
                inputs[src.name]={'source':str(src.relative_to(ROOT)), 'sha256':sha(src)}
            policy={k:read(folder/'protocol.json')[k] for k in POLICY_KEYS}
            write(dest/'policy.json',policy)
            inputs['policy.json']={'source':str((dest/'policy.json').relative_to(ROOT)), 'sha256':sha(dest/'policy.json')}
            cases.append({'case':case,'folder':str(folder.relative_to(ROOT)),'inputs':inputs})
    write(OUT/'protocol.json', {'schema':'m18-only-landscape-v1','created_unix':time.time(),
      'roster_rule':'All sealed M18 packets with an initial epoch-00 landscape in the six listed campaigns at freeze; no outcome condition.',
      'patterns':PATTERNS,'cases':cases,'thresholds_median_M17_diagonal':THRESHOLDS,
      'competitive_rule':'CVP norm <= 5/4 of minimum, inclusive; anchor entropy in bits',
      'quantiles':'nearest rank, minimum for q0; 128th absent is null',
      'quartics':'complete refined parity population only; no executed-chart sampling',
      'bounds':'64 exact CVPs per case; <=2000000 nodes each; no point search',
      'claims':'Exact CVP in retained rounded integer Gram only. Retrospective diagnostic, not predictive validation.',
      'source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),Path(__file__).with_name('visibility_lattice_v2.py')]}})
    print('FROZEN',len(cases),'M18 states',flush=True)

def features(seed, selection, folder, replay=True):
    require(seed['initial_rank']==18 and seed['generic_rank']==17 and len(seed['points'])==18,'not an M18 packet')
    require(selection['rank']==18 and selection['basis']==seed['points'],'M18 feature boundary violated')
    curve=list(map(Q,seed['curve']))
    require(curve[:3]==[0,0,0], 'short model required')
    a,b=curve[3:]
    for raw in seed['points']:
        x,y=map(Q,raw);require(y*y==x*x*x+a*x+b,'seed point off curve')
    j=Q(1728)*4*a**3/(4*a**3+27*b**2)
    g=np.array(selection['rounded_gram'],dtype=object)
    u=np.array(selection['LLL'],dtype=object)
    require(g.shape==u.shape==(18,18),'wrong Gram dimension')
    require(np.array_equal(g,g.T),'asymmetric Gram')
    scale=Q(sorted(int(g[i,i]) for i in range(17))[8])
    original=ExactParity(g.tolist())
    exact=ExactParity((u@g@u.T).tolist())
    # Fraction elimination avoids a floating determinant or inverse check.
    work=[[Q(x) for x in row] for row in u.tolist()]; det=Q(1)
    for k in range(18):
        pivot=next((i for i in range(k,18) if work[i][k]),None)
        require(pivot is not None,'singular LLL')
        if pivot!=k:work[k],work[pivot]=work[pivot],work[k];det=-det
        v=work[k][k];det*=v
        for i in range(k+1,18):
            r=work[i][k]/v
            for z in range(k+1,18):work[i][z]-=r*work[k][z]
    require(abs(det)==1,'non-unimodular LLL')
    rows=[]
    require(len(selection['anchors'])==32 and selection['extensions_per_anchor']==2 and selection['full_cosets_scored']==64,'incomparable M18 policy')
    for ai,record in enumerate(selection['anchors']):
        path=folder/f'anchor-{ai:02d}-full.npz'
        require(sha(path)==record['full_scores_sha256'],'NPZ changed')
        with np.load(path,allow_pickle=False) as data:
            residues=data['residues'];rp=data['reduced_residues'];words=data['babai_words']
            require(residues.shape==(2,18) and set(residues[:,17])=={0,1},'extension coverage')
            require(np.array_equal((rp@u)%2,residues),'parity transport')
            require(np.all(residues[:,:17]==np.array(record['anchor']['representative'][:17])%2),'anchor parity')
            bw,bn=exact.babai(rp)
            require(np.array_equal(bw,words) and np.array_equal(bn,data['babai_norms']),'Babai replay')
            require({r['extension'] for r in record['refined']}=={0,1} and len(record['refined'])==2,'exact mask omitted')
            for row in record['refined']:
                i=row['extension'];cert=row['cvp']
                if replay:
                    actual=exact.solve(rp[i],words[i],2000000)
                    require(json.dumps(actual,sort_keys=True)==json.dumps(cert,sort_keys=True),'CVP replay failed')
                w=np.array(row['representative'],dtype=object)
                require(int(w@g@w)==row['metric_norm']==cert['norm'],'CVP norm mismatch')
                require(np.array_equal(w%2,residues[i]),'representative parity')
                require(row['multiplicity']==len(cert['minima'])//2,'multiplicity mismatch')
                rows.append({'anchor':ai,'orbit':record['anchor']['orbit'],'shell':row['shell'],
                    'extension':i,'norm':row['metric_norm'],'normalized_norm':str(Q(row['metric_norm'])/scale),
                    'multiplicity':row['multiplicity'],'quartic_bits':row['quartic_bits'],
                    'quartic_max_bits':row['quartic_max_bits'],
                    'coordinate_bits':sum(bits(x) for x in row['point']),
                    'babai_ratio':str(Q(row['babai_norm'],row['metric_norm']))})
    norms=sorted(Q(r['norm'])/scale for r in rows)
    good=[r for r in rows if Q(r['norm'])/scale<=norms[0]*Q(5,4)]
    counts=Counter(r['anchor'] for r in good)
    f={'cosets':len(rows),'anchors':32,'centres':len(selection['centres']),
       'scale':int(scale),'j_bits':bits(j), 'seed_height_ratio':float(Q(int(g[17,17]))/scale),
       'seed_schur_ratio':float(original.d[-1]/scale),
       'competitive_masks':len(good),'competitive_anchors':len(counts),
       'competitive_anchor_entropy':entropy(counts.values()),
       'competitive_anchor_max_share':max(counts.values())/len(good),
       'competitive_seed_bit_entropy':entropy(Counter(r['extension'] for r in good).values()),
       'gap_8':float(norms[7]-norms[0]),'gap_32':float(norms[31]-norms[0]),'gap_128':None,
       'tail_ratio':float(norms[-1]/norms[0]),
       'multiple_minimum_masks':sum(r['multiplicity']>1 for r in rows)}
    stats('cvp',norms,f);stats('raw_cvp',[r['norm'] for r in rows],f)
    for k in ['quartic_bits','quartic_max_bits','coordinate_bits','multiplicity','babai_ratio']:
        stats(k,[Q(r[k]) for r in rows],f)
    for shell in [8,10]:
        rs=[r for r in rows if r['shell']==shell]
        f[f'shell{shell}_minima_pairs']=sum(r['multiplicity'] for r in rs)
        stats(f'shell{shell}_cvp',[Q(r['norm'])/scale for r in rs],f)
    for ext in [0,1]:stats(f'extension{ext}_cvp',[Q(r['norm'])/scale for r in rows if r['extension']==ext],f)
    for t in THRESHOLDS:f['masks_le_'+t]=sum(n<=Q(t) for n in norms)
    return {'j':str(j),'parameter':seed.get('parameter'), 'features':f,'parity_rows':rows,
            'verification':'all 64 retained exact CVPs replayed' if replay else 'NOT_REPLAYED'}

def extract():
    protocol=read(OUT/'protocol.json')
    for path,digest in protocol['source_hashes'].items():require(sha(ROOT/path)==digest,'feature source changed after protocol freeze')
    seals={}
    for case in protocol['cases']:
        name=case['case'];folder=SNAP/name
        for filename,record in case['inputs'].items():require(sha(folder/filename)==record['sha256'],'feature input hash changed')
        target=OUT/'features'/f'{name}.json'
        if not target.exists():
            require(read(folder/'seed-proof.json')['rank_lower_bound']==18,'uncertified M18 packet')
            result=features(read(folder/'seed-input.json'),read(folder/'selection.json'),folder)
            result.update(case=name, protocol_sha256=sha(OUT/'protocol.json'))
            write(target,result)
        seals[str(target.relative_to(OUT))]=sha(target)
        print('FEATURES',name,flush=True)
    write(OUT/'feature-seal.json',{'protocol_sha256':sha(OUT/'protocol.json'),'features':seals,
         'completed_unix':time.time(),'boundary':'No outcome files read by freeze/extract. Labels may now be joined.'})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['freeze','extract']);args=parser.parse_args()
    {'freeze':freeze,'extract':extract}[args.mode]()
