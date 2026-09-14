#!/usr/bin/env python3
"""Complete Q80 genus0-k0 replay, including an independent exact norm4 count."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import isqrt
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
FINITE='elkies-k3/scripts/verify_q80_genus_zero_collision_closure.py'
def require(v,msg):
    if not v:raise ValueError(msg)
def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def exact_short_count(G,bound):
    """Exact LDL traversal; integer square-root intervals cannot miss a vector."""
    n=len(G);require(n and all(len(r)==n for r in G),'square Gram matrix')
    require(all(G[i][j]==G[j][i] for i in range(n) for j in range(n)),'symmetric Gram matrix')
    L=[[Q(i==j) for j in range(n)] for i in range(n)];D=[]
    for i in range(n):
        d=Q(G[i][i])-sum(L[i][k]**2*D[k] for k in range(i))
        require(d>0,'positive definite Gram matrix');D.append(d)
        for j in range(i+1,n):L[j][i]=(G[j][i]-sum(L[j][k]*L[i][k]*D[k] for k in range(i)))/d
    w=[0]*n;hist=Counter();nodes=0
    def visit(i,remaining):
        nonlocal nodes
        nodes+=1
        if i<0:
            norm=sum(w[a]*G[a][b]*w[b] for a in range(n) for b in range(n))
            require(norm==bound-remaining and 0<=norm<=bound,'exact leaf norm')
            hist[norm]+=1;return
        centre=sum(L[j][i]*w[j] for j in range(i+1,n))
        a,b=centre.numerator,centre.denominator
        radius=remaining/D[i]*b*b
        N=isqrt(radius.numerator//radius.denominator)
        hi=(N-a)//b
        # ceil((-N-a)/b), including negative numerators.
        lo=-((N+a)//b)
        for value in range(lo,hi+1):
            term=D[i]*(Q(value)+centre)**2
            require(term<=remaining,'exact traversal interval')
            w[i]=value;visit(i-1,remaining-term)
        w[i]=0
    visit(n-1,Q(bound))
    return {'signed_norm_histogram':{str(k):v for k,v in sorted(hist.items())},'visited_nodes':nodes}
def freeze(out):
    packet={'schema':1,'source':SOURCE,'source_sha256':digest(ROOT/SOURCE),
        'checker_sha256':digest(Path(__file__)),'finite_checker':FINITE,'finite_checker_sha256':digest(ROOT/FINITE),
        'finite_input_sha256':digest(out/'input.json'),'bound':4,'rank':17,
        'cpu_limit_seconds':40,'memory_limit_bytes':4*1024**3,
        'scope':'Independent exact rational-LDL norm4 count plus the complete finite contact replay. No norm10 census or positive high-rank family.'}
    with (out/'complete-replay-input.json').open('x') as f:json.dump(packet,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({'status':'FROZEN','input_sha256':digest(out/'complete-replay-input.json')}),flush=True)
def verify(out):
    started=time.monotonic();packet=read(out/'complete-replay-input.json')
    require(packet['source']==SOURCE and packet['source_sha256']==digest(ROOT/SOURCE),'literal Gram source')
    require(packet['checker_sha256']==digest(Path(__file__)),'count checker binding')
    require(packet['finite_checker']==FINITE and packet['finite_checker_sha256']==digest(ROOT/FINITE),'finite checker binding')
    require(packet['finite_input_sha256']==digest(out/'input.json'),'finite input binding')
    require(packet['bound']==4 and packet['rank']==17 and packet['cpu_limit_seconds']==40 and packet['memory_limit_bytes']==4*1024**3,'frozen exact-count scope')
    source=read(ROOT/SOURCE);G=source['sections']['height_gram'];require(len(G)==17,'source rank')
    count=exact_short_count(G,4)
    require(count['signed_norm_histogram']=={'0':1,'4':2626},'complete exact norm4 census')
    spec=importlib.util.spec_from_file_location('q80_complete_finite_closure',ROOT/FINITE)
    finite=importlib.util.module_from_spec(spec);spec.loader.exec_module(finite)
    replay=finite.verify(out)
    require(replay['norm4_sections']*2==count['signed_norm_histogram']['4'],'independent count and point roster agree')
    return {'status':'PASS','input_sha256':digest(out/'complete-replay-input.json'),
        'finite_input_sha256':digest(out/'input.json'),'finite_result_sha256':digest(out/'result.json'),
        'independent_exact_count':count,'finite_replay':replay,
        'positive_mw17_target_complete':False,'old_norm8_singular_replay_upgraded':False,
        'elapsed_seconds':round(time.monotonic()-started,6)}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--freeze',action='store_true');parser.add_argument('--record',type=Path);args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    if args.freeze:freeze(args.input)
    else:
        result=verify(args.input)
        if args.record:
            with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
        print(json.dumps(result,sort_keys=True),flush=True)
