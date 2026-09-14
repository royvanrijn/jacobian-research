#!/usr/bin/env python3
"""Independent integer replay: count every (t,x,y) and check good reduction.

The Frobenius eigenvalue bound and specialization injectivity are written
mathematical inputs. This is not a full zeta-polynomial or exact-rank checker.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'artifacts/generated-results/elkies-k3-common-quartic-control-rank-gate-v1'
SOURCE = 'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/input.json'


def require(condition,message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def add(a,b):
    c = [0]*max(len(a),len(b))
    for i,v in enumerate(a): c[i] += v
    for i,v in enumerate(b): c[i] += v
    return c


def mul(a,b):
    c = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j] += x*y
    return c


def trim(a):
    while a and not a[-1]: a.pop()
    return a


def rem(a,b,p):
    a = trim(a[:])
    while len(a)>=len(b):
        factor,shift = a[-1]*pow(b[-1],-1,p) % p,len(a)-len(b)
        for i,c in enumerate(b): a[shift+i] = (a[shift+i]-factor*c) % p
        trim(a)
    return a


def squarefree_degree24(d,p):
    require(len(d)==25 and d[-1] % p, 'degree24 good infinity required')
    a = [v % p for v in d]
    b = trim([(i*d[i]) % p for i in range(1,25)])
    while b: a,b = b,rem(a,b,p)
    require(len(a)==1,'non-squarefree discriminant: no smooth-model certificate')


def ev(c,t,p):
    return sum(v*pow(t,i,p) for i,v in enumerate(c)) % p


def check_control(control,record,p):
    require(record['name']==control['name'] and record['prime']==p,'control attachment')
    A = add(mul(control['D'],add([2*c for c in control['s']],[1])),[-1])
    B = mul(control['D'],mul(control['s'],control['s']))
    require(record['A']==A and record['B']==B,'literal parent identity')
    d = [v % p for v in add([4*v for v in mul(A,mul(A,A))],[27*v for v in mul(B,B)])]
    require(d==record['discriminant_mod_p'],'discriminant identity')
    squarefree_degree24(d,p)
    counts = []
    for t in list(range(p))+[None]:
        a,b = (A[8] % p,B[12] % p) if t is None else (ev(A,t,p),ev(B,t,p))
        count = 1
        for x in range(p):
            for y in range(p):
                count += (y*y-x**3-a*x-b) % p == 0
        counts.append(count)
    require(counts==record['fibre_counts_finite_then_infinity'],'complete fibre counts including infinity')
    N = sum(counts)
    trace = N-1-p*p
    rho = (22*p+trace)//(2*p)
    require((N,trace,rho,rho-2)==(record['surface_count'],record['h2_trace'],record['rational_picard_rank_upper'],record['any_rational_fibration_mw_rank_upper']),'trace or rank bound')
    require(record['exact_inherited_rank']==record['geometric_picard_rank']=='UNKNOWN','upper bound promoted to exact or geometric rank')
    return {'name':control['name'],'prime':p,'surface_count':N,
            'rational_picard_rank_upper':rho,'any_rational_fibration_mw_rank_upper':rho-2}


def verify(out):
    start = time.monotonic()
    packet,result = read(out/'input.json'),read(out/'result.json')
    require(packet['source']==SOURCE and packet['source_sha256']==digest(ROOT/SOURCE),'source binding')
    require(packet['controls']==read(ROOT/SOURCE)['controls'],'control source attachment')
    require(packet['primes']==[5,7] and len(packet['controls'])==2,'scope')
    require(packet['producer_sha256']==digest(ROOT/'elkies-k3/scripts/certify_common_quartic_control_rank_gate.sage'),'producer binding')
    require(packet['checker_sha256']==digest(Path(__file__)),'checker binding')
    require(packet['preview_sha256']==digest(out/'preview.json'),'preview retained')
    require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','result binding')
    require(len(result['records'])==2,'missing or extra control')
    rows=[]
    for c,p,m in zip(packet['controls'],packet['primes'],result['records']):
        require(m['path']==c['name']+'-checkpoint.json' and digest(out/m['path'])==m['sha256'],'checkpoint binding')
        rows.append(check_control(c,read(out/m['path']),p))
    require(all(r['any_rational_fibration_mw_rank_upper']<16 for r in rows),'claimed MW16 exclusion missing')
    require(result['any_mw16_or_mw17_fibration_on_either_control']=='EXCLUDED' and result['positive_mw17_target_complete'] is False,'unsupported endpoint')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
            'controls':rows,'point_counts_independently_replayed':True,
            'written_cohomology_formally_verified':False,'positive_mw17_target_complete':False,
            'elapsed_seconds':round(time.monotonic()-start,6)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--record',type=Path)
    args=parser.parse_args(); result=verify(args.input)
    if args.record:
        with args.record.open('x') as f: json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)
