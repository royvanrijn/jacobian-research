#!/usr/bin/env sage-python
"""Prospective same-cover comparison of the retained complete Q80 norm8 layer.

The worker reads a generic-only frozen projection. It never reads the previous
product targets, exceptional points, specialization parameters or point ranks.
Each prime stage is checkpointed; absent stages never count as exclusions.
"""
import argparse
from functools import lru_cache
import csv
import gzip
from hashlib import sha256
import itertools
import json
from pathlib import Path
import resource
import runpy
import time

import numpy as np
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, gcd, lcm, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1'
DIRECT=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
PRIORITY=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json'
TABLE=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv'
HELPER=ROOT/'elkies-k3/scripts/construct_elkies_2026_bisections.sage'
INVERSION=ROOT/'elkies-k3/scripts/search_r17_norm12_11952_product_bisection_inversion.sage'


def write_new(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n').encode()
    if path.suffix=='.gz':raw=gzip.compress(raw,mtime=0)
    with path.open('xb') as f:f.write(raw)


def read(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)


def freeze(path):
    direct=read(DIRECT);priority=read(PRIORITY)
    assert sha256(TABLE.read_bytes()).hexdigest()==priority['priority_table_sha256']
    assert sha256(DIRECT.read_bytes()).hexdigest()==priority['inputs'][str(DIRECT.relative_to(ROOT))]
    with TABLE.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    words=[[int(x) for x in r['section_basis_w'].split()] for r in rows]
    gram=[[int(QQ(x)) for x in row] for row in direct['sections']['height_gram']]
    assert matrix(QQ,gram)==matrix(QQ,direct['sections']['height_gram'])
    W=np.asarray(words,dtype=np.int64);G=np.asarray(gram,dtype=np.int64)
    assert len(rows)==63917 and np.all(np.sum((W@G)*W,axis=1)==8)
    assert len({tuple(x%2 for x in w) for w in words})==63917
    packet={'schema':'q80-complete-genus-one-pairs-input-v1',
            'A':direct['weierstrass_model']['A_coefficients_low_to_high'],
            'B':direct['weierstrass_model']['B_coefficients_low_to_high'],
            'basis':[{'X':r['X'],'Y':r['Y']} for r in direct['sections']['records']],
            'gram':gram,'words':words,'orbit_masks':[int(r['orbit_mask']) for r in rows],
            'inventory_certificate':priority,
            'source_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [DIRECT,PRIORITY,TABLE,HELPER,INVERSION]},
            'primes':[509,521,523,541],
            'limits':{'cpu_seconds_per_prime':240,'address_space_gib':4,'checkpoint_size':512},
            'scope':'All63917 retained minimum-norm8 translation classes on the direct11952 alternate-Q80 parent; compare complete projective branch images over the fixed original parameter.',
            'acceptance':['same actual rational quadratic squareclass','two exact independent new directions modulo inherited MW17','smooth rational base or positive-rank genus-one base'],
            'quarantine':'Only the generic equation, generic section coordinates, generic height Gram and complete generic norm8 words are exported. Previous product targets and exceptional specializations are not read.'}
    write_new(path/'input.json.gz',packet)
    print(json.dumps({'status':'FROZEN','pencils':len(words),'pairs':len(words)*(len(words)-1)//2,'primes':packet['primes']}),flush=True)


def degree_height(point):
    if point.is_zero():return 0
    x=point[0]
    return int(max(4+x.denominator().degree(),x.numerator().degree()))


def context(packet,p):
    R=PolynomialRing(GF(p),'t');K=R.fraction_field();t=R.gen()
    cv=lambda xs:R([GF(p)(QQ(x)) for x in xs])
    A,B=cv(packet['A']),cv(packet['B']);D=4*A**3+27*B**2
    assert D.degree()==24 and not D.gcd(D.derivative()).degree(),'bad model reduction'
    E=EllipticCurve(K,[A,B])
    def rf(row):
        numerator=list(map(QQ,row['numerator_coefficients_low_to_high']))
        denominator=list(map(QQ,row['denominator_coefficients_low_to_high']))
        scale=lcm([x.denominator() for x in numerator+denominator])
        values=[int(scale*x) for x in numerator+denominator]
        content=gcd(values);values=[x//content for x in values]
        n=R(values[:len(numerator)]);d=R(values[len(numerator):])
        if not d:raise ArithmeticError('generic basis reduces to the zero section; defer this prime')
        return K(n)/d
    basis=[E(rf(r['X']),rf(r['Y'])) for r in packet['basis']]
    heights=[degree_height(P) for P in basis]
    for i,P in enumerate(basis):
        assert heights[i]==packet['gram'][i][i],'basis height changed'
        for j,Q in enumerate(basis[:i]):
            assert (degree_height(P+Q)-heights[i]-heights[j])==2*packet['gram'][i][j],'basis Gram changed'
    maximum=max(abs(c) for w in packet['words'] for c in w)
    multiples=[{n:n*P for n in range(-maximum,maximum+1)} for P in basis]
    @lru_cache(maxsize=16384)
    def partial(terms):
        if not terms:return E(0)
        i,n=terms[-1]
        return partial(terms[:-1])+multiples[i][n]
    def trace(word):
        terms=tuple((i,n) for i,n in enumerate(word) if n);cut=len(terms)//2
        return partial(terms[:cut])+partial(terms[cut:])
    return R,K,A,B,E,basis,trace


def compile_frame(ctx,word,helper,old):
    R,K,A,B,E,basis,tracer=ctx;t=R.gen();P=tracer(word)
    assert degree_height(P)==8,'trace reduction lost its height'
    X,Y=K(P[0]),K(P[1])
    for shift in (None,0,1,2):
        if shift is None:Xc,Yc,Ac=X,Y,A
        else:
            Xc=helper['invert_rational'](X(t+shift),4,R,K)
            Yc=helper['invert_rational'](Y(t+shift),6,R,K)
            Ac=old['reciprocal_polynomial'](A(t+shift),8,R)
        frame=helper['trace_chord_frame'](Xc,Yc,R)
        h,Nx,Ny,M=(frame[k] for k in ['h','Nx','Ny','M0'])
        if h.degree()==2:break
    else:raise ArithmeticError('no finite-pole chart among three shifts')
    qs=old['q_lambda_family'](h,Nx,Ny,M,Ac,R)
    if shift is not None:
        qs=[R(sum(q[i]*(t-shift)**(4-i) for i in range(5))) for q in qs]
    branch=[[int(q[i]) for q in qs] for i in range(5)]
    return {'chart_shift':shift,'h':[int(x) for x in h.list()],
            'Nx':[int(x) for x in Nx.list()],'Ny':[int(x) for x in Ny.list()],
            'M0':[int(x) for x in M.list()],'branch_matrix':branch}


def packed_images(branch,indices,p):
    matrices=np.asarray(branch,dtype=np.int64)
    values=np.arange(p,dtype=np.int64)
    ver=np.zeros((5,p+1),dtype=np.int64);ver[0,:p]=1
    for j in range(1,5):ver[j,:p]=(ver[j-1,:p]*values)%p
    ver[4,p]=1
    raw=np.matmul(matrices,ver)%p
    pivots=np.argmax(raw!=0,axis=1)
    chosen=np.take_along_axis(raw,pivots[:,None,:],axis=1)[:,0,:]
    if not np.all(chosen):raise ArithmeticError('a modular branch image is zero; stage cannot exclude pairs')
    inverses=np.asarray([0]+[pow(j,-1,p) for j in range(1,p)],dtype=np.int64)
    normalized=(raw*inverses[chosen][:,None,:])%p
    key=pivots.astype(np.uint64)*np.uint64(p**4)
    powers=[p**j for j in range(4)]
    # Drop the pivot coordinate (which equals1) and retain the other four.
    for pivot in range(5):
        mask=pivots==pivot
        sub=np.zeros_like(key)
        for j,coordinate in enumerate(i for i in range(5) if i!=pivot):
            sub+=normalized[:,coordinate,:].astype(np.uint64)*np.uint64(powers[j])
        key+=sub*mask
    packed=key*np.uint64(65536)+np.asarray(indices,dtype=np.uint64)[:,None]
    return packed.reshape(-1)


def collision_buckets(entries):
    entries.sort();keys=entries//np.uint64(65536)
    eq=keys[1:]==keys[:-1]
    starts=np.flatnonzero(eq & np.r_[True,~eq[:-1]])
    buckets=[]
    for start0 in starts:
        start=int(start0);stop=start+2
        while stop<len(entries) and keys[stop]==keys[start]:stop+=1
        members=sorted(set(int(x%np.uint64(65536)) for x in entries[start:stop]))
        if len(members)>1:buckets.append({'key':int(keys[start]),'members':members})
    return buckets


def merge_constraints(path,packet):
    n=len(packet['words']);all_bits=(1<<n)-1;stages=[]
    for p in packet['primes']:
        result=path/('prime-%d.json.gz'%p)
        if not result.exists():continue
        row=read(result);good=sum(1<<i for i in row['compiled_indices']);adj={}
        for bucket in row['collision_buckets']:
            bits=sum(1<<i for i in bucket['members'])
            for i in bucket['members']:adj[i]=adj.get(i,0)|bits
        stages.append((p,good,adj))
    survivors=[];vertices=set()
    for i in range(n):
        allowed=all_bits^((1<<(i+1))-1)
        for p,good,adj in stages:
            if (good>>i)&1:allowed&=(all_bits^good)|adj.get(i,0)
        while allowed:
            bit=allowed&-allowed;j=bit.bit_length()-1;allowed-=bit
            survivors.append([i,j]);vertices.add(i);vertices.add(j)
    return {'completed_primes':[p for p,_,_ in stages],'total_pairs':n*(n-1)//2,
            'surviving_pairs':survivors,'next_vertices':sorted(vertices),
            'status':'FINITE_CANDIDATES_REQUIRE_EXACT_RESOLUTION' if survivors else 'ALL_PAIRS_EXCLUDED_BY_PROJECTIVE_REDUCTION'}


def run_prime(path,p):
    packet=read(path/'input.json.gz');assert p in packet['primes']
    limit=packet['limits']['cpu_seconds_per_prime']
    resource.setrlimit(resource.RLIMIT_CPU,(limit,limit+5))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    start=time.process_time();ctx=context(packet,p);helper=runpy.run_path(str(HELPER));old=runpy.run_path(str(INVERSION))
    earlier=[q for q in packet['primes'] if q<p and (path/('prime-%d.json.gz'%q)).exists()]
    indices=merge_constraints(path,packet)['next_vertices'] if earlier else list(range(len(packet['words'])))
    directory=path/('prime-%d'%p);directory.mkdir(parents=True,exist_ok=True)
    packed_parts=[];all_indices=[];chunk_size=packet['limits']['checkpoint_size'];deferrals=[]
    for offset in range(0,len(indices),chunk_size):
        chosen=indices[offset:offset+chunk_size];checkpoint=directory/('%06d.json.gz'%offset)
        if checkpoint.exists():records=read(checkpoint)['records']
        else:
            records=[]
            for i in chosen:
                record=compile_frame(ctx,packet['words'][i],helper,old);record['index']=i;records.append(record)
            write_new(checkpoint,{'input_sha256':sha256((path/'input.json.gz').read_bytes()).hexdigest(),
                                  'prime':p,'indices':chosen,'records':records})
        assert [r['index'] for r in records]==chosen,'checkpoint coverage changed'
        packed_parts.append(packed_images([r['branch_matrix'] for r in records],chosen,p));all_indices.extend(chosen)
        if offset%4096==0 or offset+len(chosen)==len(indices):
            print(json.dumps({'prime':p,'compiled':len(all_indices),'selected':len(indices),'cpu_seconds':round(time.process_time()-start,3)}),flush=True)
    entries=np.concatenate(packed_parts) if packed_parts else np.asarray([],dtype=np.uint64)
    packed_parts.clear();buckets=collision_buckets(entries)
    result={'schema':'q80-complete-genus-one-prime-v1','input_sha256':sha256((path/'input.json.gz').read_bytes()).hexdigest(),
            'prime':p,'compiled_indices':all_indices,'collision_buckets':buckets,
            'projective_point_count':len(entries),'deferrals':deferrals,
            'cpu_seconds':time.process_time()-start,
            'zero_branch_vectors':0,'checkpoint_files':sorted(x.name for x in directory.glob('*.json.gz'))}
    write_new(path/('prime-%d.json.gz'%p),result)
    merged=merge_constraints(path,packet);merged['input_sha256']=result['input_sha256']
    write_new(path/('after-%d.json'%p),merged)
    print(json.dumps({'prime':p,'collision_buckets':len(buckets),'surviving_pairs':len(merged['surviving_pairs']),
                      'next_vertices':len(merged['next_vertices']),'cpu_seconds':result['cpu_seconds']}),flush=True)


def finish_embeddings(path):
    packet=read(path/'input.json.gz');audit=read(path/'matrix-rank-audit.json')
    indices=audit['remaining_indices'];p=541
    assert p in packet['primes']
    frozen={'schema':'q80-remaining-embedding-input-v1','indices':indices,'prime':p,
            'source_input_sha256':sha256((path/'input.json.gz').read_bytes()).hexdigest(),
            'audit_sha256':sha256((path/'matrix-rank-audit.json').read_bytes()).hexdigest(),
            'scope':'Certify within-pencil injectivity for the five classes whose recorded reductions have not yet supplied an invertible branch matrix.',
            'limits':{'cpu_seconds':20,'address_space_gib':4}}
    write_new(path/'embedding-input.json',frozen)
    resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    start=time.process_time();ctx=context(packet,p);helper=runpy.run_path(str(HELPER));old=runpy.run_path(str(INVERSION))
    records=[]
    for i in indices:
        r=compile_frame(ctx,packet['words'][i],helper,old);r['index']=i
        assert matrix(GF(p),r['branch_matrix']).det(),'remaining embedding not certified'
        records.append(r)
    write_new(path/'embedding-result.json.gz',{'schema':'q80-remaining-embedding-result-v1',
              'input_sha256':sha256((path/'embedding-input.json').read_bytes()).hexdigest(),
              'prime':p,'records':records,'cpu_seconds':time.process_time()-start})
    print(json.dumps({'embedding_count':len(records),'prime':p,'status':'ALL_FIVE_INVERTIBLE'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['freeze','prime','embeddings'])
    parser.add_argument('--output',type=Path,default=DEFAULT);parser.add_argument('--prime',type=int)
    args=parser.parse_args()
    if args.mode=='freeze':freeze(args.output)
    elif args.mode=='embeddings':finish_embeddings(args.output)
    else:run_prime(args.output,args.prime)
