#!/usr/bin/env python3
"""Independent finite-arithmetic replay of the complete Q80 pencil comparison.

Trace identification uses a height bound and at least15 smooth fibres.
Polynomial identities use exact interpolation with explicit degree bounds.
The written good-reduction and height arguments remain part of the proof.
"""
import argparse
from collections import Counter
import csv
from fractions import Fraction as Q
from functools import lru_cache
import gzip
from hashlib import sha256
import importlib.util
from math import comb
from pathlib import Path
import json
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1'
DIRECT=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
TABLE=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv'
SPEC=importlib.util.spec_from_file_location('q80_finite_poly',ROOT/'elkies-k3/scripts/verify_r17_mestre_shared_twist.py')
F=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(F)
require=F.require


def read(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)


def ev(f,t,p):
    out=0
    for x in reversed(f):out=(out*t+x)%p
    return out


def divexact(a,b,p):
    require(b,'zero polynomial divisor');a=list(a);q=[0]*max(0,len(a)-len(b)+1);inv=pow(b[-1],-1,p)
    while a and len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]*inv%p;q[k]=c
        for j,x in enumerate(b):a[k+j]=(a[k+j]-c*x)%p
        a=F.trim(a)
    require(not a,'inexact polynomial division');return F.trim(q)


def primitive_rf(row,p):
    nums=list(map(Q,row['numerator_coefficients_low_to_high']));dens=list(map(Q,row['denominator_coefficients_low_to_high']))
    def valuation(x):
        if not x:return 10**9
        n,d=x.numerator,x.denominator;v=0
        while n%p==0:n//=p;v+=1
        while d%p==0:d//=p;v-=1
        return v
    v=min(map(valuation,nums+dens));scale=Q(p)**(-v)
    cv=lambda xs:F.trim([int((x*scale).numerator*pow((x*scale).denominator,-1,p)%p) for x in xs])
    n,d=cv(nums),cv(dens);require(d,'zero generic-coordinate denominator reduction')
    return normalize_rf((n,d),p)


def normalize_rf(f,p):
    n,d=f;require(d,'zero rational-function denominator')
    if not n:return ([],[1])
    g=F.ff_gcd(n,d,p);n,d=divexact(n,g,p),divexact(d,g,p)
    c=pow(d[-1],-1,p);return ([c*x%p for x in n],[c*x%p for x in d])


def context(packet,p):
    require(p>25 and all(p%d for d in range(2,int(p**0.5)+1)),'invalid model prime')
    cv=lambda xs:F.trim([Q(x).numerator*pow(Q(x).denominator,-1,p)%p for x in xs])
    mul=lambda a,b:F.ff_mul(a,b,p);sub=lambda a,b:F.ff_sub(a,b,p)
    scale=lambda a,n:F.trim([n*x%p for x in a]);add=lambda a,b:sub(a,scale(b,-1))
    A,B=cv(packet['A']),cv(packet['B']);D=add(scale(mul(mul(A,A),A),4),scale(mul(B,B),27))
    require(len(A)<=9 and len(B)<=13,'K3 coefficient degree bounds')
    require(len(D)==25 and len(F.ff_gcd(D,F.trim([i*x%p for i,x in enumerate(D)][1:]),p))==1,'bad model reduction')
    rf=lambda n,d=([1]):normalize_rf((n,d),p)
    radd=lambda x,y:rf(add(mul(x[0],y[1]),mul(y[0],x[1])),mul(x[1],y[1]))
    rneg=lambda x:(scale(x[0],-1),x[1])
    rsub=lambda x,y:radd(x,rneg(y))
    rmul=lambda x,y:rf(mul(x[0],y[0]),mul(x[1],y[1]))
    rdiv=lambda x,y:rf(mul(x[0],y[1]),mul(x[1],y[0]))
    rscale=lambda x,n:rf(scale(x[0],n),x[1])
    basis=[(primitive_rf(r['X'],p),primitive_rf(r['Y'],p)) for r in packet['basis']]
    for x,y in basis:
        require(rmul(y,y)==radd(radd(rmul(rmul(x,x),x),rmul((A,[1]),x)),(B,[1])),'generic section equation over finite function field')
    height=lambda x:0 if x is None else max(4+len(x[1])-1,len(x[0])-1)
    def sum_x(P,Q):
        x,y=P;u,v=Q
        if x==u:
            if y==rneg(v):return None
            m=rdiv(radd(rscale(rmul(x,x),3),(A,[1])),rscale(y,2))
        else:m=rdiv(rsub(v,y),rsub(u,x))
        return rsub(rsub(rmul(m,m),x),u)
    hs=[height(P[0]) for P in basis]
    for i,P in enumerate(basis):
        require(hs[i]==packet['gram'][i][i],'basis height changed under reduction')
        for j,Q0 in enumerate(basis[:i]):
            require(height(sum_x(P,Q0))-hs[i]-hs[j]==2*packet['gram'][i][j],'basis Gram changed under reduction')
    sites=[t for t in range(p) if ev(D,t,p)][:17];require(len(sites)==17,'not enough smooth test fibres')
    def add_points(P,Q0,a):
        if P is None:return Q0
        if Q0 is None:return P
        x,y=P;u,v=Q0
        if x==u:
            if (y+v)%p==0:return None
            m=(3*x*x+a)*pow(2*y,-1,p)%p
        else:m=(v-y)*pow((u-x)%p,-1,p)%p
        z=(m*m-x-u)%p;return z,(m*(x-z)-y)%p
    terms=[tuple((i,n) for i,n in enumerate(w) if n) for w in packet['words']]
    parts=[(x[:len(x)//2],x[len(x)//2:]) for x in terms]
    maximum=max(abs(c) for w in packet['words'] for c in w);sums=[]
    for t in sites:
        aa=ev(A,t,p);multiples=[]
        for x,y in basis:
            xd,yd=ev(x[1],t,p),ev(y[1],t,p)
            if xd:
                require(yd,'inconsistent affine basis reduction');P=(ev(x[0],t,p)*pow(xd,-1,p)%p,ev(y[0],t,p)*pow(yd,-1,p)%p)
            else:
                require(not yd,'inconsistent zero-section reduction');P=None
            tab={0:None};value=None
            for k in range(1,maximum+1):
                value=add_points(value,P,aa);tab[k]=value;tab[-k]=None if value is None else (value[0],-value[1]%p)
            multiples.append(tab)
        def make_sum(tables,a):
            @lru_cache(maxsize=4096)
            def partial(ts):
                if not ts:return None
                i,n=ts[-1];return add_points(partial(ts[:-1]),tables[i][n],a)
            return lambda i:add_points(partial(parts[i][0]),partial(parts[i][1]),a)
        sums.append(make_sum(multiples,aa))
    chart_coefficients={None:(A,B)}
    for shift in (0,1,2):
        def change(f,bound):
            out=[0]*(bound+1)
            for i,c in enumerate(f):
                for k in range(i+1):out[bound-i+k]=(out[bound-i+k]+c*comb(i,k)*pow(shift,k,p))%p
            return F.trim(out)
        chart_coefficients[shift]=(change(A,8),change(B,12))
    return {'p':p,'A':A,'B':B,'sites':sites,'word_sums':sums,'chart_coefficients':chart_coefficients}


def full_rank(matrices,p):
    a=np.array(matrices,dtype=np.int64,copy=True);n=len(a);good=np.ones(n,dtype=bool);ix=np.arange(n)
    inverses=np.asarray([0]+[pow(i,-1,p) for i in range(1,p)],dtype=np.int64)
    for k in range(5):
        non=a[:,k:,k]!=0;good &= np.any(non,axis=1);pivot=k+np.argmax(non,axis=1)
        tmp=a[ix,k,:].copy();a[ix,k,:]=a[ix,pivot,:];a[ix,pivot,:]=tmp
        factor=a[:,k+1:,k]*inverses[a[:,k,k]][:,None]%p
        a[:,k+1:,:]=(a[:,k+1:,:]-factor[:,:,None]*a[:,k,None,:])%p
    return good


def validate_records(records,ctx,identify=True):
    p=ctx['p'];count=len(records)
    for r in records:
        require(r['chart_shift'] in (None,0,1,2),'unimplemented chart')
        require(len(r['h'])==3 and r['h'][-1]!=0 and len(r['Nx'])<=9 and len(r['Ny'])<=13 and len(r['M0'])<=4,'frame degree bounds')
        require(len(r['branch_matrix'])==5 and all(len(x)==5 for x in r['branch_matrix']),'branch matrix shape')
        require(all(isinstance(c,int) and 0<=c<p for k in ['h','Nx','Ny','M0'] for c in r[k]),'noncanonical finite coefficient')
        require(all(isinstance(c,int) and 0<=c<p for row in r['branch_matrix'] for c in row),'noncanonical branch coefficient')
        require(len(F.ff_gcd(r['h'],r['Nx'],p))==1,'cancelled trace pole divisor')
        mod=F.ff_mul(r['h'],r['h'],p)
        congruence=F.ff_sub(F.ff_mul(r['M0'],r['Nx'],p),[-x%p for x in r['Ny']],p)
        require(not F.ff_rem(congruence,mod,p),'regular chord congruence')
    padded=lambda key,n:np.asarray([r[key]+[0]*(n-len(r[key])) for r in records],dtype=np.int64)
    basis=lambda n,sites:np.asarray([[pow(int(x),j,p) for x in sites] for j in range(n)],dtype=np.int64)
    ts=list(range(25));h=padded('h',3)@basis(3,ts)%p;nx=padded('Nx',9)@basis(9,ts)%p;ny=padded('Ny',13)@basis(13,ts)%p;m=padded('M0',4)@basis(4,ts)%p
    aa=np.asarray([[ev(ctx['chart_coefficients'][r['chart_shift']][0],t,p) for t in ts] for r in records],dtype=np.int64)
    bb=np.asarray([[ev(ctx['chart_coefficients'][r['chart_shift']][1],t,p) for t in ts] for r in records],dtype=np.int64)
    h2=h*h%p;h4=h2*h2%p;h6=h4*h2%p;h8=h4*h4%p
    require(np.all(ny*ny%p==(nx*nx%p*nx+aa*nx%p*h4+bb*h6)%p),'degree24 trace Weierstrass identity')
    local=np.asarray([r['branch_matrix'] for r in records],dtype=np.int64)
    for shift in (0,1,2):
        transform=np.zeros((5,5),dtype=np.int64)
        for i in range(5):
            for k in range(i+1):transform[4-i+k,i]=comb(i,k)*pow(shift,k,p)%p
        mask=np.asarray([r['chart_shift']==shift for r in records]);local[mask]=np.matmul(transform,local[mask])%p
    q=np.matmul(local.transpose(0,2,1),basis(5,ts))%p;m2=m*m%p
    right=np.stack([(m2*m2-6*m2*nx-8*m*ny-3*nx*nx-4*aa*h4)%p,
                    ((4*m2*m-12*m*nx-8*ny)%p)*h2%p,
                    6*(m2-nx)*h4%p,4*m*h6%p,h8],axis=1)
    require(np.all((q[:,:,:17]*h6[:,None,:17]-right[:,:,:17])%p==0),'degree16 branch coefficient identities')
    if identify:
        for r in records:
            matched=0;shift=r['chart_shift']
            for t,word_sum in zip(ctx['sites'],ctx['word_sums']):
                if shift is None:
                    he=ev(r['h'],t,p);xe=ev(r['Nx'],t,p);ye=ev(r['Ny'],t,p)
                else:
                    value=(t-shift)%p
                    he=sum(c*pow(value,2-j,p) for j,c in enumerate(r['h']))%p
                    xe=sum(c*pow(value,8-j,p) for j,c in enumerate(r['Nx']))%p
                    ye=sum(c*pow(value,12-j,p) for j,c in enumerate(r['Ny']))%p
                if not he:continue
                inv=pow(he,-1,p);P=(xe*inv*inv%p,ye*inv*inv*inv%p)
                require(P==word_sum(r['index']),'height-bounded trace identification failed')
                matched+=1
            require(matched>=15,'too few smooth trace-identification fibres')
    return full_rank([r['branch_matrix'] for r in records],p)


def images(records,p):
    # Independent Horner evaluation and last-nonzero normalization.
    mats=np.asarray([r['branch_matrix'] for r in records],dtype=np.int64)
    values=np.arange(p,dtype=np.int64);out=np.repeat(mats[:,:,4,None],p,axis=2)
    for j in (3,2,1,0):out=(out*values+mats[:,:,j,None])%p
    out=np.concatenate((out,mats[:,:,4,None]),axis=2)
    pivots=4-np.argmax(out[:,::-1,:]!=0,axis=1);pivot_values=np.take_along_axis(out,pivots[:,None,:],axis=1)[:,0,:]
    require(np.all(pivot_values),'zero projective branch image')
    inv=np.asarray([0]+[pow(i,-1,p) for i in range(1,p)],dtype=np.int64)
    out=out*inv[pivot_values][:,None,:]%p;keys=pivots.astype(np.uint64)*np.uint64(p**4)
    for pivot in range(5):
        mask=pivots==pivot;key=np.zeros_like(keys)
        for coordinate in range(5):
            if coordinate!=pivot:key=key*np.uint64(p)+out[:,coordinate,:].astype(np.uint64)
        keys+=key*mask
    return (keys*np.uint64(65536)+np.asarray([r['index'] for r in records],dtype=np.uint64)[:,None]).ravel()


def buckets(entries):
    entries.sort();keys=entries//np.uint64(65536);same=np.flatnonzero(keys[:-1]==keys[1:]);result=[];last_end=0
    for raw in same:
        start=int(raw)
        if start<last_end:continue
        end=start+2
        while end<len(keys) and keys[end]==keys[start]:end+=1
        members=tuple(sorted(set(int(x%np.uint64(65536)) for x in entries[start:end])));last_end=end
        if len(members)>1:result.append(members)
    return result


def remaining_pairs(n,stages):
    all_bits=(1<<n)-1;out=[]
    for i in range(n):
        possible=all_bits^((1<<(i+1))-1)
        for good,adj in stages:
            if good>>i&1:possible&=(all_bits^good)|adj.get(i,0)
        while possible:
            bit=possible&-possible;out.append([i,bit.bit_length()-1]);possible-=bit
    return out


def require_completed_stages(path):
    for name in ['prime-521.json.gz','prime-523.json.gz','prime-541.json.gz','embedding-input.json','embedding-result.json.gz']:
        require((path/name).is_file(),'missing completed prime or embedding certificate: '+name)


def verify(path,progress=False):
    require_completed_stages(path)
    packet=read(path/'input.json.gz');require(packet['schema']=='q80-complete-genus-one-pairs-input-v1','input schema')
    require(packet['primes']==[509,521,523,541],'frozen prime panel')
    for source in (DIRECT,TABLE):require(sha256(source.read_bytes()).hexdigest()==packet['source_hashes'][str(source.relative_to(ROOT))],'generic source hash binding')
    direct=read(DIRECT)
    require(packet['A']==direct['weierstrass_model']['A_coefficients_low_to_high'] and packet['B']==direct['weierstrass_model']['B_coefficients_low_to_high'],'model projection')
    require(packet['basis']==[{'X':r['X'],'Y':r['Y']} for r in direct['sections']['records']],'basis projection')
    require(packet['gram']==[[int(Q(x)) for x in row] for row in direct['sections']['height_gram']],'generic height binding')
    with TABLE.open() as f:table=list(csv.DictReader(f,delimiter='\t'))
    require(packet['words']==[[int(x) for x in r['section_basis_w'].split()] for r in table] and packet['orbit_masks']==[int(r['orbit_mask']) for r in table],'complete generic inventory binding')
    n=len(packet['words']);require(n==63917 and len({tuple(x%2 for x in w) for w in packet['words']})==n,'complete distinct parity coverage')
    W=np.asarray(packet['words'],dtype=np.int64);G=np.asarray(packet['gram'],dtype=np.int64)
    require(np.all(np.sum((W@G)*W,axis=1)==8),'generic word norms')
    frozen=sha256((path/'input.json.gz').read_bytes()).hexdigest();stages=[];embeddings=set();summaries=[];total_frames=0;contexts={}
    for p in (521,523,541):
        row=read(path/('prime-%d.json.gz'%p));require(row['input_sha256']==frozen and row['prime']==p,'prime input binding')
        expected=list(range(n)) if not stages else sorted(set(i for pair in remaining_pairs(n,stages) for i in pair))
        require(row['compiled_indices']==expected,'adaptive prime coverage')
        files=row['checkpoint_files'];require(len(files)==(len(expected)+511)//512,'missing checkpoint')
        ctx=context(packet,p);contexts[p]=ctx;parts=[];seen=[]
        for block,name in enumerate(files):
            checkpoint=read(path/('prime-%d'%p)/name)
            subset=expected[block*512:(block+1)*512];records=checkpoint['records']
            require(checkpoint['input_sha256']==frozen and checkpoint['prime']==p and checkpoint['indices']==subset,'checkpoint binding')
            require([r['index'] for r in records]==subset,'checkpoint record coverage')
            ranks=validate_records(records,ctx)
            embeddings.update(r['index'] for r,yes in zip(records,ranks) if yes)
            parts.append(images(records,p));seen.extend(subset);total_frames+=len(records)
            if progress and (block%16==0 or len(seen)==len(expected)):print(json.dumps({'prime':p,'independently_verified':len(seen),'selected':len(expected)}),flush=True)
        require(seen==expected,'missing or repeated pencil')
        entries=np.concatenate(parts) if parts else np.asarray([],dtype=np.uint64);parts.clear();calculated=buckets(entries)
        require(Counter(calculated)==Counter(tuple(x['members']) for x in row['collision_buckets']),'projective collision buckets differ')
        require(len(entries)==row['projective_point_count'],'projective point coverage')
        adj={}
        for members in calculated:
            bits=sum(1<<i for i in members)
            for i in members:adj[i]=adj.get(i,0)|bits
        stages.append((sum(1<<i for i in expected),adj));left=remaining_pairs(n,stages)
        summary={'prime':p,'pencils':len(expected),'projective_points':len(entries),'remaining_pairs':len(left)};summaries.append(summary)
        if progress:print(json.dumps(summary),flush=True)
        del entries
    require(not remaining_pairs(n,stages),'common-cover candidates remain')
    supplement=read(path/'embedding-input.json');extra=read(path/'embedding-result.json.gz')
    require(supplement['source_input_sha256']==frozen and extra['input_sha256']==sha256((path/'embedding-input.json').read_bytes()).hexdigest(),'embedding input binding')
    missing=sorted(set(range(n))-embeddings)
    require(supplement['indices']==missing and [r['index'] for r in extra['records']]==missing and supplement['prime']==extra['prime']==541,'remaining embedding coverage')
    require(all(validate_records(extra['records'],contexts[541])),'uncertified within-pencil injectivity')
    return {'status':'PASS_COMPLETE_Q80_GENUS_ONE_PENCIL_INJECTIVITY','pencils':n,'pair_count':n*(n-1)//2,
            'prime_stages':summaries,'verified_trace_frames':total_frames+len(missing),'additional_embedding_certificates':len(missing),
            'all_pencils_injective':True,'distinct_classes_have_disjoint_rational_branch_images':True,'goal_complete':False,
            'inherited_dependencies':['published direct11952 MW17 equation and height lattice','complete retained63917 minimum-norm8 parity inventory'],
            'scope':'No two rational members of distinct pencils, and no distinct members of one pencil, in this complete norm8/pole-order-zero layer share a projective branch quartic over the fixed parameter. Other carrier layers and arbitrary quadratic-twist sections are not excluded.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input-dir',type=Path,default=DEFAULT);parser.add_argument('--output',type=Path);parser.add_argument('--progress',action='store_true');args=parser.parse_args()
    start=time.monotonic();result=verify(args.input_dir,args.progress);result['elapsed_seconds']=round(time.monotonic()-start,3)
    if args.output:
        with args.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True))
