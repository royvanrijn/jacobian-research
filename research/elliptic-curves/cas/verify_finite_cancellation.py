#!/usr/bin/env python3
"""Portable replay of exact target and residue certificates, without point search.

Uses the older binary transform and group-law implementations, not the feature
producer's evaluator or local-tree recursion. Retained rank certificates remain
the independence authority; this verifier does not issue new rank claims.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
import gzip
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import time

from finite_cancellation_corpus import OUT, canonical, digest, write
from half_lattice_pointed_sieve import binary_transform
from verify_height_model_accessibility import add

def value(c,m,n):
    return c[0]*n**4+c[1]*m*n**3+c[2]*m*m*n*n+c[3]*m**3*n+c[4]*m**4
def valuation(n,p):
    if n==0:return 10**9
    k=0
    while n%p==0:k+=1;n//=p
    return k
def constant(c,p,digits=1):
    v=valuation(c[0],p)
    return v if v<10**9 and min((valuation(x,p) for x in c[1:]),default=10**9)>=v+digits else None
def verify_tree(tree,N,D,q):
    p=tree['prime'];tries={'affine':{},'infinity':{}};mass=F(0);leaves=tree['leaves']
    for leaf in leaves:
        k=leaf['depth'];r=leaf['residue'];mass+=F(leaf['mass'])
        assert F(leaf['mass'])==F(1,(p+1)*p**(k-1)) and 0<=r<p**k
        trie=tries[leaf['kind']]
        for j in range(k):
            assert 'terminal' not in trie
            trie=trie.setdefault((r//p**j)%p,{})
        assert not trie;trie['terminal']=True
        matrix=(p**k,r,0,1) if leaf['kind']=='affine' else (0,1,p**k,r)
        if leaf['kind']=='infinity':assert r%p==0
        nn=binary_transform(N,matrix);dd=binary_transform(D,matrix);qq=binary_transform(q,matrix)
        v=min(valuation(x,p) for x in nn+dd)
        assert leaf['vg']==v
        exact=constant(nn,p)==v or constant(dd,p)==v
        assert leaf['vg_exact']==exact
        if leaf['status']=='square':
            vq=constant(qq,p,3 if p==2 else 1);assert vq is not None and vq%2==0 and exact
            unit=qq[0]//p**vq
            assert (unit%8==1 if p==2 else pow(unit%p,(p-1)//2,p)==1)
            assert leaf['vq']==vq
        elif leaf['status']=='nonsquare':
            vq=constant(qq,p);assert vq is not None
            if vq%2==0:
                assert constant(qq,p,3 if p==2 else 1)==vq
                unit=qq[0]//p**vq
                assert (unit%8!=1 if p==2 else pow(unit%p,(p-1)//2,p)!=1)
        else:assert leaf['status']=='unknown' and leaf['censored']
    assert mass==1
    # Prefix disjointness plus total measure1 certifies the full projective
    # partition, including the reciprocal chart at infinity.
    return len(leaves)

def verify():
    start=time.process_time();corpus=json.loads(gzip.decompress((OUT/'corpus.json.gz').read_bytes()))
    index={r['id']:r for r in corpus};ph=digest((OUT/'protocol.json').read_bytes());ch=digest((OUT/'corpus.json.gz').read_bytes())
    files=list(sorted((OUT/'cases').glob('*.json.gz')));seen=set();count=Counter();inputs={}
    for path in files:
        raw=path.read_bytes();inputs[path.name]=digest(raw);packet=json.loads(gzip.decompress(raw))
        assert packet['protocol_sha256']==ph and packet['corpus_sha256']==ch
        if packet['status']!='PASS_EXACT_ACCESSIBILITY':count['censored_anchors']+=1;continue
        prep=packet['prepared'];case=index[packet['case_id']];ai=prep['anchor_index']
        assert (case['id'],ai) not in seen;seen.add((case['id'],ai))
        assert prep['curve']==case['curve'] and prep['anchor']==case['generic_points'][ai]
        assert digest(canonical(prep))==packet['prepared_sha256']
        A,B=map(F,prep['curve'][3:]);a,b=map(F,prep['anchor']);assert b*b==a**3+A*a+B
        coeffs=[]
        for model in prep['models']:
            N,D,q=([int(v) for v in model[k]] for k in ['N','D','q']);matrix=tuple(map(F,model['mapping']['matrix']))
            assert gcd(*(N+D))==1
            nn=binary_transform((a**3+4*B,4*a*b,6*a*a+4*A,4*b,a),matrix)
            dd=binary_transform((-3*a*a-4*A,-8*b,-6*a,0,1),matrix)
            den=lcm(*(v.denominator for v in nn+dd));ints=[int(v*den) for v in nn+dd];content=gcd(*ints)
            assert [v//content for v in ints]==N+D
            ratio=F(model['mapping']['square_ratio']);assert ratio>0
            assert list(dd)==[ratio*v for v in q]
            assert isqrt(ratio.numerator)**2==ratio.numerator and isqrt(ratio.denominator)**2==ratio.denominator
            for tree in model['local']:count['residue_leaves']+=verify_tree(tree,N,D,q)
            coeffs.append((N,D,q,matrix));count['models']+=1
        residuals={};cofactors=defaultdict(set);observed=set()
        for row in packet['observations']:
            ti,sgn,mi=(row[k] for k in ['target_index','sign','model_index']);assert sgn in [1,-1]
            assert (ti,sgn,mi) not in observed;observed.add((ti,sgn,mi))
            x,y=map(F,case['targets'][ti]);y*=sgn;assert y*y==x**3+A*x+B
            m,n=map(int,row['coordinate']);assert gcd(m,n)==1 and n>=0
            N,D,q,(u,v,w,z)=coeffs[mi]
            assert (u*m+v*n)*(x-a)==(w*m+z*n)*(y+b)
            if (ti,sgn) not in residuals:residuals[ti,sgn]=add(add((x,y),(x,y),A),(a,-b),A)
            R=residuals[ti,sgn];nn,dd,qq=value(N,m,n),value(D,m,n),value(q,m,n)
            assert R is not None and dd and F(nn,dd)==R[0] and qq>=0 and isqrt(qq)**2==qq
            g=gcd(nn,dd);H=max(abs(m),n);Hx=max(abs(R[0].numerator),R[0].denominator)
            assert str(H)==row['H'] and str(g)==row['g'] and str(Hx)==row['Hx']
            assert max(abs(nn),abs(dd))==g*Hx and F(row['S'])==F(g*Hx,H**4)
            rest=g
            for location,tree in zip(row['local'],prep['models'][mi]['local']):
                p=location['p'];v=valuation(g,p);assert v==location['vg'] and p==tree['prime'];rest//=p**v
                leaf=tree['leaves'][location['leaf']];mod=p**leaf['depth']
                if n%p:
                    assert leaf['kind']=='affine' and (m-leaf['residue']*n)%mod==0
                else:assert leaf['kind']=='infinity' and (n-leaf['residue']*m)%mod==0
                assert leaf['status']!='nonsquare'
            assert str(rest)==row['unprocessed_g_cofactor'];cofactors[ti,sgn].add(rest);count['target_model_evaluations']+=1
        assert len(observed)==2*len(case['targets'])*len(prep['models'])
        assert all(len(values)==1 for values in cofactors.values()),'unprocessed cancellation changed between p-neighbours'
        count['anchors']+=1;count['all_other_primes_equal_comparisons']+=sum(len(prep['models'])-1 for _ in cofactors)
    assert len(files)==2*len(corpus), 'incomplete declared anchor corpus'
    result={'status':'PASS_EXACT_FINITE_CANCELLATION_REPLAY','counts':dict(count),'cpu_seconds':time.process_time()-start,
        'protocol_sha256':ph,'corpus_sha256':ch,'case_files_sha256':inputs,'verifier_sha256':digest(Path(__file__).read_bytes()),
        'boundary':'Exact maps, signed points, gcd/S identities, all residue partitions and square exclusions replayed. All unprocessed-prime cancellation is equal between the compared models. Retained independent-point certificates supply rank provenance; no new rank replay or point search.'}
    write(OUT/'replay.json',result);print(json.dumps({k:v for k,v in result.items() if k!='case_files_sha256'}))

if __name__=='__main__':verify()
