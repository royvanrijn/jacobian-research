#!/usr/bin/env python3
"""Factorization-free good-prime parity for generic linear relation roots."""
import argparse
from fractions import Fraction as Q
from math import gcd, isqrt, lcm
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('GENERIC_CHORD_RAMIFICATION_PROTOCOL.json')
PANEL=r.OUT/'rank_jump_fresh_governing_panel_inputs_v1.json'
BOUNDARY=r.OUT/'rank_jump_fresh_strict_boundary_coordinates_v1.json'
SUPPLEMENT=r.OUT/'rank_jump_fresh_retained_factor_supplement_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-generic-chord-ramification-v1'
OUT=r.OUT/'rank_jump_generic_chord_ramification_v1.json'


def model(token):
    row=next(c for c in r.read(PANEL)['cases'] if c['token']==token)
    coeff,pts=r.short(row['model'],row['generic_sections'])
    scale=lcm(*(Q(x).denominator for x in coeff))
    A=int(Q(coeff[3])*scale**4); B=int(Q(coeff[4])*scale**6)
    pts=[(Q(x)*scale**2,Q(y)*scale**3) for x,y in pts]
    old=next(c for c in r.read(BOUNDARY)['rows'] if c['token']==token)
    if old['status']!='PASS':
        old=next(c['boundary'] for c in r.read(SUPPLEMENT)['rows'] if c['token']==token)
    assert old['status']=='PASS'
    S=old['S_finite']
    disc=abs(-4*A**3-27*B**2)
    for p in S:
        while disc%p==0: disc//=p
    assert disc==1 and 2 in S
    return A,B,pts,S


def forms(A,B,pts):
    found={}
    def add(i,j,sign):
        x,y=pts[i]; u,v=pts[j]; v*=sign
        if i==j:
            slope=(3*x*x+A)/(2*y)
        else:
            assert x!=u
            slope=(v-y)/(u-x)
        intercept=y-slope*x
        if not slope: return
        a,b=intercept,slope
        d=lcm(a.denominator,b.denominator)
        a,b=int(a*d),int(b*d); g=gcd(a,b)
        a//=g; b//=g
        if b<0:a,b=-a,-b
        # Exact chord intersection identity; the third intersection is a
        # generic group sum, not a searched or claimed exceptional point.
        w=slope*slope-x-u
        assert (slope*w+intercept)**2==w**3+A*w+B
        found.setdefault((a,b),[]).append([i,j,sign])
    for i in range(len(pts)):
        add(i,i,1)
        for j in range(i+1,len(pts)):
            for sign in (1,-1):add(i,j,sign)
    return [{'a':str(a),'b':str(b),'sources':src,
             'norm':str(a**3+A*a*b*b-B*b**3)} for (a,b),src in sorted(found.items())]


def insert(basis,n):
    n=abs(n)
    if n<=1:return
    i=0
    while i<len(basis):
        a=basis[i];d=gcd(a,n)
        if d==1:i+=1;continue
        if d==a:
            while n%a==0:n//=a
            if n==1:return
            i+=1;continue
        basis.pop(i)
        insert(basis,d);insert(basis,a//d)
        insert(basis,n)
        return
    basis.append(n)
    assert len(basis)<=4096


def exponents(n,atoms):
    n=abs(n);out={}
    for j,a in enumerate(atoms):
        e=0
        while n%a==0:n//=a;e+=1
        if e:out[j]=e
    assert n==1
    return out


def nullspace(rows):
    piv={};ker=[]
    for i,v in enumerate(rows):
        word=1<<i
        while v:
            k=v.bit_length()-1
            if k not in piv:piv[k]=(v,word);break
            w,c=piv[k];v^=w;word^=c
        if not v:ker.append(word)
    return ker


def worker(token):
    A,B,pts,S=model(token);items=forms(A,B,pts);basis=[]
    for row in items:
        n=abs(int(row['norm']));assert n
        for p in S:
            while n%p==0:n//=p
        row['outside_S_norm']=str(n)
        insert(basis,n)
    refinements=0
    while True:
        basis.sort();es=[exponents(int(x['outside_S_norm']),basis) for x in items]
        split=None
        for j,c in enumerate(basis):
            active=[i for i,e in enumerate(es) if e.get(j,0)%2]
            for pos,i in enumerate(active):
                a,b=int(items[i]['a']),int(items[i]['b'])
                for k in active[pos+1:]:
                    d=gcd(c,a*int(items[k]['b'])-b*int(items[k]['a']))
                    if d not in (1,c):split=(j,d);break
                if split:break
            if split:break
        if split is None:break
        j,d=split;c=basis.pop(j);insert(basis,d);insert(basis,c//d);refinements+=1
    signatures=[0]*len(items);blocks=[]
    for j,c in enumerate(basis):
        if isqrt(c)**2==c:continue
        active=[i for i,e in enumerate(es) if e.get(j,0)%2]
        groups={}
        for i in active:
            a,b=int(items[i]['a']),int(items[i]['b']);assert gcd(b,c)==1
            root=(-a*pow(b,-1,c))%c
            assert (root**3+A*root+B)%c==0
            groups.setdefault(root,[]).append(i)
        roots=sorted(groups)
        assert len(roots)<=3
        assert all(gcd(x-y,c)==1 for i,x in enumerate(roots) for y in roots[i+1:])
        offset=3*len(blocks)
        for k,root in enumerate(roots):
            for i in groups[root]:signatures[i]^=(7^(1<<k))<<offset
        blocks.append({'atom_index':j,'roots':[str(x) for x in roots],
                       'groups':[groups[x] for x in roots]})
    kernel=nullspace(signatures)
    return {'status':'PASS','token':token,'A':str(A),'B':str(B),'S_finite':S,
        'generic_dimension':len(pts),'forms':items,'atoms':list(map(str,basis)),
        'valuations':[{str(j):e for j,e in row.items()} for row in es],
        'blocks':blocks,'resultant_refinements':refinements,'signatures':list(map(str,signatures)),
        'outside_S_parity_rank':r.rank(signatures),'coefficient_kernel_dimension':len(kernel),
        'kernel_masks':list(map(str,kernel)),
        'Selmer_dimension_beyond_G':'UNKNOWN' if kernel else 0}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in r.read(PROTOCOL)['cases']:
        path=WORK/f'{token}.json'
        if not path.exists():
            with (WORK/f'{token}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,__file__,'worker','--token',token],stdout=log,stderr=log,timeout=60)
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='60-second timeout'
            if reason:r.write_new(path,{'status':'UNKNOWN','token':token,'reason':reason})
        row=r.read(path);rows.append(row)
        print(token,row['status'],row.get('outside_S_parity_rank'),row.get('coefficient_kernel_dimension'),flush=True)
    r.write_new(OUT,{'schema':'rank-jump.generic-chord-ramification.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,PANEL,BOUNDARY,SUPPLEMENT]}})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='capture':capture()
    else:r.write_new(WORK/f'{args.token}.json',worker(args.token))
