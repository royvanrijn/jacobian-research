#!/usr/bin/env sage -python
"""Same pointed cancellation, different independence label after adjoining P.

This proves a subgroup-information boundary, not uselessness of local features
or absence of any unsearched rational point. No point search or factorization.
"""
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from math import gcd,lcm
from pathlib import Path
import time

from cancellation_cloud_training import OUT,RAW,CAS
from finite_cancellation_corpus import ROOT,canonical,digest,write


def main():
    from half_lattice_pointed_sieve import linear_combination
    from finite_cancellation_features import xmap,ev
    from cancellation_scheduler_prepare import new_write,V2
    from sage.all import EllipticCurve,QQ
    start=time.process_time();training=json.loads((OUT/'cloud-fit.json').read_text())
    row=next(r for r in training['records'] if r['directions']>=2)
    case=next(c for c in json.loads((V2/'inputs.json').read_text()) if c['id']==row['case'])
    seed=case['seed'];source=ROOT/row['source'];call=json.loads(source.read_text())
    firstpath=source.parent/'rank-input.json';first=json.loads(firstpath.read_text())
    wholepath=ROOT/row['packet'];whole=json.loads(wholepath.read_text())
    if digest(source.read_bytes())!=row['sha256'] or digest(wholepath.read_bytes())!=row['packet_sha256']:raise ArithmeticError('training witness changed')
    n=len(seed['points']);P=call['new_point']
    if first['points']!=seed['points']+[P] or first['curve']!=seed['curve'] or whole['points'][:n+1]!=first['points']:raise ArithmeticError('subgroup chain differs')
    checker=SourceFileLoader('subgroup_boundary_rank',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    ranks=[checker(packet) for packet in (seed,first,whole)]
    if [r['rank'] for r in ranks]!=[n,n+1,n+row['directions']]:raise ArithmeticError('independent rank chain failed')
    curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']];word=call['centre']
    Q=linear_combination(curve,basis,word)
    if Q!=linear_combination(curve,basis+[tuple(map(F,P))],word+[0]):raise ArithmeticError('anchor changed')
    E=EllipticCurve(QQ,list(curve));native=sum((int(c)*E(list(p)) for c,p in zip(word,basis)),E(0))
    if [str(native[0]),str(native[1])]!=list(map(str,Q)):raise ArithmeticError('independent anchor arithmetic differs')
    N,D,q=xmap(curve[3],curve[4],Q,call['mapping']);A,B,C,T=map(F,call['mapping']['matrix'])
    a,b=Q;x,y=map(F,P);coordinates=[]
    for sign in (1,-1):
        u=T*(sign*y+b)-B*(x-a);v=-C*(sign*y+b)+A*(x-a)
        den=lcm(u.denominator,v.denominator);m,k=int(u*den),int(v*den);g=gcd(m,k);m//=g;k//=g
        if m<0 or (m==0 and k<0):m,k=-m,-k
        coordinates.append((max(abs(m),abs(k)),m,k,sign))
    height,m,k,sign=min(coordinates);g=gcd(abs(int(ev(N,m,k))),abs(int(ev(D,m,k))))
    if height>call['search']['height_bound'] or g<=0:raise ArithmeticError('coordinate/cancellation witness invalid')
    new_write(OUT/'subgroup-boundary.json',{'status':'PASS_SAME_CANCELLATION_DIFFERENT_SUBGROUP_LABEL',
        'case':case['id'],'curve':seed['curve'],'anchor':list(map(str,Q)),'point':P,
        'old_word':word,'extended_word':word+[0],'rank_chain':ranks,
        'features':{'N':list(map(str,N)),'D':list(map(str,D)),'q':list(map(str,q)),
            'coordinate':[m,k],'sign':sign,'height':height,'primitive_gcd':str(g),'gcd_bits':g.bit_length()},
        'before':'P is independent of the certified starting subgroup, by the rank(n+1) certificate.',
        'after':'P belongs to the enlarged subgroup by its displayed generator list; the other certified cloud direction remains independent.',
        'source_sha256':{str(p.relative_to(ROOT)):digest(p.read_bytes()) for p in (source,firstpath,wholepath,Path(__file__))},
        'cpu_seconds':time.process_time()-start,
        'boundary':'For fixed E,Q,P, pointed forms, coordinates and all local cancellation valuations are unchanged when the known subgroup grows to contain P. Thus those features alone cannot determine whether P is a new direction for every subgroup. This is not a universal negative result about average predictive value or subgroup-aware policies.'})
    print(json.dumps({'status':'PASS_SAME_CANCELLATION_DIFFERENT_SUBGROUP_LABEL','ranks':[r['rank'] for r in ranks],'height':height,'gcd_bits':g.bit_length()}),flush=True)


if __name__=='__main__':main()
