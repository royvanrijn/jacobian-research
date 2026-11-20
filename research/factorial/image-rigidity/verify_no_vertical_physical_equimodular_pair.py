#!/usr/bin/env python3
"""Exact support for the written no-vertical physical pair proof."""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO.')
import argparse
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
PN='NO_VERTICAL_PHYSICAL_EQUIMODULAR_PAIR.md'
CN='verify_no_vertical_physical_equimodular_pair.py'
IN='no_vertical_pair_inputs.json'
LN='no_vertical_pair_source_ledger.json'
RN='no_vertical_pair_verification.json'
ROLES={'proof':('local',PN),'checker':('local',CN),'inputs':('local',IN),
 'stationary_context':('repo','research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md')}
SCOPE={'schema':1,'a_range':'real a>=1','d_range':'integer d>=2',
 'domain':'-pi/d<arg(S)<0 and its positive boundary, at defined points',
 'claim':'distinct physical points with equal Omega and equal |lambda| cannot have the same real part',
 'a_one_separate_proof':True,'global_fixed_sign_claimed':False,
 'global_branch_labeling_claimed':False,'activation_claimed':False,
 'global_SFC_claimed':False,'cutoff_claimed':False,
 'formal_verification':False,'external_review':False,'independent_replay':False}
def need(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):return sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def input_ok(d):
    need(type(d) is dict and d==SCOPE and all(type(d[k]) is type(v) for k,v in SCOPE.items()),'Exact scope/types')
def source_ok(d,repo):
    need(type(d) is dict and type(d.get('schema')) is int and d['schema']==1,'Exact schema1')
    rows=d.get('sources');need(type(rows) is list and len(rows)==4,'Four exact mandatory roles')
    seen=set()
    for row in rows:
        need(type(row) is dict and set(row)=={'role','root','path','sha256'},'Source fields')
        role=row['role'];need(type(role) is str and role in ROLES and role not in seen,'Unique role')
        need((row['root'],row['path'])==ROLES[role],'Exact role/path')
        need(type(row['sha256']) is str and re.fullmatch('[0-9a-f]{64}',row['sha256']),'SHA256 syntax')
        root=HERE if row['root']=='local' else repo;p=(root/row['path']).resolve()
        need(p.is_relative_to(root.resolve()) and p.is_file() and sha(p)==row['sha256'],'Missing/stale:'+role)
        seen.add(role)
    need(d.get('logical_inputs')==[],'Proof is elementary/self-contained; context is not a premise')
    need(d.get('stationary_context_only') is True,'Stationary context scope')
    for key in ('formal_verification','external_review','independent_replay'):need(d.get(key) is False,'Assurance flag')

NV=8
def C(c):return {} if not c else {(0,)*NV:Q(c)}
def X(i):
    v=[0]*NV;v[i]=1;return {tuple(v):Q(1)}
def clean(p):return {k:Q(v) for k,v in p.items() if v}
def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,Q(0))+v
    return clean(out)
def sc(p,c):return clean({k:v*c for k,v in p.items()})
def mul(*ps):
    out=C(1)
    for p in ps:
        nxt={}
        for a,c in out.items():
            for b,v in p.items():
                k=tuple(x+y for x,y in zip(a,b));nxt[k]=nxt.get(k,Q(0))+c*v
        out=clean(nxt)
    return out
def sq(p):return mul(p,p)
def diff(p,i):
    out={}
    for k,v in p.items():
        if k[i]:
            n=list(k);n[i]-=1;out[tuple(n)]=v*k[i]
    return out
def sub(p,i,q):
    out={}
    for k,v in p.items():
        r=list(k);power=r[i];r[i]=0;t={tuple(r):v}
        for _ in range(power):t=mul(t,q)
        out=add(out,t)
    return out
def algebra():
    x,y,a,D,R,kd,ka,z=[X(i) for i in range(NV)]
    A0=add(sq(x),sq(y));B0=add(sq(add(D,sc(x,-1))),sq(y));C0=add(sq(add(x,sc(a,-1))),sq(y))
    F1num=add(mul(D,C0),sc(A0,-1));F2num=add(A0,sc(mul(a,B0),-1))
    need(not add(F1num,sc(add(mul(add(D,C(-1)),sq(y)),mul(D,sq(add(x,sc(a,-1)))),sc(sq(x),-1)),-1)),'F1 exact numerator')
    need(not add(F2num,sc(add(mul(add(C(1),sc(a,-1)),sq(y)),sq(x),sc(mul(a,sq(add(D,sc(x,-1)))),-1)),-1)),'F2 exact numerator')
    theta=add(sc(mul(add(D,sc(a,-1)),x,B0,C0),-1),mul(add(D,sc(x,-1)),A0,C0),mul(add(x,sc(a,-1)),A0,B0))
    combination=add(mul(add(D,sc(x,-1)),F2num,C0),sc(mul(add(x,sc(a,-1)),F1num,B0),-1))
    need(theta==combination,'Full cleared theta derivative')
    eqC=add(sq(a),sc(mul(a,x),-2),mul(R,kd));eqE=add(sq(D),sc(mul(D,x),-2),mul(R,ka))
    rtarget=add(mul(R,add(mul(D,kd),sc(mul(a,ka),-1))),sc(mul(a,D,add(D,sc(a,-1))),-1))
    need(not add(rtarget,sc(mul(D,eqC),-1),mul(a,eqE)),'R formula from both modulus equations')
    xtarget=add(sc(mul(x,add(mul(D,kd),sc(mul(a,ka),-1))),2),sc(mul(sq(D),kd),-1),mul(sq(a),ka))
    need(not add(xtarget,mul(kd,eqE),sc(mul(ka,eqC),-1)),'X formula from both modulus equations')
    numerator=add(sc(mul(sq(D),add(mul(D,z),sc(a,-1))),2),sc(mul(D,add(mul(sq(D),z),sc(sq(a),-1))),-2))
    need(not add(numerator,sc(mul(a,D,add(D,sc(a,-1))),2)),'Exact negative derivative X(z)')
    zn=mul(a,add(D,C(-1)));zd=mul(D,add(a,C(-1)))
    lower=add(mul(add(D,a,C(-1)),add(mul(D,zn),sc(mul(a,zd),-1))),sc(mul(sq(D),zn),-1),mul(sq(a),zd))
    need(not lower,'Lower X endpoint')
    need(not add(mul(add(D,a),add(D,sc(a,-1))),sc(sq(D),-1),sq(a)),'Upper X endpoint')
    # Differential ring H=coth(xv), H_x=-v(H^2-1).
    xx,v,H=X(0),X(1),X(2)
    def deriv(p):return add(diff(p,0),sc(mul(v,add(sq(H),C(-1)),diff(p,2)),-1))
    second=deriv(deriv(mul(xx,H)))
    expected=sc(mul(v,add(sq(H),C(-1)),add(mul(xx,v,H),C(-1))),2)
    need(second==expected,'Full coth convexity identity')
    # a=1: F2 numerator equals D(2X-D), zero at X=D/2.
    need(not add(sub(F2num,2,C(1)),sc(mul(D,add(sc(x,2),sc(D,-1))),-1)),'a1 F2 factor')
    return 10

def K(t,r):return (t**r-t)/(t**r-1)
def candidates():
    total=positive=physical=endpoint=distinct_crossings=0
    for d in range(2,13):
        values=sorted(set([1,2,3,d*d,2*d*d,(d-1)**2+1]))
        for a in values:
            D=a+d
            for j in range(1,7):
                t=1+Q(1,2**j*a);total+=1
                kd=K(t,D)
                if a==1:
                    ka=Q(0);x=Q(D,2);R=Q(d)/kd
                    need(1<x<D,'a1 strict strip')
                else:
                    ka=K(t,a);den=D*kd-a*ka
                    need(kd>ka>0 and den>0,'Positive K and denominator')
                    x=(D*D*kd-a*a*ka)/(2*den);R=Q(a*D*d)/den
                    need(a+Q(d-1,2)<x<Q(D+a,2),'Exact strict X strip')
                Cc=a*a-2*a*x;Ee=D*D-2*D*x
                need((t*R+Cc)/(R+Cc)==t**D,'First complete modulus identity')
                need((t*R+Ee)/(R+Ee)==t**a,'Second complete modulus identity')
                y2=R-x*x;y1=t*R-x*x
                if y2<=0:continue
                positive+=1
                c1=(x*x-D*(x-a)**2)/(D-1)
                need(y2<c1<y1,'F1 strict crossing')
                if a>1:
                    c2=(x*x-a*(D-x)**2)/(a-1)
                    need(y2<c2<y1,'F2 strict crossing')
                    distinct_crossings+=int(c1!=c2)
                else:
                    need(Ee==0 and x*x-(D-x)**2==0,'a1 F2 identically zero');endpoint+=1
                # Sufficient physical wedge certificate: atan(y/X)<y/X<1/d<pi/d.
                if y1*d*d<x*x:physical+=1
    need(positive>100 and physical>30 and endpoint>0 and distinct_crossings>0,'Nonvacuous exact candidate domain')
    return dict(triples=total,positive_squared_ordinate_candidates=positive,
                rationally_certified_physical_candidates=physical,a_one_candidates=endpoint,
                distinct_crossing_pairs=distinct_crossings)
def hyperbolic_series():
    for n in range(1,17):
        coefficient=Q(1,factorial(2*n))-Q(1,factorial(2*n+1))
        need(coefficient==Q(2*n,factorial(2*n+1))>0,'x cosh x-sinh x positive series')
    return 16
def negatives():
    t=Q(5,4);a,d=2,3;D=a+d;kd,ka=K(t,D),K(t,a)
    den=D*kd-a*ka;x=(D*D*kd-a*a*ka)/(2*den)
    need(x>a+Q(d-1,2),'Wrong reversed X bound rejected')
    x0,y0,a0,D0=Q(4),Q(2),Q(2),Q(6)
    f1=D0/(x0*x0+y0*y0)-1/((x0-a0)**2+y0*y0)
    f2=1/((D0-x0)**2+y0*y0)-a0/(x0*x0+y0*y0)
    actual=-(D0-a0)*x0/(x0*x0+y0*y0)+(D0-x0)/((D0-x0)**2+y0*y0)+(x0-a0)/((x0-a0)**2+y0*y0)
    need(actual==(D0-x0)*f2-(x0-a0)*f1 and actual!=(D0-x0)*f2+(x0-a0)*f1,'Wrong theta sign rejected')
    need(K(t,1)==0,'a1 cannot divide by K_a')
    return 3
def source_controls(ledger,repo):
    names=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError,TypeError,OSError):names.append(name);return
        raise ValueError('Negative accepted:'+name)
    for role in ROLES:
        z=deepcopy(ledger);z['sources']=[r for r in z['sources'] if r['role']!=role]
        reject('omit:'+role,lambda z=z:source_ok(z,repo))
        z=deepcopy(ledger);next(r for r in z['sources'] if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda z=z:source_ok(z,repo))
    for key in ['global_fixed_sign_claimed','global_branch_labeling_claimed','activation_claimed','global_SFC_claimed','cutoff_claimed']:
        z=deepcopy(SCOPE);z[key]=True
        reject('false-scope:'+key,lambda z=z:input_ok(z))
    for v in [True,'1',2]:
        z=deepcopy(ledger);z['schema']=v
        reject('schema:'+str(v),lambda z=z:source_ok(z,repo))
    return names
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=HERE/RN);args=ap.parse_args();repo=args.repo_root.resolve()
    data=load(HERE/IN);input_ok(data);ledger=load(HERE/LN);source_ok(ledger,repo)
    identities=algebra();cases=candidates();series=hyperbolic_series();neg=negatives();controls=source_controls(ledger,repo)
    out={'schema':1,'status':'PASS','scope':data,'proof_sha256':sha(HERE/PN),'checker_sha256':sha(HERE/CN),
      'inputs_sha256':sha(HERE/IN),'source_ledger_sha256':sha(HERE/LN),'sources':ledger['sources'],
      'counts':{'universal_polynomial_differential_identities':identities,'positive_hyperbolic_series_coefficients':series,
       **cases,'mathematical_negative_controls':neg,'source_scope_controls':len(controls)},
      'negative_names':controls,'meaning':'Exact controls support the elementary universal strict-inequality proof; no activation/global sign/SFC consequence is inferred.',
      'assurance_flags':{'formal_verification':False,'external_review':False,'independent_replay':False}}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','counts':out['counts']},sort_keys=True))
if __name__=='__main__':main()
