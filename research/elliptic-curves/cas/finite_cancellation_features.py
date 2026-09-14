#!/usr/bin/env sage -python
"""Target-blind pointed-model features and separately callable target evaluation.

No Bezout divisor, C/L optimization, factorization, or elliptic point search.
Residue leaves are proved constant or explicitly censored. The finite panel
cannot account for cancellation at unprocessed primes.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import gzip
import json
from math import gcd, isqrt, lcm, log2
from pathlib import Path
import signal
import time

from finite_cancellation_corpus import ROOT, OUT, canonical, digest, write
from half_lattice_pointed_sieve import invariants
from search_observability import transform, multiply, primitive
from pointed_box_equivalence import box_key

CAS=Path(__file__).resolve().parent
PRIMES=(2,3,5,7,11,13,17,19,23,29,31)

def ev(c,m,n=1): return sum(a*m**i*n**(4-i) for i,a in enumerate(c))
def integral_pair(n,d):
    values=tuple(map(F,n))+tuple(map(F,d));den=lcm(*(x.denominator for x in values))
    z=[int(x*den) for x in values];g=gcd(*z)
    return z[:5] if g==1 else [x//g for x in z[:5]],z[5:] if g==1 else [x//g for x in z[5:]]
def xmap(A,B,anchor,mapping):
    a,b=map(F,anchor);raw=[a**3+4*B,4*a*b,6*a*a+4*A,4*b,a]
    matrix=tuple(map(F,mapping['matrix']))
    n,d=integral_pair(transform(raw,matrix),transform(tuple(map(F,mapping['raw_coefficients'])),matrix))
    q=list(map(int,map(F,mapping['discriminant_quartic'])))
    assert all(F(v).denominator==1 for v in mapping['discriminant_quartic'])
    return n,d,q
def vp(z,p):
    if not z:return 10**9
    v=0
    while z%p==0:z//=p;v+=1
    return v
def cval(z,p): return min(vp(v,p) for v in z)
def fixed_val(z,p,digits=1):
    v=vp(z[0],p)
    return v if v<10**9 and all(vp(a,p)>=v+digits for a in z[1:]) else None
def q_condition(q,p):
    v=fixed_val(q,p)
    if v is not None and v%2:return 'nonsquare',v
    digits=3 if p==2 else 1
    v=fixed_val(q,p,digits)
    if v is None:return 'unknown',None
    unit=q[0]//p**v
    square=unit%8==1 if p==2 else pow(unit%p,(p-1)//2,p)==1
    return ('square' if square else 'nonsquare'),v
def affine(c,p,r):
    # Tiny-degree integer substitution, without Fraction arithmetic overhead.
    from math import comb
    return [p**j*sum(c[i]*comb(i,j)*r**(i-j) for i in range(j,5)) for j in range(5)]

def local_tree(n,d,q,p,depth=None,cap=4096):
    depth=depth or (8 if p==2 else 4);leaves=[];nodes=0
    def visit(nn,dd,qq,kind,a,k):
        nonlocal nodes
        nodes+=1
        status,vq=q_condition(qq,p);v=min(cval(nn,p),cval(dd,p))
        vg_fixed=(fixed_val(nn,p)==v or fixed_val(dd,p)==v)
        resolved=status=='nonsquare' or (status=='square' and vg_fixed)
        if resolved or k>=depth or nodes>=cap:
            leaves.append({'kind':kind,'residue':a,'depth':k,'mass':str(F(1,(p+1)*p**(k-1))),
                'status':status if resolved else 'unknown','vg':v,'vq':vq,
                'vg_exact':vg_fixed,'censored':not resolved})
            return
        for r in range(p):visit(affine(nn,p,r),affine(dd,p,r),affine(qq,p,r),kind,a+r*p**k,k+1)
    for r in range(p):visit(affine(n,p,r),affine(d,p,r),affine(q,p,r),'affine',r,1)
    visit(affine(n[::-1],p,0),affine(d[::-1],p,0),affine(q[::-1],p,0),'infinity',0,1)
    assert sum(F(x['mass']) for x in leaves)==1
    good=[x for x in leaves if x['status']=='square']
    unknown=sum(float(F(x['mass'])) for x in leaves if x['status']=='unknown')
    mass=sum(float(F(x['mass'])) for x in good)
    weights=[float(F(x['mass']))*p**(x['vq']/2) for x in good]
    denom=sum(weights)
    return {'prime':p,'nodes':nodes,'leaves':leaves,'soluble_mass':mass,'unknown_mass':unknown,
        'expected_vg_uniform':sum(float(F(x['mass']))*x['vg'] for x in good)/mass if mass else None,
        'expected_vg_differential':sum(w*x['vg'] for w,x in zip(weights,good))/denom if denom else None,
        'differential_resolved_mass':denom}

def features(n,d,q):
    locals=[local_tree(n,d,q,p) for p in PRIMES]
    samples=[];seen=set()
    for m,n0 in [(i,32) for i in range(-32,33)]+[(32,i) for i in range(-32,33)]:
        a,b=primitive(m,n0)
        if (a,b) in seen:continue
        seen.add((a,b))
        if ev(q,a,b)>=0:
            samples.append(log2(max(abs(ev(n,a,b)),abs(ev(d,a,b))))-4*log2(max(abs(a),abs(b))))
    if not samples:raise ArithmeticError('fixed real sample has no soluble address')
    f={'real_logS':sum(samples)/len(samples),'coefficient_bits':log2(max(map(abs,n+d))),
        'real_sample_count':len(samples),'unknown_mass_sum':sum(z['unknown_mass'] for z in locals),
        'log_soluble_mass':sum(log2(max(z['soluble_mass'],1e-15)) for z in locals),
        'local_g_uniform':sum((z['expected_vg_uniform'] or 0)*log2(z['prime']) for z in locals),
        'local_g_differential':sum((z['expected_vg_differential'] or 0)*log2(z['prime']) for z in locals),
        'missing_local_expectations':sum(z['expected_vg_uniform'] is None for z in locals)}
    return f,locals

def neighbour_maps(old,pari):
    """Same integral neighbour formula as retained pilot; no Curve302 constant."""
    P=list(map(F,old['reduced_P']));Q=list(map(F,old['reduced_Q']))+[F(0)]*2
    disc=list(map(F,old['discriminant_quartic']));inv=invariants(disc)
    seen={tuple(map(str,box_key(old['matrix'])))};out=[];attempts=[]
    polynomial=lambda f:'+'.join(f'({v})*x^{i}' for i,v in enumerate(f))
    for p in PRIMES[1:]:
        for r in range(p):
            if ev(disc,r)%p or sum(i*disc[i]*r**(i-1) for i in range(1,5))%p:continue
            s=(-int(ev(Q,r))*pow(2,-1,p))%p
            pp=affine(P,p,r);qq=affine(Q,p,r)
            pp=[(v-s*w-(s*s if i==0 else 0))/p**2 for i,(v,w) in enumerate(zip(pp,qq))]
            qq=[(v+(2*s if i==0 else 0))/p for i,v in enumerate(qq)]
            if any(v.denominator!=1 for v in pp+qq):
                attempts.append([p,r,'nonintegral']);continue
            ret=pari('my(m,C);C=hyperellred(['+polynomial(pp)+','+polynomial(qq[:3])+'],&m);[C,m]')
            RP=[F(str(ret[0][0].polcoef(i))) for i in range(5)]
            RQ=[F(str(ret[0][1].polcoef(i))) for i in range(3)]
            q=[4*RP[i]+sum(RQ[j]*RQ[i-j] for j in range(3) if 0<=i-j<3) for i in range(5)]
            change=tuple(F(str(ret[1][1][i,j])) for i in range(2) for j in range(2))
            second=multiply((p,r,0,1),change);matrix=multiply(tuple(map(F,old['matrix'])),second)
            key=tuple(map(str,box_key(matrix)))
            if key in seen:attempts.append([p,r,'duplicate']);continue
            tr=transform(tuple(map(F,old['raw_coefficients'])),matrix);j=next(i for i,z in enumerate(q) if z);ratio=tr[j]/q[j]
            assert inv==invariants(q) and all(a==ratio*b for a,b in zip(tr,q)) and ratio>0
            assert isqrt(ratio.numerator)**2==ratio.numerator and isqrt(ratio.denominator)**2==ratio.denominator
            m=deepcopy(old);m.update(matrix=list(map(str,matrix)),first_matrix=old['matrix'],second_matrix=list(map(str,second)),
                reduced_P=list(map(str,RP)),reduced_Q=list(map(str,RQ)),discriminant_quartic=list(map(str,q)),
                square_ratio=str(ratio),coordinate_policy={'kind':'raw','matrix':list(map(str,matrix))})
            seen.add(key);out.append({'name':f'neighbour-p{p}-r{r}','mapping':m,'prime':p,'residue':r})
            attempts.append([p,r,'retained'])
            if len(out)==2:return out,attempts
    return out,attempts

def prepare(curve,generic,anchor_index,mapper):
    """This function cannot receive a target, label or target-derived height."""
    start=time.process_time();A,B=map(F,curve[3:]);points=[tuple(map(F,p)) for p in generic]
    centre={'representative':[int(j==anchor_index) for j in range(len(points))]}
    t=time.process_time();old=mapper.mapping(tuple(map(F,curve)),points,centre)
    base_seconds=time.process_time()-t
    others,attempts=neighbour_maps(old,mapper.pari)
    maps=[{'name':'factor_free','mapping':old}]+others
    for m in maps:
        n,d,q=xmap(A,B,points[anchor_index],m['mapping'])
        m['N']=list(map(str,n));m['D']=list(map(str,d));m['q']=list(map(str,q))
        m['features'],m['local']=features(n,d,q)
    return {'curve':curve,'anchor':generic[anchor_index],'anchor_index':anchor_index,
        'models':maps,'attempts':attempts,'base_map_cpu_seconds':base_seconds,'preparation_cpu_seconds':time.process_time()-start}

def add(P,Q,A):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    slope=(3*x*x+A)/(2*y) if P==Q else (v-y)/(u-x)
    xx=slope*slope-x-u
    return xx,slope*(x-xx)-y
def leaf_for(tree,m,n):
    p=tree['prime'];kind='affine' if n%p else 'infinity'
    for z in tree['leaves']:
        if z['kind']!=kind:continue
        mod=p**z['depth'];a=m*pow(n,-1,mod)%mod if kind=='affine' else n*pow(m,-1,mod)%mod
        if a==z['residue']:return z
    raise ArithmeticError('projective residue partition missed target')
def evaluate(prepared,targets):
    """Oracle stage: receives sealed models and target points, never changes them."""
    out=[];A=F(prepared['curve'][3]);a,b=map(F,prepared['anchor'])
    for ti,P in enumerate(targets):
        for sign in (1,-1):
            x,y=F(P[0]),sign*F(P[1]);assert x!=a
            slope=(y+b)/(x-a);R=add(add((x,y),(x,y),A),(a,-b),A);assert R is not None
            for mi,M in enumerate(prepared['models']):
                u,v,w,z=map(F,M['mapping']['matrix']);m,n=primitive(z*slope-v,-w*slope+u)
                H=max(abs(m),abs(n));N=ev(list(map(int,M['N'])),m,n);D=ev(list(map(int,M['D'])),m,n)
                q=ev(list(map(int,M['q'])),m,n);g=gcd(N,D)
                assert D and F(N,D)==R[0] and q>=0 and isqrt(q)**2==q
                Hx=max(abs(R[0].numerator),R[0].denominator)
                assert max(abs(N),abs(D))==g*Hx
                primes=[];rest=g
                for tree in M['local']:
                    p=tree['prime'];val=vp(g,p);rest//=p**val;leaf=leaf_for(tree,m,n)
                    assert leaf['status']!='nonsquare'
                    if leaf['vg_exact']:assert leaf['vg']==val
                    primes.append({'p':p,'vg':val,'leaf':tree['leaves'].index(leaf)})
                S=F(max(abs(N),abs(D)),H**4)
                out.append({'target_index':ti,'sign':sign,'model_index':mi,'coordinate':[str(m),str(n)],
                    'H':str(H),'Hx':str(Hx),'g':str(g),'S':str(S),'logH':log2(H),'logg':log2(g),
                    'logS':log2(S.numerator)-log2(S.denominator),'local':primes,
                    'unprocessed_g_cofactor':str(rest),'cofactor_log2':log2(rest)})
    return out

def alarm(signum,frame):raise TimeoutError('declared per-anchor wall cap')
def build(limit=None):
    mapper=SourceFileLoader('finite_cancellation_factor_free',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    corpus=json.loads(gzip.decompress((OUT/'corpus.json.gz').read_bytes()))
    ph=digest((OUT/'protocol.json').read_bytes());ch=digest((OUT/'corpus.json.gz').read_bytes())
    source_hashes={p.name:digest(p.read_bytes()) for p in [Path(__file__),CAS/'lean_factor_free_pari_mapping.sage',CAS/'half_lattice_pointed_sieve.py',CAS/'search_observability.py']}
    signal.signal(signal.SIGALRM,alarm);start=time.process_time();done=0
    for case in corpus[:limit]:
        for ai in range(2):
            path=OUT/'cases'/f'{case["id"]}-{ai}.json.gz'
            if path.exists():continue
            before=time.process_time();signal.alarm(10)
            try:
                prepared=prepare(case['curve'],case['generic_points'],ai,mapper)
                # Serialize and hash before any target coordinate is consulted.
                sealed=canonical(prepared);seal=digest(sealed)
                measured=evaluate(prepared,case['targets'])
                payload={'status':'PASS_EXACT_ACCESSIBILITY','case_id':case['id'],'family':case['family'],
                    'j_group':case['j_group'],'prepared':prepared,'prepared_sha256':seal,'observations':measured}
            except (TimeoutError,ArithmeticError,ValueError,AssertionError) as e:
                payload={'status':'CENSORED_PREPARATION_OR_EVALUATION','case_id':case['id'],'anchor_index':ai,'reason':type(e).__name__+': '+str(e)}
            finally:signal.alarm(0)
            payload.update(protocol_sha256=ph,corpus_sha256=ch,source_sha256=source_hashes,cpu_seconds=time.process_time()-before)
            path.parent.mkdir(exist_ok=True,parents=True)
            with gzip.GzipFile(str(path),'wb',mtime=0) as f:f.write(canonical(payload))
            done+=1
            if done%100==0:print(json.dumps({'new_anchors':done,'case':case['id'],'cpu_seconds':time.process_time()-start}),flush=True)
            if time.process_time()-start>1800:print('DECLARED_CPU_CAP_CHECKPOINT',flush=True);return
    print(json.dumps({'status':'COMPLETE_REQUESTED_CORPUS','new_anchors':done,'cpu_seconds':time.process_time()-start}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--limit',type=int)
    a=p.parse_args();build(a.limit)
