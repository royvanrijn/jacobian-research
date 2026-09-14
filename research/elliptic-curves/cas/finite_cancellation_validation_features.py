#!/usr/bin/env python3
"""Pure q residue features with separately charged optional gcd refinement.

Both policies share exactly the same real and q features. The q tree has no
N,D inputs. Gcd constancy may refine soluble leaves but cannot alter q masses.
"""
from fractions import Fraction as F
from math import log2

from finite_cancellation_features import (
    PRIMES, affine, cval, ev, fixed_val, neighbour_maps, q_condition, xmap)
from search_observability import primitive

REAL = ['real_logS', 'coefficient_bits']
Q_NAMES = ['q_log_soluble_mass', 'q_unknown_mass_sum', 'q_missing_local_expectations']
G_NAMES = ['local_g_uniform', 'local_g_differential', 'g_unknown_mass_sum', 'g_missing_local_expectations']
FEATURES = {'real_q': REAL + Q_NAMES, 'real_q_g': REAL + Q_NAMES + G_NAMES}


def q_tree(q, p, depth=None, cap=4096):
    depth = depth or (8 if p==2 else 4)
    leaves=[]; nodes=0
    def visit(qq, kind, a, k):
        nonlocal nodes
        nodes += 1
        status,vq=q_condition(qq,p)
        if status!='unknown' or k>=depth or nodes>=cap:
            leaves.append({'kind':kind,'residue':a,'depth':k,
                           'mass':str(F(1,(p+1)*p**(k-1))), 'status':status,'vq':vq})
            return
        for r in range(p):visit(affine(qq,p,r),kind,a+r*p**k,k+1)
    for r in range(p):visit(affine(q,p,r),'affine',r,1)
    visit(affine(q[::-1],p,0),'infinity',0,1)
    assert sum(F(x['mass']) for x in leaves)==1
    good=[x for x in leaves if x['status']=='square']
    mass=sum(float(F(x['mass'])) for x in good)
    return {'prime':p,'nodes':nodes,'leaves':leaves,'soluble_mass':mass,
            'unknown_mass':sum(float(F(x['mass'])) for x in leaves if x['status']=='unknown')}


def q_features(q):
    trees=[q_tree(q,p) for p in PRIMES]
    return {'q_log_soluble_mass':sum(log2(max(t['soluble_mass'],1e-15)) for t in trees),
            'q_unknown_mass_sum':sum(t['unknown_mass'] for t in trees),
            'q_missing_local_expectations':sum(t['soluble_mass']==0 for t in trees)}, trees


def gcd_tree(n, d, qt, depth=None, cap=4096):
    """Refine only resolved soluble q leaves; retain all uncertainty explicitly."""
    p=qt['prime'];depth=depth or (8 if p==2 else 4);leaves=[];nodes=0
    def visit(nn,dd,kind,a,k,vq):
        nonlocal nodes
        nodes+=1
        vg=min(cval(nn,p),cval(dd,p))
        exact=fixed_val(nn,p)==vg or fixed_val(dd,p)==vg
        if exact or k>=depth or nodes>=cap:
            leaves.append({'kind':kind,'residue':a,'depth':k,'mass':str(F(1,(p+1)*p**(k-1))),
                           'status':'square' if exact else 'unknown','vq':vq,'vg':vg,'vg_exact':exact})
            return
        for r in range(p):visit(affine(nn,p,r),affine(dd,p,r),kind,a+r*p**k,k+1,vq)
    for leaf in qt['leaves']:
        if leaf['status']!='square':
            leaves.append({**leaf,'vg':None,'vg_exact':False});continue
        nn=n if leaf['kind']=='affine' else n[::-1]
        dd=d if leaf['kind']=='affine' else d[::-1]
        visit(affine(nn,p**leaf['depth'],leaf['residue']),affine(dd,p**leaf['depth'],leaf['residue']),
              leaf['kind'],leaf['residue'],leaf['depth'],leaf['vq'])
    assert sum(F(x['mass']) for x in leaves)==1
    good=[x for x in leaves if x['status']=='square']
    mass=sum(float(F(x['mass'])) for x in good)
    weights=[float(F(x['mass']))*p**(x['vq']/2) for x in good]
    denom=sum(weights)
    return {'prime':p,'nodes':nodes,'leaves':leaves,
            'unknown_mass':sum(float(F(x['mass'])) for x in leaves if x['status']=='unknown'),
            'expected_vg_uniform':sum(float(F(x['mass']))*x['vg'] for x in good)/mass if mass else None,
            'expected_vg_differential':sum(w*x['vg'] for w,x in zip(weights,good))/denom if denom else None}


def gcd_features(n,d,q_trees):
    trees=[gcd_tree(n,d,t) for t in q_trees]
    return {'local_g_uniform':sum((t['expected_vg_uniform'] or 0)*log2(t['prime']) for t in trees),
            'local_g_differential':sum((t['expected_vg_differential'] or 0)*log2(t['prime']) for t in trees),
            'g_unknown_mass_sum':sum(t['unknown_mass'] for t in trees),
            'g_missing_local_expectations':sum(t['expected_vg_uniform'] is None for t in trees)}, trees


def real_features(n,d,q):
    samples=[];seen=set()
    for m,n0 in [(i,32) for i in range(-32,33)]+[(32,i) for i in range(-32,33)]:
        a,b=primitive(m,n0)
        if (a,b) in seen:continue
        seen.add((a,b))
        if ev(q,a,b)>=0:
            samples.append(log2(max(abs(ev(n,a,b)),abs(ev(d,a,b))))-4*log2(max(abs(a),abs(b))))
    if not samples:raise ArithmeticError('fixed real sample has no soluble address')
    return {'real_logS':sum(samples)/len(samples),'coefficient_bits':log2(max(map(abs,n+d)))}


def model_features(n,d,q,with_g):
    qf,qt=q_features(q);f={**real_features(n,d,q),**qf};local={'q':qt}
    if with_g:
        gf,gt=gcd_features(n,d,qt);f.update(gf);local['g']=gt
    return f,local


def prepare(curve,anchor,mapper,with_g):
    curve=tuple(map(F,curve));point=tuple(map(F,anchor))
    assert curve[:3]==(0,0,0), 'The retained factor-free worker expects a short model.'
    old=mapper.mapping(curve,[point],{'representative':[1]})
    other,attempts=neighbour_maps(old,mapper.pari)
    models=[{'name':'factor_free','mapping':old}]+other
    for model in models:
        n,d,q=xmap(curve[3],curve[4],point,model['mapping'])
        model.update(N=list(map(str,n)),D=list(map(str,d)),q=list(map(str,q)))
        model['features'],model['local']=model_features(n,d,q,with_g)
    return {'curve':list(map(str,curve)),'anchor':list(map(str,point)),
            'models':models,'attempts':attempts}


def select(prepared,fit):
    base=prepared['models'][0]['features'];scores=[]
    for m in prepared['models']:
        z=[(m['features'][k]-base[k]-mu)/sd for k,mu,sd in zip(fit['features'],fit['means'],fit['scales'])]
        scores.append(fit['coefficients'][0]+sum(a*b for a,b in zip(fit['coefficients'][1:],z)))
    return min(range(len(scores)),key=lambda i:scores[i]),scores
