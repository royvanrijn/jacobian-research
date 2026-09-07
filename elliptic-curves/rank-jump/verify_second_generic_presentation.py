#!/usr/bin/env python3
"""Rational coefficient replay, root enumeration and global duplicate identities."""
import argparse
from fractions import Fraction as Q
from math import isqrt,gcd,lcm
from pathlib import Path
import retrospective as r
import second_generic_presentation as source
from verify_fresh_symbolic_discriminant import trim,mul,value

OUTPUT=r.OUT/'rank_jump_second_generic_presentation_verification_v1.json'


def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([x*c for x in a])
def power(a,n):
    out=[Q(1)]
    for _ in range(n):out=mul(out,a)
    return out
def prime(p):return p>1 and all(p%i for i in range(2,isqrt(p)+1))
def mod(a,p):return [x.numerator*pow(x.denominator,-1,p)%p for x in a]
def roots_mod(a,p):
    assert prime(p);v=mod(a,p);assert v[-1]
    roots=[]
    for x in range(p):
        y=0
        for c in reversed(v):y=(y*x+c)%p
        if not y:roots.append(x)
    return v,roots
def root_free(a):
    for p in range(3,252,2):
        if not prime(p) or any(x.denominator%p==0 for x in a):continue
        v,roots=roots_mod(a,p)
        if not roots:return {'prime':p,'polynomial_mod_p':v}
    raise AssertionError('No bounded independent exclusion witness')
def square_root(q):
    if q<0:return None
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    return Q(a,b) if a*a==q.numerator and b*b==q.denominator else None


def homogeneous_substitute(f,n,d,degree):
    out=[Q(0)]
    for i,c in enumerate(f):out=add(out,scale(mul(power(n,i),power(d,degree-i)),c))
    return out


def duplicates(inputs,rows):
    from sage.all import QQ,matrix
    families={x['family']:x for x in inputs['families']};cases={x['token']:x for x in inputs['cases']}
    pairs=[('074d9','historic-lineage-074d9'),('a1-fibration-01','curve398-p16875')];out=[]
    for left,right in pairs:
        matches=[]
        for token in cases:
            a=[x for x in rows if x['token']==token and x['family']==left][0]['preimages']
            b=[x for x in rows if x['token']==token and x['family']==right][0]['preimages']
            if not a or not b:continue
            assert len(a)==len(b)==1 and a[0]['rationally_isomorphic'] and b[0]['rationally_isomorphic']
            matches.append((Q(a[0]['parameter']),Q(b[0]['parameter'])))
        assert len(matches)>=3
        M=matrix(QQ,[[str(t),1,str(-t*u),str(-u)] for t,u in matches[:3]])
        kernel=M.right_kernel().basis();assert len(kernel)==1
        v=list(map(lambda x:Q(str(x)),kernel[0]));den=lcm(*(x.denominator for x in v));v=[int(x*den) for x in v]
        g=gcd(*v);v=[x//g for x in v]
        if next(x for x in v if x)<0:v=[-x for x in v]
        a,b,c,d=v;assert a*d-b*c
        assert all(c*t+d and (a*t+b)/(c*t+d)==u for t,u in matches)
        n=[Q(b),Q(a)];denom=[Q(d),Q(c)]
        A=list(map(Q,families[left]['A']));B=list(map(Q,families[left]['B']))
        AA=homogeneous_substitute(list(map(Q,families[right]['A'])),n,denom,8)
        BB=homogeneous_substitute(list(map(Q,families[right]['B'])),n,denom,12)
        ca=AA[-1]/A[-1];cb=BB[-1]/B[-1];u2=cb/ca;u=square_root(u2)
        assert u is not None and AA==scale(A,u**4) and BB==scale(B,u**6)
        out.append({'left_family':left,'right_family':right,'mobius_matrix':v,'weierstrass_numerator':str(u),
            'weierstrass_scaling':'numerator/(c*t+d)^2 from left to right',
            'A_identity_degree':8,'B_identity_degree':12,'specialized_matches':len(matches)})
    return out


def compute():
    data=r.read(source.INPUT);out=r.read(source.OUTPUT)
    for obj in (data,out):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha
    families={x['family']:x for x in data['families']};cases={x['token']:x for x in data['cases']};verified=[]
    assert len(out['rows'])==112 and len({(x['token'],x['family']) for x in out['rows']})==112
    for row in out['rows']:
        assert row['status']=='PASS';c=cases[row['token']];f=families[row['family']]
        A=list(map(Q,f['A']));B=list(map(Q,f['B']));a,b=Q(c['a']),Q(c['b'])
        F=add(scale(power(A,3),b*b),scale(power(B,2),-a**3))
        assert F==list(map(Q,row['j_equation_ascending'])) and len(F)-1==row['finite_degree']
        assert row['infinity_multiplicity']==25-len(F)
        remaining=list(map(Q,row['remaining_monic_polynomial']));product=scale(remaining,Q(row['remaining_unit']));found=[]
        control=row['control']
        if control:
            assert c['control_family']==f['family'] and control['parameter']==c['control_parameter']
            t=Q(control['parameter']);assert value(remaining,t)!=0
            product=mul(product,power([-t,Q(1)],control['multiplicity']));found.append((str(t),control['multiplicity']))
        assert F==product
        witnesses=[]
        if row['root_free_modular_certificate']:
            w=row['root_free_modular_certificate'];v,roots=roots_mod(remaining,w['prime'])
            assert v==w['polynomial_mod_p'] and roots==[];witnesses.append(w)
        elif row['rational_factorization']:
            product=[Q(1)]
            for fac in row['rational_factorization']:
                q=list(map(Q,fac['coefficients_ascending']));assert q[-1]==1 and len(q)-1==fac['degree']
                product=mul(product,power(q,fac['multiplicity']))
                if len(q)==2:found.append((str(-q[0]),fac['multiplicity']))
                else:witnesses.append(root_free(q))
            assert product==remaining
        else:assert remaining==[1]
        if len(F)<25:found.append(('infinity',25-len(F)))
        assert sorted(found)==sorted((x['parameter'],x['multiplicity']) for x in row['preimages'])
        for x in row['preimages']:
            av,bv=(A[8],B[12]) if x['parameter']=='infinity' else (value(A,Q(x['parameter'])),value(B,Q(x['parameter'])))
            assert [str(av),str(bv)]==[x['source_a'],x['source_b']] and 4*av**3+27*bv**2
            u2=(b/bv)/(a/av);assert str(u2)==x['scaling_square'] and u2**2==a/av and u2**3==b/bv
            u=square_root(u2);assert (u is not None)==x['rationally_isomorphic']
            assert x['rational_scaling']==(str(u) if u is not None else None)
        verified.append({'token':row['token'],'family':row['family'],'rational_preimages':len(found),
            'rational_isomorphisms':sum(x['rationally_isomorphic'] for x in row['preimages']),
            'root_exclusion_certificates':witnesses})
    duplicate=duplicates(data,out['rows'])
    classes={name:name for name in families}
    for d in duplicate:classes[d['right_family']]=d['left_family']
    for row in out['rows']:
        if row['preimages']:assert classes[row['family']]==classes[cases[row['token']]['control_family']]
    return {'schema':'rank-jump.second-generic-presentation-verification.v1','status':'PASS','rows':verified,
        'global_base_equivalences':duplicate,'target_family_cells':112,'rational_preimages':sum(x['rational_preimages'] for x in verified),
        'non_base_equivalent_rational_presentations':0,
        'bindings':source.bindings([Path(__file__),source.INPUT,source.OUTPUT,
            Path(__file__).with_name('verify_fresh_symbolic_discriminant.py')]),
        'boundary':'Every rational preimage excluded or verified by rational coefficients and direct finite-field root enumeration. Both duplicate families globally identified by exact weighted polynomial substitutions. No point or rank computation.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['target_family_cells'],'cells;',result['rational_preimages'],'rational preimages;',len(result['global_base_equivalences']),'global duplicate identities')
