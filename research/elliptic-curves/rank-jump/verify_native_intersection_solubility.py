#!/usr/bin/env python3
"""Independent exact group law and rational short-coset enumeration."""
import argparse
from math import isqrt
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_native_pair_collapse_locus_inputs_v1.json'
LOCUS=r.OUT/'rank_jump_native_pair_collapse_locus_v1.json'
NORM_INPUT=r.OUT/'rank_jump_norm_six_carrier_solubility_inputs_v1.json'
NORM=r.OUT/'rank_jump_norm_six_carrier_solubility_v1.json'
PRIOR=r.OUT/'rank_jump_paired_character_moments_v1.json'
PRIOR_VERIFIED=r.OUT/'rank_jump_paired_character_moments_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_native_intersection_solubility_verification_v1.json'


def add(A,P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    m=(v-y)/(u-x) if x!=u else (3*x*x+A)/(2*y)
    xx=m*m-x-u;return xx,m*(x-xx)-y


def short_vectors(G,bound):
    """Exact LDL enumeration: all pruning and endpoint rounding are rational."""
    n=G.nrows();L=matrix(QQ,n,n);ds=[]
    for i in range(n):
        L[i,i]=1;d=G[i,i]-sum(L[i,k]**2*ds[k] for k in range(i));assert d>0;ds.append(d)
        for j in range(i+1,n):L[j,i]=(G[j,i]-sum(L[j,k]*L[i,k]*ds[k] for k in range(i)))/d
    x=[0]*n;out=[];nodes=0
    def visit(i,left):
        nonlocal nodes
        nodes+=1
        if i<0:
            if any(x):
                v=vector(ZZ,x);assert v*G*v<=bound;out.append(tuple(x))
            return
        c=QQ(sum(L[j,i]*x[j] for j in range(i+1,n)));den=c.denominator();num=c.numerator()
        rad=isqrt(int((left*den**2/ds[i]).floor()))
        lo=int((QQ(-rad-num)/den).ceil());hi=int((QQ(rad-num)/den).floor())
        for z in range(lo,hi+1):
            used=ds[i]*(z+c)**2;assert used<=left
            x[i]=z;visit(i-1,left-used)
        x[i]=0
    visit(n-1,QQ(bound));return out,nodes


def verify():
    inp=r.read(INPUT);old=r.read(LOCUS);ni=r.read(NORM_INPUT);no=r.read(NORM)
    for data in (inp,old,ni,no):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');K=R.fraction_field();A=R(inp['A']);B=R(inp['B'])
    def dec(c):return K(R(c['numerator']))/R(c['denominator'])
    basis=[]
    for c in inp['sections']:
        x=K(R(c['x_coefficients_low_to_high']))
        if 'y_coefficients_low_to_high' in c:y=K(R(c['y_coefficients_low_to_high']))
        else:
            ch=c['chord'];ref=basis[ch['reference_basis_index']]
            y=ref[1]+R(ch['slope_coefficients_low_to_high'])*(x-ref[0])
        assert y*y==x*x*x+A*x+B;basis.append((x,y))
    S=None
    for n,P in zip(inp['generic_word'],basis,strict=True):
        n=int(n);assert abs(n)<=1
        if n:S=add(A,S,(P[0],n*P[1]))
    assert S==(dec(old['generic_translate']['x']),dec(old['generic_translate']['y']))
    c,d=inp['covers'];q=R(c['residual_chord']['q_coefficients']);q2=R(d['residual_chord']['q_coefficients'])
    Z=PolynomialRing(K,'z');L=K.extension(Z.gen()**2-q,'u');u=L.gen()
    def lift(c):return [R(c['lifted_section'][key+'_coefficients']) for key in ('x0','x1','y0','y1')]
    x0,x1,y0,y1=lift(c);P=(L(x0)+u*x1,L(y0)+u*y1);Q=add(L(A),(L(S[0]),L(S[1])),(P[0],-P[1]))
    h=R(d['trace_section']['h_coefficients']);m=K(R(d['residual_chord']['M_coefficients']))/h
    tx=K(R(d['trace_section']['Nx_coefficients']))/h**2;ty=K(R(d['trace_section']['Ny_coefficients']))/h**3
    b=-ty-m*tx;residual=Q[1]-m*Q[0]-b
    aa=dec(old['chord_residual']['constant']);bb=dec(old['chord_residual']['root_coefficient'])
    assert residual==aa+u*bb and bb
    maps={key:dec(c) for key,c in old['rational_maps'].items()};ur=maps['u'];assert ur==-aa/bb
    xx=K(Q[0][0])+ur*K(Q[0][1]);yy=K(Q[1][0])+ur*K(Q[1][1])
    dx0,dx1,dy0,dy1=lift(d)
    assert maps['qx']==xx and maps['qy']==yy and maps['v']==(xx-dx0)/dx1
    assert maps['px']==x0+ur*x1 and maps['py']==y0+ur*y1
    norm=aa*aa-q*bb*bb
    f=R(d['quadratic_cover']['leading_coefficients'])*xx**2+R(d['quadratic_cover']['linear_coefficients'])*xx+R(d['quadratic_cover']['constant_coefficients'])
    common=norm.numerator().gcd(f.numerator()).monic();good=R(old['intersection_polynomial']);ex=R(old['excluded_polynomial'])
    while common.gcd(ex).degree()>0:common=common//common.gcd(ex)
    assert common==good and good.degree()==1 and ex.gcd(good)==1
    root=-good[0]/good[1];assert root==QQ(-4112)/1937
    # Verify both native square equations and the group relation independently.
    P0=(maps['px'](root),maps['py'](root));Q0=(xx(root),yy(root));S0=(S[0](root),S[1](root))
    assert maps['u'](root)**2==q(root) and maps['v'](root)**2==q2(root)
    assert Q0==(dx0(root)+maps['v'](root)*dx1(root),dy0(root)+maps['v'](root)*dy1(root))
    assert add(A(root),P0,Q0)==S0
    assert all(y*y==x*x*x+A(root)*x+B(root) for x,y in (P0,Q0,S0))
    checks=[];G=matrix(ZZ,ni['gram']);prior=r.read(PRIOR)
    for case,row in zip(ni['pairs'],no['rows'],strict=True):
        assert row['id']==case['id'];V=matrix(ZZ,row['reduced_basis']);w=vector(ZZ,case['traces'][0])+vector(ZZ,case['traces'][1])
        assert abs(V.det())==2**16
        for v in V.rows():assert all(z%2==0 for z in v) or all((v[i]-w[i])%2==0 for i in range(17))
        vs,nodes=short_vectors(V*G*V.transpose(),6)
        actual={tuple(vector(ZZ,v)*V) for v in vs};expected=set()
        for v in row['unoriented_vectors']:
            z=vector(ZZ,v['vector']);assert z*G*z==v['norm']
            expected.add(tuple(z));expected.add(tuple(-z))
        assert actual==expected and len(actual)==row['signed_short_vector_count']
        inherited=next(x for x in prior['paired_carrier_comparisons'] if set(x['labels'])==set(case['labels']))
        assert inherited['carrier_global_solubility']==case['previously_certified_global_solubility']
        if row['id']=='FG':
            z=vector(ZZ,row['known_intersection_vector']);ss=vector(ZZ,[ZZ(x) for x in inp['generic_word']])
            assert z==2*ss-w and z*G*z==6 and row['intersection_number']==1
        checks.append({'id':row['id'],'signed_vectors':len(actual),'exact_enumeration_nodes':nodes,'has_norm_six':row['has_norm_six']})
    return {'schema':'rank-jump.native-intersection-solubility-verification.v1','status':'PASS',
        'unique_open_intersection_parameter':str(root),'rational_group_relation_verified':True,
        'coset_checks':checks,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,LOCUS,NORM_INPUT,NORM,PRIOR,PRIOR_VERIFIED,Path(__file__),HERE/'retrospective.py')},
        'boundary':'Independent rational group arithmetic and exact LDL enumeration; no PARI short-vector enumeration or numerical heights in this verifier.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=verify()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert r.read(OUTPUT)==data
    print('PASS unique rational translated intersection; independent norm-six coset checks',data['coset_checks'])
