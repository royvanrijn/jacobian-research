#!/usr/bin/env sage-python
"""Separate exact replay of the finite marked bisection experiment.

Does not import the constructor or enumerate another shell. It verifies the
marked lattice, group-law traces, line elimination, coefficient maps, and
rational nonsplitting. No exceptional point data are read.
"""
import hashlib
import json
from math import isqrt, prod
from pathlib import Path
import resource
import time
from sage.all import QQ, ZZ, PolynomialRing, matrix, vector, block_diagonal_matrix, prime_range

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/local/elliptic-curves/marked-two-class-v1'


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def add(P,Q,A):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    slope=(3*x*x+A)/(2*y) if x==u else (v-y)/(u-x)
    xx=slope*slope-x-u
    return xx,-y+slope*(x-xx)


def mul(n,P,A):
    if n<0:return mul(-n,(P[0],-P[1]),A)
    result=None
    while n:
        if n&1:result=add(result,P,A)
        n//=2
        if n:P=add(P,P,A)
    return result


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(1200,1205))
    start=time.process_time()
    g=json.loads((OUT/'geometry.json').read_text());proto=json.loads((OUT/'protocol.json').read_text())
    assert digest(OUT/'geometry.json')==proto['geometry_sha256']
    for rel,h in proto['sources'].items():assert digest(OUT/'source-snapshots'/rel)==h
    ns=matrix(ZZ,g['source_ns_gram']);L=matrix(ZZ,g['frame_gram'])
    transport=matrix(ZZ,g['unimodular_transport']);inverse=matrix(ZZ,g['inverse_transport'])
    assert transport*inverse==matrix.identity(ZZ,19) and abs(transport.det())==1
    assert transport*ns*transport.transpose()==block_diagonal_matrix(matrix(ZZ,[[0,1],[1,0]]),-L)
    assert L.det()==948 and L.is_positive_definite()
    F,O,Rc=[vector(ZZ,g[k]) for k in ['fibre','zero','nonidentity_I2_component']]
    assert F*ns*F==0 and O*ns*O==-2 and F*ns*O==1
    assert Rc*ns*Rc==-2 and Rc*ns*F==Rc*ns*O==0
    S=[vector(ZZ,s) for s in g['generic_divisors']]
    projections=[s-O-(2+s*ns*O)*F+QQ(s*ns*Rc)/2*Rc for s in S]
    assert -matrix(QQ,projections)*ns*matrix(QQ,projections).transpose()==matrix(QQ,g['family']['generic_height_gram'])
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen();family=g['family']
    A,B=[R(family[k]) for k in ['A_coefficients_low_to_high','B_coefficients_low_to_high']]
    assert (K(1728*4*A**3)/(4*A**3+27*B**2)).derivative()!=0
    def fun(r):return K(R(r['numerator']))/R(r['denominator'])
    def generic_fun(r):return K(R(r['numerator_coefficients_low_to_high']))/R(r['denominator_coefficients_low_to_high'])
    points=[(generic_fun(s['X']),generic_fun(s['Y'])) for s in family['sections']]
    for x,y in points:assert y*y==x**3+A*x+B
    trace_cache={};checked=[];constant_residuals=[];nonsplit=[];unknown=[]
    receipts=json.loads((OUT/'compilation-v2.json').read_text())['receipts']
    candidates=json.loads((OUT/'frozen-candidates.json').read_text())['rows']
    for receipt in receipts:
        index=receipt['index'];path=OUT/'compiled-export-v2'/f'{index:03d}.json';r=json.loads(path.read_text())
        if 'sha256' in receipt:assert digest(path)==receipt['sha256']
        assert r['divisor']==candidates[index]
        D=vector(ZZ,r['divisor']['native_ns_coordinates']);h=(r['divisor']['norm']-10)//4
        w=vector(ZZ,r['divisor']['reduced_frame_coordinates'])
        assert w*L*w==r['divisor']['norm']
        assert D==2*O+(h+4)*F+w*matrix(ZZ,g['frame_basis_in_source_ns'])
        assert D*ns*D==-2 and D*ns*F==2 and D*ns*O==h
        assert all(D*ns*vector(ZZ,s)>=0 for s in g['effective_wall_divisors'])
        if r['status']=='UNKNOWN_INTERRUPTED_BUDGET_CHARGED':unknown.append(index);continue
        word=tuple(r['trace_word'])
        if word not in trace_cache:
            P=None
            for n,G in zip(word,points):P=add(P,mul(int(n),G,A),A)
            trace_cache[word]=P
        trace=trace_cache[word];assert trace is not None
        tx,ty=trace[0],-trace[1]
        f0,f1,f2=map(R,r['RR_relation']);k=r['RR_vertical_degree']
        assert all(f.degree()<=b for f,b in zip([f0,f1,f2],[k,k-4,k-6]))
        assert f0+f1*tx+f2*ty==0 and f2!=0
        c=r['trace_zero_intersection'];m=r['I2_intersection'];eps=m%2
        pole=max(tx.denominator().degree()//2,
                 (tx.numerator().degree()-tx.denominator().degree()-4+1)//2,
                 (ty.numerator().degree()-ty.denominator().degree()-6+2)//3)
        assert pole==c and D*ns*Rc==m and 0<=m<=2
        phi=sum((n*p for n,p in zip(word,projections)),vector(QQ,[0]*19))
        negative=O+(c+2)*F-phi-QQ(eps)/2*Rc
        assert D+negative==3*O+k*F-r['subtract_I2_multiplicity']*Rc
        if r['subtract_I2_multiplicity']:
            location=r['I2_location'];nx=QQ(location['node_x'])
            if location['chart']=='finite':
                bad=QQ(location['t']);assert f0(bad)+f1(bad)*nx==0
                assert 3*nx**2+A(bad)==0 and nx**3+A(bad)*nx+B(bad)==0
            else:
                assert f0[k]+f1[k-4]*nx==0
                assert 3*nx**2+A[8]==0 and nx**3+A[8]*nx+B[12]==0
        X=PolynomialRing(K,'x');x=X.gen()
        cubic=(K(f0)+K(f1)*x)**2-K(f2)**2*(x**3+K(A)*x+K(B))
        residual,remainder=cubic.quo_rem(x-tx)
        assert remainder==0 and residual.degree()==2
        residual=residual.monic();disc=K(residual[1]**2-4*residual[0])
        if r['status']=='RESIDUAL_NOT_GEOMETRIC_GENUS_ZERO':
            # In this actual prefix all73 discarded branches are rational
            # constant squares; the residual splits over Q(t).
            assert len(r['branch'])==1 and QQ(r['branch'][0]).is_square()
            assert disc.is_square()
            constant_residuals.append(index);continue
        assert r['status']=='EXACT_GEOMETRIC_GENUS_ZERO_BISECTION'
        branch=R(r['branch']);assert branch.degree() in [1,2]
        assert branch.gcd(branch.derivative()).degree()==0
        sf=fun(r['discriminant_square_factor']);assert disc==sf*sf*branch
        x0,x1,y0,y1=[fun(r['point_map'][key]) for key in ['x0','x1','y0','y1']]
        assert 2*x0==-residual[1] and 2*x1==sf and x1!=0
        assert y0*y0+y1*y1*branch==x0**3+3*x0*x1*x1*branch+A*x0+B
        assert 2*y0*y1==3*x0*x0*x1+x1**3*branch+A*x1
        assert f0+f1*x0+f2*y0==0 and f1*x1+f2*y1==0
        value=QQ(branch(QQ(3)/17));assert value!=0
        # No floating approximation or factorization decides rational splitting.
        num,den=int(value.numerator()),int(value.denominator())
        a,b=isqrt(abs(num)),isqrt(den)
        assert num<0 or a*a!=num or b*b!=den
        witness=None
        for p in prime_range(3,1010):
            p=int(p)
            if p==17 or any(int(q.denominator())%p==0 for q in branch):continue
            xx=3*pow(17,-1,p)%p
            vv=sum((int(q.numerator())*pow(int(q.denominator()),-1,p)*pow(xx,j,p) for j,q in enumerate(branch)))%p
            if pow(vv,(p-1)//2,p)==p-1:
                witness={'prime':p,'control_residue':xx,'branch_value_residue':vv};break
        assert witness is not None
        nonsplit.append({'index':index,'branch_value_at_control':str(value),
            'absolute_numerator_floor_sqrt':str(a),'denominator_floor_sqrt':str(b),
            'local_nonsplitting':witness})
        checked.append(index)
    primes=sorted({r['local_nonsplitting']['prime'] for r in nonsplit})
    result={'schema':'marked-two-class.independent-replay.v1','status':'PASS_BOUNDED_VISIBILITY_MISS',
        'marking':'full integral rank17 frame, determinant948, exact U and A1 transports',
        'nonconstant_j':True,'generic_sections_checked':16,'candidate_divisors_checked':len(candidates),
        'exact_smooth_geometrically_rational_bisections':checked,
        'geometrically_reducible_residuals':constant_residuals,'unknown_candidates':unknown,
        'nonsplit_at_control':nonsplit,'nonzero_control_labels':0,
        'common_base':'NOT_REACHED_NO_SPLIT_CARRIER','prospective_fibres':'NOT_REACHED',
        'strict_quotient_increment':'UNKNOWN_NO_FRESH_FIBRES','ideal_class_increment':'UNKNOWN_NO_FRESH_FIBRES',
        'nonsplitting_progression':{'primes':primes,'M':str(prod(primes)),
            'formula':'t=3/17+M*z, z in Z','scope':'All182 completed bisections are nonsplit at every parameter in this explicit progression. This is a visibility exclusion for these curves only, not a rank or class-block exclusion.'},
        'cpu_seconds':time.process_time()-start,
        'inputs':{str(p.relative_to(ROOT)):digest(p) for p in [OUT/'geometry.json',OUT/'protocol.json',OUT/'frozen-candidates.json',OUT/'compilation-v2.json',Path(__file__)]},
        'boundary':'Checks all256 retained divisor classes and independently replays all255 completed equation computations. Does not regenerate the incomplete shell traversal or the upstream surface realization; the pre-existing exact model/marking certificates are inherited.'}
    output=OUT/'independent-replay.json'
    with output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','nonsplit_at_control','exact_smooth_geometrically_rational_bisections','geometrically_reducible_residuals']}),flush=True)


if __name__=='__main__':main()
