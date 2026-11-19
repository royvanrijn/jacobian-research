#!/usr/bin/env python3
"""Exact standard-library replay of the low-degree odd-cubic certificates.

The separate proof supplies the unrestricted Galois reconstruction and
degree-case coverage. This checker derives the actual normalized moment
generators and replays exact rational Bezout identities without a CAS.
"""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def clean(p):return{k:Q(v)for k,v in p.items()if v}
def add(p,q,scale=1):
    p=dict(p)
    for k,v in q.items():p[k]=p.get(k,Q(0))+scale*v
    return clean(p)
def scale(p,v):return clean({k:v*c for k,c in p.items()})
def mul(p,q):
    out={}
    for a,v in p.items():
        for b,w in q.items():
            k=tuple(i+j for i,j in zip(a,b));out[k]=out.get(k,Q(0))+v*w
    return clean(out)
def power(p,n):
    out={(0,)*len(next(iter(p))):Q(1)}
    for _ in range(n):out=mul(out,p)
    return out
def evaluate(p,values):
    return sum((v*prod(a**k for a,k in zip(values,powers))for powers,v in p.items()),Q(0))
def prod(values):
    out=1
    for v in values:out*=v
    return out
def decode(rows):return clean({tuple(r['powers']):Q(r['coefficient'])for r in rows})
def decode_uni(rows):return clean({(r['power'],):Q(r['coefficient'])for r in rows})

def h01_moment(m,linear_H):
    out={};N=3*m+1
    for h in(range(m+1)if linear_H else[0]):
        for j in range(m-h+1):
            for k in range((m-h-j)//2+1):
                ell=m-h-j-2*k;mon=j,k,ell
                v=comb(m,h)*comb(N,j)*comb(N-j,k)*(-1)**ell*comb(m+ell,ell)
                out[mon]=out.get(mon,0)+v
    return clean(out)

def h3_moment(m,linear_J):
    out={}
    for j in range(m+1):
        for k in range((m-j)//2+1):
            for ell in range((m-j-2*k)//3+1):
                if j+k+ell>m:continue
                r=m-j-2*k-3*ell
                if not linear_J and r:continue
                v=comb(m,j)*comb(m-j,k)*comb(m-j-k,ell)*(comb(3*m+1,r)if linear_J else 1)
                mon=j,k,ell;out[mon]=out.get(mon,0)+v
    return clean(out)

def solve_c(f2):
    alpha=f2[(0,1,0)]
    assert all(j==0 or(i,j,k)==(0,1,0)for i,j,k in f2)
    return clean({(i,k):-v/alpha for(i,j,k),v in f2.items()if not j})
def substitute_c(p,cstar):
    out={}
    for(i,j,k),v in p.items():out=add(out,mul({(i,k):v},power(cstar,j)))
    return out
def set_d_one(p):
    out={}
    for(i,j),v in p.items():out[(i,)]=out.get((i,),Q(0))+v
    return clean(out)

def check_unit_certificate(case,fn):
    cstar=solve_c(fn(2));generators=[]
    for m,factor,encoded in zip(case['generator_moments'],case['generator_scalars_relative_to_substituted_moment'],case['primitive_generators']):
        got=scale(substitute_c(fn(m),cstar),Q(factor))
        assert got==decode(encoded);generators.append(got)
    multipliers=list(map(decode,case['multipliers']));total={}
    for a,g in zip(multipliers,generators):total=add(total,mul(a,g))
    assert total=={(0,0):Q(1)}
    altered=multipliers[0].copy();altered[(0,0)]=altered.get((0,0),Q(0))+1
    bad=mul(altered,generators[0])
    for a,g in zip(multipliers[1:],generators[1:]):bad=add(bad,mul(a,g))
    assert bad!={(0,0):Q(1)}
    return {'degree_bound':case['degree_bound'],'multiplier_terms':list(map(len,multipliers)),
            'identity':'exactly1','corrupted_multiplier_rejected':True}

def check_univariate_certificate(case,expected):
    gens=list(map(decode_uni,case['generators']));multipliers=list(map(decode_uni,case['multipliers']))
    assert gens==expected
    total={}
    for a,g in zip(multipliers,gens):total=add(total,mul(a,g))
    assert total=={(0,):Q(1)}
    return {'generator_degrees':[max(m[0]for m in g)for g in gens],
            'multiplier_degrees':[max(m[0]for m in a)for a in multipliers],'identity':'exactly1'}

def h0_controls(case):
    fn=lambda m:h01_moment(m,False)
    cstar=solve_c(fn(2));g4=substitute_c(fn(4),cstar);g6=substitute_c(fn(6),cstar)
    result=check_univariate_certificate(case,[set_d_one(g4),set_d_one(g6)])
    assert{(i,):v for(i,j),v in g4.items()if not j}=={(4,):-Q(1157)}
    result['d_zero_restriction']='-1157*b^4'
    return result

def constant_J_controls():
    fn=lambda m:h3_moment(m,False)
    cstar=solve_c(fn(2));g4=substitute_c(fn(4),cstar);g6=substitute_c(fn(6),cstar)
    assert g4=={(4,0):-Q(7,2),(1,1):Q(12)}
    assert{(j,):v for(i,j),v in g6.items()if not i}=={(2,):Q(15)}
    assert evaluate(g6,[Q(1),Q(7,24)])==Q(1397,192)
    return {'first_two_equations':'b=-a^2/2; a*(24c-7a^3)/2=0',
            'a_zero_sixth_equation':'15c^2=0','a_nonzero_normalized_sixth':str(Q(1397,192))}

def zero_center_controls(case):
    expected=[]
    for m in[2,4]:expected.append({(m-j,):Q(comb(m,j),m+2*j+1)for j in range(m+1)})
    return check_univariate_certificate(case,expected)

def unary_mul(a,b,order):
    out=[Q(0)]*min(len(a)+len(b)-1,order+1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):
            if i+j<len(out):out[i+j]+=v*w
    return out
def unary_power(a,n,order):
    out=[Q(1)]
    for _ in range(n):out=unary_mul(out,a,order)
    return out
def raw_even_moment(C,D,m):
    return sum((Q(comb(m,j),m+2*j+1)*(
        unary_mul(unary_power(C,m-j,m),unary_power(D,j,m),m)+[Q(0)]*(m+1))[m]
        for j in range(m+1)),Q(0))
def kappa(m):
    if m%2:return Q(0)
    return sum((Q((-1)**j*comb(m,j),m+2*j+1)for j in range(m+1)),Q(0))
def raw_controls():
    count=0
    for linearH in[False,True]:
        for b,c,d in[(0,0,0),(1,2,3),(-2,1,2),(Q(1,2),-1,Q(2,3))]:
            H=[1,1]if linearH else[1];J=[1,b,c];K=[1,d]
            C=unary_mul(H,unary_power(J,2,5),5);D=[-v for v in unary_mul(H,unary_power(K,2,3),3)]
            for m in[2,4,6,8]:
                assert raw_even_moment(C,D,m)==kappa(m)*evaluate(h01_moment(m,linearH),[b,c,d]);count+=1
    for b,c,d in[(0,0,0),(1,2,3),(-2,1,2),(Q(1,2),-1,Q(2,3))]:
        H=[1,b,c,d];C=unary_mul(H,[1,2,1],5);D=[-v for v in H]
        for m in[2,4,6,8]:
            assert raw_even_moment(C,D,m)==kappa(m)*evaluate(h3_moment(m,True),[b,c,d]);count+=1
    return count

def bmul(a,b):
    out={}
    for(i,j),v in a.items():
        for(k,l),w in b.items():
            mon=i+k,j+l;out[mon]=out.get(mon,0)+v*w
    return{k:v for k,v in out.items()if v}
def badd(a,b,scale=1):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,0)+scale*v
    return{k:v for k,v in out.items()if v}
def bpower(a,n):
    out={(0,0):1}
    for _ in range(n):out=bmul(out,a)
    return out
def original_sharp_six():
    R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}
    factor=badd(bmul(bpower(R,2),W),{},1)
    factor={k:2*v for k,v in factor.items()}
    factor=badd(factor,bmul(bpower(T,2),badd({k:2*v for k,v in R.items()},Z)),-1)
    quotient=bmul(bpower(badd(R,Z),2),factor);F=bmul(T,quotient)
    assert sum(quotient.values())==144
    assert sum(v*(-1)**(6-j)for(i,j),v in F.items())==8
    P={(0,0):1};rows=[]
    for m in range(1,9):
        P=bmul(P,F);n=6*m
        pure=sum(v*factorial(i)*factorial(n-i)for(i,j),v in P.items()if j==i)
        mixed=sum(v*factorial(i+1)*factorial(n-i)for(i,j),v in P.items()if j==i+1)
        assert pure==0 and Q(mixed,factorial(n+2))==kappa(m)
        rows.append({'m':m,'support':len(P),'pure':0,'normalized_marker_Z':str(kappa(m))})
    return {'degree':6,'radial_zero_value':8,'T_quotient_value_at_T_zero':144,'controls':rows}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path)
    args=parser.parse_args()
    path=HERE/'odd_cubic_height_low_degree_certificates.json';data=json.loads(path.read_text())
    pins=[HERE/'ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md',Path(__file__),path,
          HERE/'odd_cubic_height_low_degree_source_ledger.json']
    ledger=json.loads(pins[-1].read_text())
    for resource in ledger['retained_resources']:
        p=HERE/resource['retained_path']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==resource['sha256']
        pins.append(p)
    out={'status':'passed','python':platform.python_version(),
         'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in pins},
         'h1_unit_ideal':check_unit_certificate(data['h1'],lambda m:h01_moment(m,True)),
         'h3_J_nonconstant_unit_ideal':check_unit_certificate(data['h3_j1'],lambda m:h3_moment(m,True)),
         'h0_constant_squarefree_part':h0_controls(data['h0_d_nonzero']),
         'constant_J_nonconstant_H':constant_J_controls(),
         'zero_center':zero_center_controls(data['zero_center_core']),
         'independent_raw_interval_controls':raw_controls(),'sharp_degree_six':original_sharp_six(),
         'certificates_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'scope':'Exact finite ideal-one identity and low-degree case controls, independent raw interval and original four-variable contractions. Full rational-u reduction and universal classification implications use the separate unrestricted written proof.'}
    if args.write_receipt:args.write_receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'raw_controls':out['independent_raw_interval_controls'],
                      'unit_ideal_cases':2,'original_F6_controls':8,'all_low_degree_branches_pass':True}))

if __name__=='__main__':main()
