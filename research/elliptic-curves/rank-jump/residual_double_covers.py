#!/usr/bin/env python3
"""Frozen conic parametrizations and simultaneous residual quartic tests."""
import argparse
from itertools import product
from math import gcd,isqrt
from pathlib import Path
import subprocess
import retrospective as r
import small_quotient_covers as old
from retrospective_secant_pencils import rational_square

PROTOCOL=Path(__file__).with_name('RESIDUAL_DOUBLE_COVER_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_residual_double_covers_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-residual-double-covers-v1'
PAIRS=[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]


def evaluate(coeff,s,t):
    return sum(r.F(c)*s**(4-i)*t**i for i,c in enumerate(coeff))


def parameters(H):
    yield (1,0)
    for t in range(1,H+1):
        for s in range(-H,H+1):
            if gcd(s,t)==1:yield (s,t)


def local_test(coeffs,p,k,sign):
    modulus=p**k;squares={x*x%modulus for x in range(modulus)}
    params=[(s,1) for s in range(modulus)]+[(1,t) for t in range(0,modulus,p)]
    surviving=[]
    for s,t in params:
        values=[int(evaluate(h,s,t))*sign%modulus for h in coeffs]
        if all(x in squares for x in values):surviving.append([s,t])
    return {'prime':p,'exponent':k,'modulus':modulus,'sign':sign,'projective_classes_tested':len(params),
        'survivor_count':len(surviving),'first_survivor':surviving[0] if surviving else None,
        'conclusion':'LOCALLY_INSOLUBLE' if not surviving else 'FINITE_CONGRUENCE_SURVIVAL_ONLY'}


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,pari
    from sage.version import version
    source=r.read(old.OUTPUT);spec=r.read(PROTOCOL)
    R=PolynomialRing(QQ,['s','t']);s,t=R.gens();records=[];quartics=[]
    for row in source['records']:
        forms=[]
        for name in ['Q1','Q2']:
            M=matrix(QQ,3)
            for (i,j),c in zip(PAIRS,row[name]):
                M[i,j]=M[j,i]=QQ(c)/(1 if i==j else 2)
            forms.append(M)
        Q1,Q2=forms;assert all(x in ZZ for x in Q2.list())
        raw=pari.qfsolve(pari(Q2));assert raw.type()=='t_COL'
        base=vector(QQ,list(raw));assert base*Q2*base==0
        pivot=next(i for i in range(3) if base[i]);others=[i for i in range(3) if i!=pivot]
        v=vector(R,[0,0,0]);v[others[0]]=s;v[others[1]]=t
        z=(v*Q2*v)*base-2*(base*Q2*v)*v
        assert z*Q2*z==0
        h=R(z*Q1*z);coeff=[h.monomial_coefficient(s**(4-i)*t**i) for i in range(5)]
        assert all(x in ZZ for x in coeff)
        U=PolynomialRing(QQ,'x');x=U.gen();poly=U(h(x,1))
        assert poly.degree()==4 and poly.gcd(poly.derivative())==1
        quartics.append(poly)
        records.append({'beta_index':row['beta_index'],'Q1_Gram':[list(map(str,x)) for x in Q1.rows()],
            'Q2_Gram':[list(map(str,x)) for x in Q2.rows()],
            'qfsolve_point':list(map(str,base)),'pivot':pivot,'other_axes':others,
            'parametrization_coefficients':[[str(zi.monomial_coefficient(s**(2-j)*t**j)) for j in range(3)] for zi in z],
            'quartic_coefficients':list(map(str,coeff)),'quartic_discriminant':str(poly.discriminant())})
    f,g=quartics;res=f.resultant(g);assert res and f.gcd(g)==1
    coeffs=[x['quartic_coefficients'] for x in records]
    hits=[];tested=0
    for a,b in parameters(spec['limits']['rational_parameter_height']):
        tested+=1;values=[evaluate(h,a,b) for h in coeffs]
        if rational_square(values[0]*values[1]):
            root=r.F(isqrt((values[0]*values[1]).numerator),isqrt((values[0]*values[1]).denominator))
            assert root*root==values[0]*values[1]
            hits.append({'parameter':[a,b],'quartic_values':list(map(str,values)),
                'product_square_root':str(root),
                'signed_lifts':{str(e):all(rational_square(e*v) for v in values) for e in [1,-1]}})
    local=[local_test(coeffs,p,k,e) for p,k in spec['limits']['local_moduli'] for e in [1,-1]]
    return {'schema':'rank-jump.residual-double-covers.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,old.OUTPUT)},
        'software':{'sage':version,'pari':str(pari.version())},'status':'PASS','quartics':records,
        'resultant':str(res),'product_coefficients_ascending':list(map(str,f*g)),
        'product_quotient_genus':3,'simultaneous_cover_genus':5,'diagonal_quotient_degree':2,
        'diagonal_quotient_unramified':True,'rational_parameters_tested':tested,
        'product_quotient_hits':hits,'local_tests':local,
        'boundary':'No rational miss is an absence theorem. Local survival is not Q_p solubility. The product quotient is independent of sign, but its rational points need not lift.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'checkpoint.json'
    if not path.exists():
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--destination',str(path)],
                    cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                reason=None if p.returncode==0 else 'worker failure'
            except subprocess.TimeoutExpired:reason='30-second timeout'
            if reason and not path.exists():r.write_new(path,{'status':'UNKNOWN','reason':reason})
    data=r.read(path);r.write_new(OUTPUT,data)
    print(data.get('status'),len(data.get('product_quotient_hits',[])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);p.add_argument('--destination',type=Path)
    args=p.parse_args()
    if args.mode=='capture':capture()
    elif args.mode=='worker':r.write_new(args.destination,compute())
    else:assert r.read(OUTPUT)==compute();print('PASS residual double-cover replay')
