#!/usr/bin/env python3
"""Replay one prospective Q80 genus-one pencil member and its exact rank gate."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import resource
import time
import sympy as S
from verify_q80_branch_trace_specialization import quotient, rank, val

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus1-single-trace-control-v1/input.json'

def verify(path):
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(20,20))
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    d=json.loads(path.read_text());source=ROOT/d['source']
    assert sha256(source.read_bytes()).hexdigest()==d['source_sha256']
    s=json.loads(source.read_text());t=S.symbols('t')
    def pol(cs):return S.Poly(sum(S.Rational(c)*t**i for i,c in enumerate(cs)),t,domain=S.QQ)
    model=s['weierstrass_model'];records=s['sections']['records']
    A=pol(model['A_coefficients_low_to_high']);B=pol(model['B_coefficients_low_to_high'])
    assert len(records)==17 and d['trace_index']==1 and d['parameter']==0
    r=records[1]
    X,D=[pol(r['X'][k+'_coefficients_low_to_high']) for k in ['numerator','denominator']]
    Y,DY=[pol(r['Y'][k+'_coefficients_low_to_high']) for k in ['numerator','denominator']]
    h=S.sqf_part(D).monic();Nx=(X*h*h).exquo(D);Ny=(Y*h**3).exquo(DY)
    M=(-Ny*S.invert(Nx,h*h))%(h*h)
    q=(M**4-6*M*M*Nx-8*M*Ny-3*Nx*Nx-4*A*h**4).exquo(h**6)
    for name,f in [('h',h),('Nx',Nx),('Ny',Ny),('M',M),('q',q)]:assert f==pol(d[name])
    assert h.degree()==2 and q.degree()==4
    assert S.gcd(q,q.diff()).degree()==0 and S.gcd(q,4*A**3+27*B**2).degree()==0
    x0=(M*M-Nx).exquo(2*h*h);y0=(M*(M*M-3*Nx)-2*Ny).exquo(2*h**3)
    x1=h.mul_ground(S.Rational(1,2));y1=M.mul_ground(S.Rational(1,2))
    assert (y0*y0+y1*y1*q-x0**3-3*x0*x1*x1*q-A*x0-B).is_zero
    assert (2*y0*y1-3*x0*x0*x1-x1**3*q-A*x1).is_zero
    assert not x1.is_zero
    # At the branch the residual chord is tangent: twice its point equals T.
    inv=lambda f:S.invert(f,q)
    xb=(x0%q);yb=(y0%q);slope=((3*xb*xb+A)*inv(2*yb))%q
    xx=(slope*slope-2*xb)%q;yy=(slope*(xb-xx)-yb)%q
    assert ((xx*D-X)%q).is_zero and ((yy*DY-Y)%q).is_zero
    p=d['irreducibility_prime'];assert S.isprime(p)
    def reduced(f,p):
        cs=[int(c.p)*pow(int(c.q),-1,p)%p for c in f.all_coeffs()]
        assert cs[0]
        return S.Poly.from_list(cs,t,modulus=p).monic()
    qp=reduced(q,p);z=S.Poly(t,t,modulus=p)
    def power(f,n):
        out=S.Poly(1,t,modulus=p)
        while n:
            if n&1:out=(out*f)%qp
            f=(f*f)%qp;n//=2
        return out
    assert (power(z,p**4)-z).rem(qp).is_zero
    assert S.gcd(qp,power(z,p**2)-z).degree()==0
    def ev(cs,v,p):
        out=0
        for c in reversed(cs):
            c=F(c);out=(out*v+c.numerator*pow(c.denominator,-1,p))%p
        return out
    groups=[];rows=[]
    for place in d['places']+[d['no_2_torsion_place']]:
        p,v=place['prime'],place['t'];assert S.isprime(p) and p>3
        assert ev([d['q'][-1]],0,p) and ev(d['q'],v,p)==0
        assert ev([str((i+1)*F(d['q'][i+1])) for i in range(4)],v,p)
        a,b=[ev(model[k],v,p) for k in ['A_coefficients_low_to_high','B_coefficients_low_to_high']]
        assert (4*a**3+27*b*b)%p
        points=[]
        for r in records:
            point=tuple(ev(r[k]['numerator_coefficients_low_to_high'],v,p)*pow(ev(r[k]['denominator_coefficients_low_to_high'],v,p),-1,p)%p for k in ['X','Y'])
            x,y=point;assert (y*y-x**3-a*x-b)%p==0;points.append(point)
        labels,dim,order,double_order=quotient(a,b,p)
        local=[sum(((labels[P]>>i)&1)<<j for j,P in enumerate(points)) for i in range(dim)]
        assert all(row&2==0 for row in local)
        rows+=local;groups.append({'prime':p,'t':v,'order':order,'double_order':double_order,'rows':local})
        if place==d['no_2_torsion_place']:assert order%2 and order==place['order']
    assert rank(rows)==16
    # Freeze an arithmetic progression, not a finite sample of pencil members.
    qj=[q,(4*M**3-12*M*Nx-8*Ny).exquo(h**4),
        (6*(M*M-Nx)).exquo(h*h),4*M,h*h]
    assert all(f.degree()<=4 for f in qj)
    balls=[]
    for p in sorted({d['irreducibility_prime']}|{g['prime'] for g in groups}):
        assert all(val(F(str(c)),p)>=0 for f in qj for c in f.all_coeffs())
        exponent=1
        for j,f in enumerate(qj[1:],1):
            for c in f.all_coeffs():
                if c:exponent=max(exponent,(1-val(F(str(c)),p)+j-1)//j)
        assert all(val(F(str(c)),p)+j*exponent>=1 for j,f in enumerate(qj[1:],1) for c in f.all_coeffs())
        balls.append({'prime':p,'exponent':exponent})
    return {'status':'PASS','input_sha256':sha256(path.read_bytes()).hexdigest(),
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'finite_quotient_helper_sha256':sha256((Path(__file__).parent/'verify_q80_branch_trace_specialization.py').read_bytes()).hexdigest(),
            'branch_degree':4,'cover_genus':1,'trace_index':1,'parameter':0,
            'irreducibility_prime':d['irreducibility_prime'],'groups':groups,
            'matrix_rank':16,'branch_halving_kernel_basis':[2],
            'branch_field_2_torsion':False,'exact_lift_identity':True,'exact_branch_half':True,
            'section_coefficient_degrees':[x0.degree(),x1.degree(),y0.degree(),y1.degree()],
            'congruence_balls':balls,'parameter_progression_modulus':str(prod(b['prime']**b['exponent'] for b in balls)),
            'rank_gain':1,'function_field_rank':18,'base_rational_points_certified':False,
            'positive_correlated_target_complete':False,'formal_verification':False,
            'elapsed_seconds':time.monotonic()-started}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT);parser.add_argument('--record',type=Path)
    args=parser.parse_args();result=verify(args.input)
    if args.record:
        with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True))
