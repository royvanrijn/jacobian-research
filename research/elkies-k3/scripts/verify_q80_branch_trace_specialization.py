#!/usr/bin/env python3
"""Independent complete finite-group quotients for the branch trace gate."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-branch-trace-specialization-v1/result.json'


def val(q,p):
    if not q:
        return 10**9
    n,d=q.numerator,q.denominator
    v=0
    while n%p==0:n//=p;v+=1
    while d%p==0:d//=p;v-=1
    return v


def mod(q,p):return q.numerator*pow(q.denominator,-1,p)%p
def ev(c,t):return sum(F(x)*t**i for i,x in enumerate(c))


def rank(vectors):
    basis={}
    for v in vectors:
        while v:
            i=(v&-v).bit_length()-1
            if i not in basis:
                basis[i]=v;break
            v^=basis[i]
    return len(basis)


def add(P,Q,a,p):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    slope=((3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
    z=(slope*slope-x-u)%p
    return z,(slope*(x-z)-y)%p


def quotient(a,b,p):
    squares={}
    for y in range(p):squares.setdefault(y*y%p,[]).append(y)
    points=[None]+[(x,y) for x in range(p) for y in squares.get((x**3+a*x+b)%p,[])]
    group=set(points)
    doubled={add(P,P,a,p) for P in points}
    assert doubled<=group
    labels={P:0 for P in doubled};dimension=0
    while len(labels)<len(group):
        P=next(P for P in points if P not in labels)
        old=list(labels.items());bit=1<<dimension
        for Q,label in old:
            R=add(P,Q,a,p)
            assert R in group and R not in labels
            labels[R]=label^bit
        dimension+=1
    assert dimension<=2 and len(group)==len(doubled)*2**dimension
    return labels,dimension,len(group),len(doubled)


def verify(path):
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(20,20))
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    data=json.loads(path.read_text());source_path=ROOT/data['source']
    assert sha256(source_path.read_bytes()).hexdigest()==data['source_sha256']
    producer=ROOT/'elkies-k3/scripts/audit_q80_branch_trace_specialization.py'
    assert sha256(producer.read_bytes()).hexdigest()==data['script_sha256']
    source=json.loads(source_path.read_text());model=source['weierstrass_model'];sections=source['sections']['records']
    assert data['branch_parameter']==0 and data['rank']==17
    A=F(model['A_coefficients_low_to_high'][0]);B=F(model['B_coefficients_low_to_high'][0])
    coords=[]
    for section in sections:
        point=[]
        for k in ['X','Y']:
            f=section[k];point.append(F(f['numerator_coefficients_low_to_high'][0])/F(f['denominator_coefficients_low_to_high'][0]))
        coords.append(point)
    assert [str(A),str(B)]==[data['A'],data['B']]
    assert [[str(x),str(y)] for x,y in coords]==data['points']
    assert all(y*y==x**3+A*x+B for x,y in coords)
    selected_primes=sorted({data['rows'][r['row']]['prime'] for r in data['selected']})
    matrices=[];groups=[];balls=[]
    for p in selected_primes:
        assert p>=5 and all(p%d for d in range(2,int(p**0.5)+1))
        a,b=mod(A,p),mod(B,p);assert (4*a**3+27*b*b)%p
        labels,dim,order,double_order=quotient(a,b,p)
        point_labels=[labels[(mod(x,p),mod(y,p))] for x,y in coords]
        rows=[sum(((label>>i)&1)<<j for j,label in enumerate(point_labels)) for i in range(dim)]
        matrices+=rows
        oldrows=[r['code'] for r in data['rows'] if r['prime']==p]
        assert rank(rows)==rank(oldrows)==rank(rows+oldrows)
        groups.append({'prime':p,'group_order':order,'double_subgroup_order':double_order,
                       'quotient_dimension':dim,'point_labels':point_labels,'rows':rows})
        # Normalize each rational-function denominator to constant coefficient one.
        polynomials=[list(map(F,model[k])) for k in ['A_coefficients_low_to_high','B_coefficients_low_to_high']]
        for section in sections:
            for k in ['X','Y']:
                f=section[k];d0=F(f['denominator_coefficients_low_to_high'][0]);assert d0
                for part in ['numerator_coefficients_low_to_high','denominator_coefficients_low_to_high']:
                    polynomials.append([F(c)/d0 for c in f[part]])
        exponent=1
        for poly in polynomials:
            assert val(poly[0],p)>=0
            for i,c in enumerate(poly[1:],1):
                if c:exponent=max(exponent,(1-val(c,p)+i-1)//i)
        assert all(val(c,p)+i*exponent>=1 for poly in polynomials for i,c in enumerate(poly[1:],1))
        balls.append({'prime':p,'exponent':exponent})
    assert rank(matrices)==17
    # A fixed, prospective quartic realizes a compatible branch code at131.
    branches=[0,2,76,103];root_free=[]
    for t in branches:
        At,Bt=[ev(model[k],t) for k in ['A_coefficients_low_to_high','B_coefficients_low_to_high']]
        witness=None
        for p in data['prime_pool']:
            if At.denominator%p==0 or Bt.denominator%p==0:continue
            a,b=mod(At,p),mod(Bt,p)
            if (4*a**3+27*b*b)%p==0:continue
            if not any((x**3+a*x+b)%p==0 for x in range(p)):
                witness={'branch':t,'prime':p,'A':a,'B':b};break
        assert witness is not None
        root_free.append(witness)
    # Quartic C: w^2=t(t-2)(t-76)(t-103); inverse-coordinate cubic.
    k=-2*76*103;a2=2*76+2*103+76*103;a4=-(2+76+103)*k;a6=k*k
    slope=F(a4,2*k);X2=slope*slope-a2;Y2=-k-slope*X2
    assert Y2*Y2==X2**3+a2*X2**2+a4*X2+a6
    short_X2=9*X2+3*a2
    assert short_X2.denominator>1
    return {'status':'PASS','result_sha256':sha256(path.read_bytes()).hexdigest(),
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'independent_complete_finite_group_quotients':True,'rank':17,'groups':groups,
            'congruence_balls':balls,'integer_progression_modulus':str(prod(b['prime']**b['exponent'] for b in balls)),
            'control_branch_parameters':branches,'control_root_free_witnesses':root_free,
            'control_base_cubic':{'a2':a2,'a4':a4,'a6':a6,'point':[0,k]},
            'control_doubled_point':[str(X2),str(Y2)],'control_short_doubled_x':str(short_X2),
            'control_rank_gain':0,'positive_correlated_target_complete':False,
            'written_descent_formally_verified':False,'elapsed_seconds':time.monotonic()-started}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--record',type=Path);args=parser.parse_args();result=verify(args.input)
    if args.record:
        with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='groups'},sort_keys=True))
