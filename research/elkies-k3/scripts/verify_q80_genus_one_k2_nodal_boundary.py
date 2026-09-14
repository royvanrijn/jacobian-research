#!/usr/bin/env python3
"""Independent finite replay for the Q80 genus1 k2 nodal boundary theorem."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path
import resource
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
NORM4='artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
PRODUCER='elkies-k3/scripts/certify_q80_genus_one_k2_nodal_boundary.sage'
CPP='elkies-k3/scripts/certify_q80_nodal_quadratic_census.cpp'
OLD='elkies-k3/scripts/verify_q80_genus_zero_closure.py'
P=131
def require(v,message):
    if not v:raise ValueError(message)
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def trim(a):
    a=list(a)
    while a and a[-1]==0:a.pop()
    return a
def add(a,b,m=P):return trim([((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%m for i in range(max(len(a),len(b)))])
def scale(a,c,m=P):return trim([v*c%m for v in a])
def sub(a,b,m=P):return add(a,scale(b,-1,m),m)
def mul(a,b,m=P):
    if not a or not b:return []
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim([x%m for x in c])
def divrem(a,b,m=P):
    a=trim(a);b=trim(b);require(bool(b),'division by zero polynomial')
    quotient=[0]*max(0,len(a)-len(b)+1);inverse=pow(b[-1],-1,m)
    while a and len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]*inverse%m;quotient[k]=c
        for j,v in enumerate(b):a[k+j]=(a[k+j]-c*v)%m
        a=trim(a)
    return trim(quotient),a
def exact(a,b):
    q,r=divrem(a,b);require(not r,'exact polynomial division');return q
def ev(a,b):
    value=0
    for x in reversed(a):value=(value*b+x)%P
    return value
def gcd(a,b):
    while b:a,b=b,divrem(a,b)[1]
    return scale(a,pow(a[-1],-1,P)) if a else []
def powmod(a,n,f):
    b=[1]
    while n:
        if n&1:b=divrem(mul(b,a),f)[1]
        a=divrem(mul(a,a),f)[1];n>>=1
    return b
def irreducible(f):
    n=len(f)-1;require(n in (2,3,19),'prime-degree factor certificate')
    h=[0,1]
    for k in range(1,n+1):
        h=powmod(h,P,f)
        if k==1 and len(gcd(sub(h,[0,1]),f))!=1:return False
    return h==[0,1]
def square_up_to_scalar(f):
    f=trim(f)
    if not f:return True
    offset=next(i for i,v in enumerate(f) if v)
    if offset%2 or (len(f)-1)%2:return False
    f=scale(f[offset:],pow(f[offset],-1,P));n=(len(f)-1)//2
    y=[1]+[0]*n
    for k in range(1,n+1):y[k]=(f[k]-sum(y[j]*y[k-j] for j in range(1,k)))*66%P
    return mul(y,y)==f
def rank(rows):
    a=[[v%P for v in row] for row in rows];r=0
    for c in range(len(a[0]) if a else 0):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];inverse=pow(a[r][c],-1,P)
        for i in range(r+1,len(a)):
            factor=a[i][c]*inverse%P
            for j in range(c,len(a[i])):a[i][j]=(a[i][j]-factor*a[r][j])%P
        r+=1
        if r==len(a):break
    return r
def fraction_poly(values,m):
    return trim([Fraction(v).numerator*pow(Fraction(v).denominator,-1,m)%m for v in values])
def nodal_census(A,B,root):
    """Check all 131^3 x values; use constant-term rather than leading-term roots."""
    z0=np.tile(np.arange(P,dtype=np.int64),P);z1=np.repeat(np.arange(P,dtype=np.int64),P)
    inverse=np.array([0]+[pow(i,-1,P) for i in range(1,P)],dtype=np.int64);found=[]
    for z2 in range(P):
        x=np.zeros((5,P*P),dtype=np.int64)
        x[0]=(root[0]+88*z0)%P;x[1]=(root[1]+62*z0+88*z1)%P
        x[2]=(z0+62*z1+88*z2)%P;x[3]=(z1+62*z2)%P;x[4]=z2
        square=np.zeros((9,P*P),dtype=np.int64)
        for i in range(5):
            for j in range(5):square[i+j]+=x[i]*x[j]
        square%=P;f=np.zeros((13,P*P),dtype=np.int64)
        for i in range(9):
            for j in range(5):f[i+j]+=(square[i]+A[i])*x[j]
        for i in range(13):f[i]+=B[i]
        f%=P;s=np.zeros((11,P*P),dtype=np.int64)
        for i in range(12,1,-1):
            s[i-2]=f[i];f[i-1]=(f[i-1]-62*f[i])%P;f[i-2]=(f[i-2]-88*f[i])%P
        require(not np.any(f[:2]),'nodal-root divisibility of every coefficient tuple')
        require(not np.any(s[0]==0),'nonzero constant term from the root-free fibre at zero')
        normalized=s*inverse[s[0]]%P;y=np.zeros((6,P*P),dtype=np.int64);y[0]=1
        for i in range(1,6):
            coefficient=np.zeros(P*P,dtype=np.int64)
            for j in range(1,i):coefficient+=y[j]*y[i-j]
            y[i]=(normalized[i]-coefficient)*66%P
        good=np.ones(P*P,dtype=bool)
        for i in range(6,11):
            coefficient=np.zeros(P*P,dtype=np.int64)
            for j in range(max(0,i-5),6):coefficient+=y[j]*y[i-j]
            good &= coefficient%P==normalized[i]
        for k in np.flatnonzero(good):found.append(([int(z0[k]),int(z1[k]),z2],[int(v) for v in x[:,k]]))
    return found
def verify(out):
    started=time.monotonic();packet=read(out/'input.json');result=read(out/'result.json')
    expected={SOURCE,PRODUCER,CPP,OLD,str(Path(__file__).relative_to(ROOT)),
              *[DEGREE+'/'+n for n in ['input.json','result.json','rational-fibres.json','quadratic-fibres.json','independent-replay.json']],
              *[NORM4+'/'+n for n in ['input.json','result.json','norm4-sections.json','complete-replay-input.json','complete-independent-replay.json']]}
    require(set(packet['bindings'])==expected,'exact source bindings')
    for s,h in {**packet['bindings'],**packet['preserved_preflight']}.items():require(digest(ROOT/s)==h,'source or preflight binding '+s)
    require(packet['prime']==P and packet['cpu_seconds']==40 and packet['memory_bytes']==4*1024**3,'frozen resource and arithmetic scope')
    require((packet['nodal_census_count'],packet['etale_coefficient_count'],packet['split_contact_pairs'])==(P**3,P,21),'complete finite ranges')
    require(result['status']=='PASS' and result['input_sha256']==digest(out/'input.json'),'result binding')
    require(set(result['records'])=={'character-gate.json','etale-coefficients.json','split-lifts.json','nodal-census.json'},'complete result records')
    for s,h in result['records'].items():require(digest(out/s)==h,'result record '+s)
    spec=importlib.util.spec_from_file_location('q80_k2_prior',ROOT/OLD);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    old_replay=old.verify(ROOT/NORM4)
    require(old_replay['finite_replay']['norm4_sections']==1313,'full reduced norm4 lattice and points')
    source=read(ROOT/SOURCE)['weierstrass_model']
    A=fraction_poly(source['A_coefficients_low_to_high'],P);B=fraction_poly(source['B_coefficients_low_to_high'],P)
    A2=fraction_poly(source['A_coefficients_low_to_high'],P*P);B2=fraction_poly(source['B_coefficients_low_to_high'],P*P)
    gate=read(out/'character-gate.json');delta=add(scale(mul(mul(A,A),A),4),scale(mul(B,B),27))
    product=[1]
    require([len(h)-1 for h,e in gate['discriminant_factors']]==[2,3,19],'all singular-fibre residue degrees')
    for h,e in gate['discriminant_factors']:
        require(e==1 and irreducible(h),'squarefree irreducible discriminant factor');product=mul(product,h)
    require(product==scale(delta,pow(delta[-1],-1,P)),'complete degree24 discriminant')
    q=[88,62,1];node=gate['node_root'];simple=gate['simple_root']
    require(gate['discriminant_factors'][0][0]==q,'unique quadratic nodal orbit')
    for x in [node,simple]:require(not divrem(add(add(mul(mul(x,x),x),mul(A,x)),B),q)[1],'selected nodal-fibre root')
    require(not divrem(add(scale(mul(node,node),3),A),q)[1],'double root')
    require(bool(divrem(add(scale(mul(simple,simple),3),A),q)[1]),'simple root')
    require(ev(node,100)==122 and ev(simple,100)==18,'encoded extension root real parts')
    require(node[1]*31%P==38 and simple[1]*31%P==55,'encoded extension root imaginary parts')
    rational=read(ROOT/DEGREE/'rational-fibres.json')['rows'];quadratic=read(ROOT/DEGREE/'quadratic-fibres.json')['rows']
    rp=[(r['t'],e,c) for r in rational for e,c in zip(r['rational_roots'],r['rational_codes'])]
    require(len(rp)==len({r[2] for r in rp})==110 and all(r[2] for r in rp),'all110 character codes nonzero and distinct')
    zeros=[r for r in quadratic if 0 in r['norm_codes']]
    require(zeros==gate['quadratic_zero_rows'] and len(zeros)==1 and zeros[0]['t']==4161,'unique degree2 zero norm orbit')
    pairs=[[list(a[:2]),list(b[:2])] for a,b in itertools.combinations_with_replacement(rp,2) if a[2]^b[2]==13412]
    norm=[{'t':r['t'],'root':e,'smooth':r['smooth']} for r in rational+quadratic for e,c in zip(r['roots'],r['norm_codes']) if c==13412]
    require(pairs==gate['rational_pairs_norm13412']==[[[35,114],[75,11]]],'complete rational norm13412 pairs')
    require(norm==gate['quadratic_norm13412']==[{'t':4161,'root':5100,'smooth':False}]*2,'complete quadratic norm13412 support')
    split=[r['t'] for r in rational if len(r['rational_roots'])==3];require(split==gate['split_sites'] and len(split)==16,'full split fibre roster')
    etale=read(out/'etale-coefficients.json');D=mul(q,mul([-35%P,1],[-75%P,1]));x0=etale['x0']
    require(D==etale['D'] and len(x0)<=4 and divrem(x0,q)[1]==node and ev(x0,35)==114 and ev(x0,75)==11,'unique degree-below4 interpolation')
    require([r['c'] for r in etale['rows']]==list(range(P)),'all131 coefficient choices')
    for row in etale['rows']:
        x=add(x0,scale(D,row['c']));s=exact(add(add(mul(mul(x,x),x),mul(A,x)),B),D)
        require(s==row['quotient'] and not square_up_to_scalar(s) and row['square_up_to_scalar'] is False,'etale coefficient exclusion')
    points=read(ROOT/NORM4/'norm4-sections.json')['records'];incidence=defaultdict(list)
    for point in points:
        for b in list(range(P))+[None]:
            x=(point['x'][4] if len(point['x'])>4 else 0) if b is None else ev(point['x'],b)
            y=(point['y'][6] if len(point['y'])>6 else 0) if b is None else ev(point['y'],b)
            if not y:incidence[b].append((point['index'],x))
    by_pair=defaultdict(dict)
    for b,rows in incidence.items():
        for (i,e),(j,f) in itertools.combinations(rows,2):by_pair[(i,j)][b]=(e==f)
    complete={(ij,b,c) for ij,sites in by_pair.items() for b,equal in sites.items() if not equal for c,agree in sites.items() if agree}
    rows=read(out/'split-lifts.json')['rows'];require(len(rows)==len(complete)==21,'complete mixed-contact pair roster')
    require({(tuple(r['pair']),r['disagree'],r['agree']) for r in rows}==complete,'no omitted mixed-contact pair')
    for row in rows:
        b,c=row['disagree'],row['agree'];require(b in split and b!=c,'smooth distinct contact sites')
        H=mul([-b%P,1],[1] if c is None else [-c%P,1]);D=mul(H,H);fixed=2 if c is None else 4
        require(row['fixed_D_coefficient']==fixed and D[fixed]==1,'literal scalar coefficient chart including infinity')
        J=[[0]*24 for _ in range(26)];target=[]
        for n,index in enumerate(row['pair']):
            x=points[index]['x'];s=exact(points[index]['y'],H);require(len(s)<=5,'degree4 node factor')
            local=[[0]*k+add(scale(mul(x,x),3),A) for k in range(5)]
            local += [[0]*k+scale(mul(D,s),-2) for k in range(5)]
            shared=[[0]*k+scale(mul(s,s),-1) for k in range(5) if k!=fixed]
            for j,f in enumerate(local):
                for k,v in enumerate(f):J[13*n+k][10*n+j]=v
            for j,f in enumerate(shared):
                for k,v in enumerate(f):J[13*n+k][20+j]=v
            error=sub(add(add(mul(mul(x,x,P*P),x,P*P),mul(A2,x,P*P),P*P),B2,P*P),mul(D,mul(s,s,P*P),P*P),P*P)
            error += [0]*(13-len(error));require(all(v%P==0 for v in error),'first-lift error divisible by131')
            target += [(-v//P)%P for v in error]
        require(target==row['error_target'],'literal rational coefficients modulo131 squared')
        require(rank(J)==row['jacobian_rank']==24 and rank([r+[v] for r,v in zip(J,target)])==row['augmented_rank']==25,'inconsistent first lift')
    census=read(out/'nodal-census.json');found=nodal_census(A,B,simple)
    require(census['status']=='PASS' and census['tried']==P**3 and census['accepted']==len(found)==2,'independent complete131-cubed census')
    require(found==[(r['z'],r['x']) for r in census['records']],'independent constant-term and leading-term censuses agree')
    for row in census['records']:
        x,s=row['x'],row['s']
        require(divrem(x,q)[1]==simple and add(add(mul(mul(x,x),x),mul(A,x)),B)==scale(mul(q,mul(s,s)),row['scalar']),'exact residual section equation')
        roots=[b for b in range(P) if not ev(s,b)]
        require(roots==row['rational_s_roots'] and not set(roots)&set(split),'no required full-splitting contact')
    square=mul(q,q);I=(12*square[4]*square[0]-3*square[3]*square[1]+square[2]**2)%P
    require(I==60 and pow((62*62-4*88)%P,65,P)==P-1,'remaining quartic invariant and anisotropic branch quadratic')
    require(result['necessary_D_reduction']==square and result['both_abscissas_require_negative_gauss_valuation'] is True,'written theorem boundary')
    require(result['remaining_nodal_denominator_boundary']=='UNKNOWN' and result['positive_mw17_target_complete'] is False and result['old_norm8_replay_upgraded'] is False,'open endpoint and assurance')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
            'independent_finite_replay':True,'discriminant_factor_degrees':[2,3,19],
            'etale_coefficient_exclusions':P,'split_pair_lift_exclusions':21,
            'nodal_coefficient_tuples':P**3,'nodal_points':2,'valid_nodal_split_contacts':0,
            'remaining_D_mod131':square,'remaining_quartic_I_mod131':I,
            'written_stable_map_and_branch_classification_formally_verified':False,
            'remaining_boundary':'UNKNOWN','positive_mw17_target_complete':False,
            'old_norm8_replay_upgraded':False,'elapsed_seconds':round(time.monotonic()-started,6)}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT);parser.add_argument('--record',type=Path);args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    result=verify(args.input)
    if args.record:
        with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)
