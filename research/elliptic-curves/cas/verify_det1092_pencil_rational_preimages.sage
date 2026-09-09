#!/usr/bin/env sage-python
"""Factorization-free rational-root certificates for five frozen cover fibres.

At one prime from the unchanged12-prime pool, enumerate all simple roots,
lift to modulus M>8*max(|f(0)|,|lc(f)|)^2, and exactly Gauss-reduce each
congruence lattice. This proves rational-root absence with elementary
integer identities, independently of the Q-factorization constructor.
One --case per process, <=25 seconds. No elliptic point search.
"""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_pencil_exact_incidence_v1';OLD=ART/'det1092_pencil_multiples_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,obj):
    path=OUT/name
    if path.exists():assert json.loads(path.read_text())==obj
    else:
        with path.open('x') as stream:json.dump(obj,stream,indent=2,sort_keys=True);stream.write('\n')
def evaluate(coeff,x,m):
    value=0
    for c in reversed(coeff):value=(value*x+int(c))%m
    return value
def gauss(a,b):
    a=vector(ZZ,a);b=vector(ZZ,b)
    while True:
        if b*b<a*a:a,b=b,a
        q=(2*(a*b)+a*a)//(2*(a*a))
        if q==0:return a,b
        b-=q*a
def check_reduced(a,b,M,r):
    assert abs(matrix(ZZ,[a,b]).det())==M
    assert all((x[0]-r*x[1])%M==0 for x in [a,b])
    assert a*a<=b*b and 2*abs(a*b)<=a*a
def root_certificate(f,pool):
    assert f[0]!=0 and f.leading_coefficient()!=0 and all(c in ZZ for c in f.list())
    coeff=list(map(int,f.list()));der=[i*c for i,c in enumerate(coeff)][1:]
    A=abs(ZZ(f[0]));B=abs(ZZ(f.leading_coefficient()));H=max(A,B)
    chosen=None;attempts=[]
    for p in pool:
        if coeff[-1]%p==0:attempts.append({'prime':p,'gate':'leading_coefficient_zero'});continue
        roots=[x for x in range(p) if evaluate(coeff,x,p)==0]
        simple=all(evaluate(der,x,p)!=0 for x in roots)
        attempts.append({'prime':p,'gate':'PASS' if simple else 'multiple_root','roots':roots})
        if simple:chosen=(p,roots);break
    assert chosen is not None,'No eligible prime in the unchanged pool; retain UNKNOWN.'
    p,roots=chosen;reports=[];rational=[]
    for root in roots:
        r=int(root);M=int(p)
        while M<=8*H*H:
            nextM=M*M
            r=(r-evaluate(coeff,r,nextM)*pow(evaluate(der,r,nextM),-1,nextM))%nextM
            M=nextM
        assert r%p==root and evaluate(coeff,r,M)==0
        a,b=gauss([M,0],[r,1]);check_reduced(a,b,M,r)
        record={'initial_root':root,'modulus':str(M),'lift':str(r),
                'reduced_basis':[[str(c) for c in v] for v in [a,b]]}
        if a*a>2*H*H:record['decision']='NO_BOUNDED_NUMERATOR_DENOMINATOR_VECTOR'
        else:
            assert a[1]!=0
            q=QQ(a[0])/a[1];is_root=bool(f(q)==0)
            record.update(decision='RATIONAL_ROOT' if is_root else 'ONLY_POSSIBLE_RATIO_IS_NOT_A_ROOT',candidate=str(q))
            if is_root:rational.append(str(q))
        reports.append(record)
    return {'prime':p,'prime_attempts':attempts,'all_roots_mod_p':roots,
            'numerator_bound':str(A),'denominator_bound':str(B),'common_bound':str(H),
            'lifts':reports,'rational_roots':rational}

parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,required=True);args=parser.parse_args()
assert 0<=args.case<5
oldprotocol=json.loads((OLD/'protocol.json').read_text());pool=oldprotocol['ramification_certificate_primes']
save('root-replay-protocol.json',{'classification':'independent elementary rational-root proof',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OLD/'protocol.json',OUT/'protocol.json',Path(__file__)]},
    'limits':{'seconds_per_case':25,'cases':5,'prime_pool':pool,'factorizations':0,
              'new_parameters':0,'point_searches':0},
    'method':'Rational-root theorem bounds; simple Hensel lifts; exact2D Gauss certificate at M>8H^2. First eligible prime from the old pool.'})
R=PolynomialRing(QQ,'z');z=R.gen()
# Complete positive and negative regressions, without using factorization.
tests=[(3*z-2,['2/3']),(5*z+7,['-7/5']),((3*z-2)*(5*z+7),['-7/5','2/3']),
       (z*z-2,[]),(z**3-z+1,[])]
for f,expected in tests:assert sorted(root_certificate(f,pool)['rational_roots'])==sorted(expected)
inp=json.loads((OUT/('input-%02d.json'%args.case)).read_text());f=R(inp['primitive_incidence_polynomial'])
case=inp['case'];row=json.loads((OLD/('multiple%d.json'%case['n'])).read_text())
N=R(row['t_of_z']['numerator']);D=R(row['t_of_z']['denominator']);g=N-QQ(case['parameter'])*D
g*=lcm([q.denominator() for q in g.list()]);g/=gcd([ZZ(q) for q in g.list()])
assert f==g or f==-g
assert f.degree()==inp['map_degree'] and not inp['infinity_preimage'] and N.gcd(D).degree()==0
proof=root_certificate(f,pool);assert proof['rational_roots']==[]
report={'status':'PASS_NO_RATIONAL_PROJECTIVE_PREIMAGE','classification':'verified application with independent elementary proof',
        'case':case,'degree':int(f.degree()),'root_certificate':proof,
        'synthetic_regressions':5,'checker_sha256':sha(Path(__file__)),
        'boundary':'Excludes this frozen cover only, not other covers, rational elliptic points, or rank jumps.'}
save('root-replay-%02d.json'%args.case,report)
print(args.case,case['label'],'degree',f.degree(),'prime',proof['prime'],
      'lifted roots',len(proof['lifts']),report['status'],flush=True)
