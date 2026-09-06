#!/usr/bin/env python3
"""Exact common-value rational block: two fixed small families, no point search."""
import argparse
from math import gcd
from pathlib import Path
import subprocess
import retrospective as r

PROTOCOL=Path(__file__).with_name('SHARED_VALUE_SOLUBLE_BLOCK_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_shared_value_soluble_block_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-shared-value-soluble-block-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL)}


def count_points(coeffs,p):
    return 1+sum(1 if (v:=sum(int(a)*pow(x,i,p) for i,a in enumerate(coeffs))%p)==0
                 else 2 if pow(v,(p-1)//2,p)==1 else 0 for x in range(p))


def compute(c):
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,EllipticCurve,QuadraticField
    from sage.version import version
    spec=r.read(PROTOCOL);c=QQ(c);D=c*c-2
    R=PolynomialRing(QQ,'m');m=R.gen();K=R.fraction_field()
    X=PolynomialRing(K,'x');x=X.gen()
    f=x**3+m*x*x-(m+3)*x+c*c
    q=m*m+3*m+9;disc=R(f.discriminant());Delta=16*disc;c4=16*q
    assert f(-1)==f(2)==2*m+2+c*c
    assert f(0)==c*c and f(1)==D
    assert f(c*c/2)==c*c*D*(2*m+2+c*c)/8
    assert Delta.degree()==4 and c4.degree()==2
    if c==1:
        assert disc==q*q
        finite={'type':'II','count':2,'discriminant_order':2,'c4_order':1,'root_lattice_rank':0}
    else:
        assert disc.gcd(disc.derivative())==1 and disc.gcd(c4)==1
        finite={'type':'I1','count':4,'discriminant_order':1,'c4_order':0,'root_lattice_rank':0}
    # Infinity model x=v^-2 X, y=v^-3 Y: v(c4)=2, v(Delta)=8, type I2*.
    parent_geo=10-(2+6)
    N=PolynomialRing(QQ,'n');n=N.gen();mn=(n*n-c*c-2)/2
    an=[N(0),mn,N(0),-mn-3,N(c*c)]
    En=EllipticCurve(N.fraction_field(),an)
    delta_n=N(En.discriminant());c4_n=N(En.c4())
    assert delta_n==Delta(mn) and c4_n==c4(mn)
    assert delta_n.degree()==8 and c4_n.degree()==4
    # The branch point d=0 is smooth for these fixed c; no omitted singularity.
    assert disc(-(c*c+2)/2)!=0
    if c==1:
        assert delta_n==(n**4+27)**2
        base_finite={'type':'II','count':4,'root_lattice_rank':0}
    else:
        assert delta_n.gcd(delta_n.derivative())==1 and delta_n.gcd(c4_n)==1
        base_finite={'type':'I1','count':8,'root_lattice_rank':0}
    # Infinity now has v(c4)=0, v(Delta)=4, hence I4.
    base_geo=10-(2+3)
    points=[(0,c),(-1,n),(2,n)]
    for a,b in points:assert b*b==a**3+mn*a*a-(mn+3)*a+c*c
    assert D==1+mn-(mn+3)+c*c
    assert c*c*n*n*(2*D)/16==(c*c/2)**3+mn*(c*c/2)**2-(mn+3)*(c*c/2)+c*c

    # A single frozen specialization certifies the rational three-section independence.
    n0=QQ(spec['specializations']['independence_anchor_n']);m0=QQ(mn(n0))
    coefficients=[c*c,-m0-3,m0,QQ(1)]
    E=EllipticCurve(QQ,[0,m0,0,-m0-3,c*c]);ps=[E(0,c),E(-1,n0),E(2,n0)]
    P=PolynomialRing(QQ,'z');pol=P(coefficients)
    assert P.change_ring(GF(2))(pol).is_irreducible() # no rational two-torsion
    signatures=[0,0,0];offset=0;fingerprints=[]
    realroots=pol.roots(AA,multiplicities=False)
    for j,point in enumerate(ps):signatures[j]|=r.pack(int(point[0]<a) for a in realroots)
    offset=len(realroots)
    for p in r.primes(spec['limits']['largest_fingerprint_prime']):
        if pol.discriminant()%p==0:continue
        fp=PolynomialRing(GF(p),'z')(pol);roots=sorted(int(t) for t in fp.roots(multiplicities=False))
        if len(roots)!=3:continue
        cols=[]
        for j,point in enumerate(ps):
            bits=[]
            for root in roots:
                value=(int(point[0])-root)%p
                if value==0:value=int(fp.derivative()(root))
                bits.append(int(pow(value,(p-1)//2,p)==p-1))
            sig=r.pack(bits);cols.append(sig);signatures[j]|=sig<<offset
        fingerprints.append({'prime':p,'roots':roots,'point_signatures':cols});offset+=3
    assert r.rank(signatures)==3

    # Certify the two remaining geometric directions are nontorsion over
    # their respective constant quadratic fields, using finite reduction bounds.
    nonrational=[]
    for rad,px,ymult,label in [(D,QQ(1),QQ(1),'inherited'),(2*D,c*c/2,c*n0/4,'new')]:
        F=QuadraticField(rad,'a');EF=E.change_ring(F);point=EF(px,ymult*F.gen())
        pool=[]
        for p in r.primes(spec['limits']['largest_torsion_bound_prime']):
            if pol.discriminant()%p==0 or ZZ(rad)%p==0 or pow(int(rad)%p,(p-1)//2,p)!=1:continue
            pool.append((p,count_points(coefficients,p)))
        pair=next(((p,np,q,nq) for p,np in pool for q,nq in pool if p<q and np%q and nq%p),None)
        assert pair is not None
        p,np,q,nq=pair;bound=gcd(np,nq);multiple=bound*point
        assert not multiple.is_zero()
        nonrational.append({'role':label,'quadratic_radicand':str(rad),
            'point_x':str(px),'point_y_multiplier':str(ymult),
            'split_good_prime_counts':[[p,np],[q,nq]],'torsion_order_divides':bound,
            'bound_multiple_is_nonzero':True,'bound_multiple_coordinates':list(map(str,multiple))})
    # The geometric basis has three rational vectors and distinct constant characters D, 2D.
    assert not D.is_square() and not (2*D).is_square() and not QQ(2).is_square()
    galois_diagonals=[[1,1,1,(-1)**a,(-1)**b] for a in range(2) for b in range(2)]
    collapse=None
    if c==1:
        E1=EllipticCurve(QQ,[0,-1,0,-2,1]);R1=E1(0,1)
        assert 2*R1==E1(2,1) and -3*R1==E1(-1,1)
        counts=[[p,count_points([1,-2,-1,1],p)] for p in [11,13]]
        assert counts==[[11,15],[13,16]]
        assert gcd(15,16)==1 and 15%13 and 16%11
        collapse={'n':1,'m':-1,'relations':['P_2 = 2 R','P_-1 = -3 R'],
                  'retained_specialized_subgroup_rank':1,'good_prime_counts':counts,
                  'whole_curve_rank':'UNKNOWN here; only the retained subgroup is certified'}
    return {'bindings':bindings(),'c':int(c),'software':{'sage':version},
        'parent_discriminant_coefficients':list(map(str,Delta)),
        'parent_c4_coefficients':list(map(str,c4)),'parent_finite_fibres':finite,
        'parent_infinity':{'type':'I2*','c4_order':2,'discriminant_order':8,'root_lattice_rank':6},
        'parent_geometric_generic_rank':parent_geo,'parent_arithmetic_generic_rank':1,
        'base_discriminant_coefficients':list(map(str,delta_n)),'base_c4_coefficients':list(map(str,c4_n)),
        'base_finite_fibres':base_finite,'base_infinity':{'type':'I4','c4_order':0,'discriminant_order':4,'root_lattice_rank':3},
        'base_geometric_generic_rank':base_geo,'base_arithmetic_generic_rank':3,'new_arithmetic_generic_rank':2,
        'constant_galois_diagonals':galois_diagonals,'constant_character_radicands':[str(D),str(2*D)],
        'anchor':{'n':int(n0),'m':str(m0),'rational_points':[list(map(str,point)) for point in ps],
            'cubic_discriminant':str(pol.discriminant()),'rational_two_torsion':False,
            'joint_fingerprint_masks':signatures,'independent_point_rank':3,'finite_fingerprints':fingerprints},
        'nonrational_nontorsion_witnesses':nonrational,'collapse_control':collapse,
        'claim_boundary':'Exact generic ranks for two fixed small families. No assertion that every specialization retains the block, no rank result for a production curve.'}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for c in r.read(PROTOCOL)['fixed_c_values']:
        if check:
            assert compute(c)==next(x for x in r.read(OUTPUT)['rows'] if x['c']==c)
            print('PASS symbolic and arithmetic replay',c,flush=True);continue
        path=WORK/f'c-{c}.json'
        if not path.exists():
            with (WORK/f'c-{c}.log').open('x') as log:
                proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--c',str(c),
                    '--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                if proc.returncode:raise RuntimeError(f'failed worker; see {path.with_suffix(".log")}')
        row=r.read(path);assert row['bindings']==bindings();rows.append(row);print('checkpoint',c,'PASS',flush=True)
    if not check:r.write_new(OUTPUT,{'schema':'rank-jump.shared-value-soluble-block.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check'])
    p.add_argument('--c',type=int);p.add_argument('--destination',type=Path);args=p.parse_args()
    if args.mode=='worker':r.write_new(args.destination,compute(args.c))
    else:capture(args.mode=='check')
