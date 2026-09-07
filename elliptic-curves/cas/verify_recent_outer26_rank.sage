#!/usr/bin/env python3
"""Standalone independent finite-group ranks on the two retained outer26 clouds."""
import json,hashlib,argparse
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,Matrix,lcm,gcd,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUTS=[ART/('recent_outer26_specialized_'+n+'_mod2_v1.json') for n in ('189','192')]
OUT=ART/'recent_outer26_sage_replay_v1.json'

def check_one(INPUT):
    d=json.loads(INPUT.read_bytes())
    E=EllipticCurve(QQ,d['curve'])
    points=[E([QQ(x),QQ(y)]) for x,y in d['points']]
    if len(points)<26:raise ArithmeticError('complete point cloud required')
    projective=[]
    for P in points:
        den=lcm([c.denominator() for c in P]);v=[ZZ(c*den) for c in P];g=gcd(v)
        projective.append([c//g for c in v])
    rows=[];groups=[]
    for prime in [s['prime'] for s in d['signatures']]:
        F=GF(prime);e=EllipticCurve(F,[F(c) for c in E.a_invariants()])
        if e.discriminant()==0:raise ArithmeticError('bad reduction in proof')
        elements=e.points();key=lambda P:tuple(int(c) for c in P)
        doubles={key(2*P):2*P for P in elements};reps=[e(0)]
        masks={key(P):0 for P in doubles.values()}
        for P in elements:
            if key(P) in masks:continue
            bit=len(reps);extra=[R+P for R in reps]
            for i,R in enumerate(extra):
                for T in doubles.values():masks[key(R+T)]=i+bit
            reps+=extra
        if len(masks)!=len(elements) or len(reps) not in (1,2,4):raise ArithmeticError('quotient enumeration incomplete')
        dimension=ZZ(len(reps)).valuation(2)
        reduced=[e([F(c) for c in P]) for P in projective]
        rows.extend([[(masks[key(P)]>>i)&1 for P in reduced] for i in range(dimension)])
        groups.append({'prime':prime,'order':len(elements),'double_subgroup_order':len(doubles),'quotient_dimension':int(dimension)})
    M=Matrix(GF(2),rows);chosen=d['independent_column_indices']
    union_rank=int(M.rank());chosen_rank=int(M.matrix_from_columns(chosen).rank())
    if union_rank!=d['rank_lower_bound'] or chosen_rank!=len(chosen):
        raise ArithmeticError('independent enumerated quotient ranks differ')
    p=d['rank_certificate']['no_rational_2_torsion_prime'];F=GF(p);R=PolynomialRing(F,'x');x=R.gen()
    if not (x**3+F(E.a4())*x+F(E.a6())).is_irreducible():raise ArithmeticError('rational2-torsion exclusion failed')
    return {"input":str(INPUT.relative_to(ROOT)),"point_count":len(points),"union_finite_rank":union_rank,"exported_basis_rank":chosen_rank,"groups":groups,"no_two_torsion_prime":p}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    rows=[check_one(p) for p in INPUTS]
    result={'schema':'elliptic-curves.recent-outer26-sage.v1','status':'PASS','rows':rows,
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),*INPUTS]},
        'scope':'Standalone Sage point membership, complete finite groups and quotients by doubles, exported-basis independence and rational2-torsion exclusion. Only rank lower bounds; source search history is separately verified.'}
    if args.check:assert result==json.loads(OUT.read_text())
    else:
        with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print('PASS',[(r['point_count'],r['union_finite_rank']) for r in rows])
