#!/usr/bin/env sage-python
"""Standalone complete finite-group replay of retained26 full point clouds."""
import json,argparse
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,Matrix,lcm,gcd,PolynomialRing
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
    return {"input":INPUT.name,"point_count":len(points),"union_finite_rank":union_rank,"exported_basis_rank":chosen_rank,"groups":groups,"no_two_torsion_prime":p}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,nargs='+',required=True);args=p.parse_args()
 rows=[check_one(path) for path in args.input]
 print('PASS independent retained26 full clouds',[(r['point_count'],r['union_finite_rank']) for r in rows])
