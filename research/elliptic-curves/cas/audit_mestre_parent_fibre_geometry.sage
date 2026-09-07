#!/usr/bin/env sage-python
"""Exact degree-two quotient geometry of the six already constructed parents."""
import argparse,json,hashlib
from pathlib import Path
from sage.all import QQ,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
SOURCE=ART/'mestre_parent_portfolio_intake_v1.json'
OUT=ART/'mestre_parent_fibre_geometry_v1.json'

def compute():
    R=PolynomialRing(QQ,'T');T=R.gen();Sring=PolynomialRing(QQ,'s');s=Sring.gen();rows=[]
    for row in json.loads(SOURCE.read_text())['rows']:
        A=R(row['A_coefficients']);B=R(row['B_coefficients'])
        if any(A[i] for i in range(1,9,2)) or any(B[i] for i in range(1,13,2)):raise ArithmeticError('quadratic base-change parity failed')
        a=Sring([A[2*i] for i in range(5)]);b=Sring([B[2*i] for i in range(7)])
        if R(a(T*T))!=A or R(b(T*T))!=B:raise ArithmeticError('exact s=T^2 identity failed')
        delta=-16*(4*a**3+27*b**2);c4=-48*a
        profile=sorted([int(f.degree()),int(m)] for f,m in delta.squarefree_decomposition())
        if a.degree()!=4 or b.degree()!=6 or delta.degree()!=10 or profile!=[[1,2],[8,1]] or c4.gcd(delta).degree()!=0:raise ArithmeticError('rational quotient semistable profile differs')
        double=next(f.monic() for f,m in delta.squarefree_decomposition() if m==2)
        base=-double[0];pullback=R(T*T-base)
        if not base or delta(0)==0:raise ArithmeticError('branch at zero must be smooth')
        # Both branches of the multiplicative nodes must be distinguished
        # before counting rational fibre-component classes.
        P=PolynomialRing(QQ,'x');x=P.gen()
        finite_cubic=x**3+a(base)*x+b(base);node=finite_cubic.gcd(finite_cubic.derivative()).monic()
        if node.degree()!=1:raise ArithmeticError('finite I2 node not unique')
        node_x=-node[0];finite_split=bool((3*node_x).is_square())
        infinity_cubic=x**3+a[4]*x+b[6];infnode=infinity_cubic.gcd(infinity_cubic.derivative()).monic()
        if infnode.degree()!=1 or not (3*(-infnode[0])).is_square():raise ArithmeticError('split I2 infinity node required')
        rational_base_values=2 if base.is_square() else 0
        # An I2 fibre has one nonidentity geometric component, independent
        # of its split sign. Its Galois orbit contributes one Q divisor.
        finite_orbits=2 if rational_base_values else 1
        rows.append({'id':row['id'],'outer_u':row['outer_u'],
          'quotient_A_coefficients':list(map(str,a.list())),'quotient_B_coefficients':list(map(str,b.list())),
          'quotient_discriminant_profile':profile,'quotient_fibres':'8 I1 + I2 + split I2 at infinity',
          'quotient_geometric_mw_rank':6,'finite_double_base_s':str(base),'finite_I2_split_over_Q':finite_split,
          'parent_double_fibre_polynomial':list(map(str,pullback.list())),
          'parent_double_fibres_rational':bool(rational_base_values),'parent_nonidentity_component_Q_orbits':3+finite_orbits,
          'parent_geometric_NS_rank_lower_bound':18,'parent_arithmetic_NS_rank_lower_bound':2+11+3+finite_orbits})
    return {'schema':'elliptic-curves.mestre-parent-fibre-geometry.v1','status':'PASS','rows':rows,
      'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),SOURCE)},
      'scope':'Exact s=T^2 descent of six existing K3 fibrations to degree(4,6) rational elliptic surfaces with8I1+2I2, hence geometric MW rank6 by Shioda-Tate. The base change has16I1+2I2+split I4. Together with the already certified11 rational sections this gives NS lower bounds only. It is not a new parent, new fibration, full NS calculation, generic K3 rank upper bound, or arithmetic foundry admission. No rootless-frame or coefficient search is performed.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute()
    if a.check:assert result==json.loads(OUT.read_text())
    else:
        with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print('PASS6 QUADRATIC BASE CHANGES',[(r['outer_u'],r['parent_double_fibres_rational'],r['parent_arithmetic_NS_rank_lower_bound']) for r in result['rows']])
