#!/usr/bin/env sage-python
"""Portable exact specialized spans, including rational2-torsion cases."""
import argparse,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,PolynomialRing,matrix,lcm,gcd
def quotient(E,points,p,ell):
    F=GF(p);e=EllipticCurve(F,[F(c) for c in E.a_invariants()]);assert e.discriminant()!=0
    key=lambda P:tuple(int(c) for c in P)
    elements=e.points();multiples={key(ell*P):ell*P for P in elements};labels={key(P):0 for P in multiples.values()};reps=[e(0)]
    for P in elements:
        if key(P) in labels:continue
        size=len(reps);new=[]
        for k in range(1,ell):
            for j,R in enumerate(reps):
                Q=R+k*P;new.append(Q)
                for T in multiples.values():labels[key(Q+T)]=j+k*size
        reps+=new
    assert len(labels)==len(elements) and len(reps) in (1,ell,ell*ell)
    reduced=[]
    for P in points:
        den=lcm(c.denominator() for c in P);v=[ZZ(c*den) for c in P];div=gcd(v);reduced.append(e([F(c//div) for c in v]))
    dimension=ZZ(len(reps)).valuation(ell)
    return [[(labels[key(P)]//ell**j)%ell for P in reduced] for j in range(dimension)]
def main(path,out):
    data=json.loads(path.read_text());R=PolynomialRing(QQ,'T');K=R.fraction_field();results=[]
    for row in data['rows']:
        seed=row['original'];parent=data['parents'][seed['parent_index']];t=QQ(seed['fibre_T']);u=QQ(seed['scale_from_parent']);E=EllipticCurve(QQ,seed['curve'])
        assert E.a4()==u**4*R(parent['raw_A'])(t) and E.a6()==u**6*R(parent['raw_B'])(t)
        points=[E([QQ(c) for c in P]) for P in seed['points']]
        for P,j in zip(points,seed['generic_indices']):
            source=parent['sections'][j];assert P==E([u*u*K(source[0])(t),u**3*K(source[1])(t),K(source[2])(t)])
        proof=row['proof'];ell=proof['modulus'];indices=row['basis_indices'];rows=[]
        for sig in proof['signatures']:rows+=quotient(E,points,sig['prime'],ell)
        M=matrix(GF(ell),rows);assert M.rank()==len(indices) and M.matrix_from_columns(indices).rank()==len(indices)
        p=proof['torsion_prime'];ee=EllipticCurve(GF(p),[GF(p)(c) for c in E.a_invariants()]);assert ee.discriminant()!=0 and ee.cardinality()%ell!=0
        basis=[points[j] for j in indices]
        if row['relations'] is not None:
            assert len(row['relations'])==len(points)
            for r in row['relations']:
                residual=E([QQ(c) for c in r['torsion_remainder']]);assert 2*residual==E(0)
                assert r['multiplier']>0 and r['multiplier']*points[r['point_index']]==sum((ZZ(c)*P for c,P in zip(r['basis_word'],basis)),E(0))+residual
        else:assert indices==list(range(len(points)))
        results.append({'id':seed['id'],'status':'PASS','exact_displayed_span_rank':len(indices),'modulus':ell,'finite_primes':[s['prime'] for s in proof['signatures']]})
        print(seed['id'],'PASS EXACT DISPLAYED SPAN',len(indices),flush=True)
    result={'status':'PASS','rows':results,'scope':'Exact model/section specialization, all point equations, independently enumerated finite quotient groups and prime-to-ell torsion exclusion. Exact elliptic relations close the two smaller displayed spans. No whole-curve rank upper bounds or new point search.'}
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.input,a.output)
