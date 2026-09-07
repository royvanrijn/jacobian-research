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
    data=json.loads(path.read_text());results=[]
    for cloud in data['clouds']:
        E=EllipticCurve(QQ,cloud['curve']);points=[E([QQ(c) for c in P]) for P in cloud['points']];audits=[]
        for audit in cloud['audits']:
            ell=audit['modulus'];rows=[]
            for sig in audit['signatures']:rows+=quotient(E,points,sig['prime'],ell)
            M=matrix(GF(ell),rows);indices=audit['independent_indices'];rank=M.rank();assert rank==audit['rank_lower_bound']==len(indices) and M.matrix_from_columns(indices).rank()==rank
            p=audit['no_rational_ell_torsion_prime'];e=EllipticCurve(GF(p),[GF(p)(c) for c in E.a_invariants()]);assert e.discriminant()!=0 and e.cardinality()%ell!=0
            audits.append({'modulus':ell,'rank_lower_bound':int(rank),'finite_primes':[s['prime'] for s in audit['signatures']]})
        assert max(a['rank_lower_bound'] for a in audits)==cloud['rank_lower_bound']
        results.append({'id':cloud['id'],'point_count':len(points),'initial_rank':cloud['initial_rank'],'rank_lower_bound':cloud['rank_lower_bound'],'discovered_rank_gain_lower_bound':cloud['rank_lower_bound']-cloud['initial_rank'],'audits':audits})
        print('PASS full cloud',cloud['id'],len(points),'points; rank >=',cloud['rank_lower_bound'],flush=True)
    with out.open('x') as f:json.dump({'status':'PASS','rows':results,'scope':'Every retained rational point checked on its exact equation; complete finite groups and quotients at odd primes independently enumerated, torsion excluded at each proof modulus, and exported independent columns verified. Rank lower bounds only; exposure, provenance and exact initial-span upper bounds are separate dependencies.'},f,indent=2);f.write('\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.input,a.output)
