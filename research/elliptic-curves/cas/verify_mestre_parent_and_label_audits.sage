#!/usr/bin/env sage-python
"""Independent full finite groups for parent seeds and coherent generic labels."""
import sys,json,argparse
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,PolynomialRing,Matrix,lcm,gcd
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
import verify_mestre_fermigier_two_section_generic_rank13 as old
from mestre_root_tuples import SixRootMestreConstruction
ART=ROOT/'artifacts/generated-results/elliptic-curves';PARENTS=ART/'mestre_parent_portfolio_intake_v1.json';LABELS=ART/'mestre_component_label_audit_v1.json';OUT=ART/'mestre_parent_and_label_independent_v1.json'

def quotient(curve,points,prime,ell):
    E=EllipticCurve(QQ,[QQ(str(c)) for c in curve]);rational=[E([QQ(str(x)),QQ(str(y))]) for x,y in points];F=GF(prime);e=EllipticCurve(F,[F(c) for c in E.a_invariants()]);assert e.discriminant()!=0
    key=lambda P:tuple(int(c) for c in P)
    elements=e.points();multiples={key(ell*P):ell*P for P in elements};reps=[e(0)];labels={key(P):0 for P in multiples.values()}
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
    for P in rational:
        den=lcm([c.denominator() for c in P]);v=[ZZ(c*den) for c in P];g=gcd(v);reduced.append(e([F(c//g) for c in v]))
    dimension=int(ZZ(len(reps)).valuation(ell));rows=[[(labels[key(P)]//ell**j)%ell for P in reduced] for j in range(dimension)]
    return rows,{'prime':prime,'group_order':len(elements),'multiple_subgroup_order':len(multiples),'quotient_dimension':dimension}

def expected():
    parents=cert.read(PARENTS);labels=cert.read(LABELS);assert parents['status']==labels['status']=='PASS';parent_rows=[]
    for r in parents['rows']:
        rows=[];groups=[]
        for sig in r['signatures']:
            bits,g=quotient(r['curve'],r['points'],sig['prime'],2);rows+=bits;groups.append(g)
        M=Matrix(GF(2),rows);indices=r['independent_column_indices'];assert M.rank()==r['rank_lower_bound'] and M.matrix_from_columns(indices).rank()==len(indices)
        p=r['rank_certificate']['no_rational_2_torsion_prime'];F=GF(p);R=PolynomialRing(F,'x');x=R.gen();assert (x**3+F(QQ(r['curve'][3]))*x+F(QQ(r['curve'][4]))).is_irreducible()
        parent_rows.append({'id':r['id'],'points':len(r['points']),'rank_lower_bound':int(M.rank()),'groups':groups})
    U=PolynomialRing(QQ,'u');K=U.fraction_field();roots=[K(s) for s in labels['normalized_root_functions']];ys=[[K(c) for c in row] for row in labels['global_extra_ordinate_coefficients']]
    def ev(f,value):return QQ(f.numerator()(value))/QQ(f.denominator()(value))
    original_rows=[];coherent_rows=[];probes=[]
    for record,(pu,pt,prime) in zip(labels['probes'],old.PROBES):
        assert record['outer_u']==str(pu) and record['T']==str(pt) and record['prime']==prime
        u,t=QQ(str(pu)),QQ(str(pt));specialized=[cert.F(str(ev(r,u))) for r in roots];construction=SixRootMestreConstruction(tuple(specialized))
        indices=[2*construction.roots.index(specialized[j])+s for j in labels['fixed_source_root_order'] for s in (0,1)]
        assert indices[:11]==record['visible_indices_in_old_sorted_roster']
        model,points=old.specialized_points(pu,pt);coherent=[points[i] for i in indices[:11]]
        for j,y in enumerate(ys):
            Y=sum(ev(a,u)*t**i for i,a in enumerate(y))/QQ(str(construction.quartic_square_scale));sign=1 if Y>0 else -1;assert Y!=0 and sign==record['extra_ordinate_signs_relative_to_old_positive_sqrt'][j]
            P=points[12+j];coherent.append((P[0],sign*P[1]))
        a,g=quotient(model,old.basis_points(points),prime,3);b,h=quotient(model,coherent,prime,3);assert g==h
        original_rows+=a;coherent_rows+=b;probes.append(g)
    old_rank=int(Matrix(GF(3),original_rows).rank());new_rank=int(Matrix(GF(3),coherent_rows).rank());assert old_rank==labels['original_stacked_rank']==13 and new_rank==labels['coherent_stacked_rank']==11
    curve,points=old.specialized_points(cert.F(-5),cert.F(1));_,torsion=quotient(curve,[],19,3);assert torsion['group_order']==28
    paths=[Path(__file__).resolve(),PARENTS,LABELS,Path(old.__file__),CAS/'mestre_root_tuples.py']
    return {'status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'parents':parent_rows,'old_mixed_column_rank':old_rank,'coherent_column_rank':new_rank,'generic_three_torsion_exclusion':torsion,'label_probe_groups':probes,'scope':'Independent complete finite-group enumeration and quotient cosets, exact point membership and torsion exclusion. Six14-point seed clouds have rank lower bound11. The same old13 generic probes give13 with mixed labels but11 with fixed labels and globally rational ordinate branches. Surface-count and symbolic-square identities are separately certified; no generic-rank upper bound.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=expected()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve independent replay')
        checkpoint(OUT,result)
    print('PASS six parent11 seeds; mixed/coherent13/11')
