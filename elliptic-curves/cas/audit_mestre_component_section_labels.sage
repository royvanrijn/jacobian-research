#!/usr/bin/env sage-python
"""Replay the old thirteen probes with coherent generic root labels and signs."""
import sys
from pathlib import Path
from sage.all import QQ,PolynomialRing,prod,GF,matrix
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import verify_mestre_fermigier_two_section_generic_rank13 as old
from icarm_curve245_mestre import fermigier_roots
from probe_mestre_fermigier_two_section_local_continuation import reconstructed_second_line,normalized_data
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import primitive_visible_points,quartic_point_to_short_jacobian
from research_runtime.store import checkpoint
OUT=ROOT/'artifacts/generated-results/elliptic-curves/mestre_component_label_audit_v1.json'
D=ROOT/'artifacts/local/elliptic-curves/mestre-parent-portfolio-intake-v2'

def main():
    if OUT.exists():raise FileExistsError('preserve section-label audit')
    U=PolynomialRing(QQ,'u');K=U.fraction_field();u=K.gen();R=PolynomialRing(K,'T');T=R.gen();Sring=PolynomialRing(R,'x');x=Sring.gen()
    v,c2,m2=reconstructed_second_line(u);c1,m1=normalized_data(u,v)[4:6]
    raw=fermigier_roots(u,v);roots=tuple((r-raw[0])/(raw[1]-raw[0]) for r in raw)
    print('SYMBOLIC ROOTS',flush=True)
    F=prod((x-r-T)*(x-r+T) for r in roots);S=x**6
    for j in range(5,-1,-1):S+=R((F[6+j]-(S*S)[6+j])/2)*x**j
    rem=S*S-F;assert rem.degree()<=4
    quartic=[]
    for j in range(5):
        q,res=rem[j].quo_rem(T*T);assert not res;quartic.append(q)
    def ev(f,value):return QQ(f.numerator()(value))/QQ(f.denominator()(value))
    ys=[]
    for c,m in ((c1,m1),(c2,m2)):
        xx=c+m*T;square=sum(q*xx**j for j,q in enumerate(quartic));degree=square.degree();assert degree%2==0
        y=R(square.leading_coefficient().sqrt())*T**(degree//2)
        for j in range(degree//2-1,-1,-1):y+=((square[degree//2+j]-(y*y)[degree//2+j])/(2*y.leading_coefficient()))*T**j
        assert y*y==square
        if sum(ev(a,QQ(-5)) for a in y)<0:y=-y
        ys.append(y)
    print('EXACT GLOBAL ORDINATE SQUARES',flush=True)
    order=sorted(range(6),key=lambda j:ev(roots[j],QQ(-5)))
    rows=[];original_rows=[];records=[]
    for pu,pt,p in old.PROBES:
        uu,tt=QQ(str(pu)),QQ(str(pt));normalized=[cert.F(str(ev(r,uu))) for r in roots];construction=SixRootMestreConstruction(tuple(normalized));visible=primitive_visible_points(construction,pt)
        indices=[2*construction.roots.index(normalized[j])+s for j in order for s in (0,1)]
        labelled=[visible[i] for i in indices]
        model,original=old.specialized_points(pu,pt);points=[quartic_point_to_short_jacobian(construction,pt,P) for P in labelled[:11]];signs=[]
        for i,y in enumerate(ys):
            value=sum(ev(a,uu)*tt**j for j,a in enumerate(y))/QQ(str(construction.quartic_square_scale))
            assert value!=0;sign=1 if value>0 else -1;signs.append(sign)
            P=original[12+i];points.append((P[0],sign*P[1]))
        original_sig=old.mod_l_reduction_signature(model,old.basis_points(original),p,modulus=3)
        sig=old.mod_l_reduction_signature(model,tuple(points),p,modulus=3)
        original_rows.extend(original_sig.rows);rows.extend(sig.rows)
        records.append({'outer_u':str(pu),'T':str(pt),'prime':p,'visible_indices_in_old_sorted_roster':indices[:11],'extra_ordinate_signs_relative_to_old_positive_sqrt':signs,'original_rows':[list(r) for r in original_sig.rows],'coherent_rows':[list(r) for r in sig.rows]})
    original_rank=int(matrix(GF(3),original_rows).rank());coherent_rank=int(matrix(GF(3),rows).rank());assert original_rank==13
    result={'schema':'elliptic-curves.mestre-component-label-audit.v1','status':'PASS','original_stacked_rank':original_rank,'coherent_stacked_rank':coherent_rank,'fixed_source_root_order':order,'normalized_root_functions':list(map(str,roots)),'global_extra_ordinate_coefficients':[[str(c) for c in y.list()] for y in ys],'probes':records,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),Path(old.__file__),CAS/'screen_mestre_fermigier_two_section_height_triage.py',CAS/'mestre_root_tuples.py',CAS/'probe_mestre_fermigier_two_section_local_continuation.py',CAS/'icarm_curve245_mestre.py',CAS/'nagao_1994.py']},'boundary':'Same13 historical finite probes, with one fixed generic root order and exact global rational ordinate branches. Different rank after restoring coherent labels invalidates the mixed-specialization independence argument; it does not alone prove a generic-rank upper bound. No elliptic point or parameter search.'}
    checkpoint(OUT,result);print('PASS old/coherent',original_rank,coherent_rank,flush=True)

if __name__=='__main__':main()
