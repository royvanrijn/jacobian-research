#!/usr/bin/env python3
"""Fixed algebraic descent gate; no parameter search or Jacobian descent."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,matrix,prod
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_native_genus_five_lift_inputs_v1.json'
VERIFIED=r.OUT/'rank_jump_native_genus_five_lift_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_native_genus_two_lift_gate_v1.json'


def compute():
    inp=r.read(INPUT);verified=r.read(VERIFIED)
    for data in (inp,verified):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');qs=[R(c['form']) for c in inp['covers']];h,f,g=qs;sextic=prod(qs)
    assert sextic.degree()==6 and sextic.is_squarefree()
    assert all(not q.discriminant().is_square() for q in qs)
    # Native deck group is F2^3; finite inertia generators are 001,010,100.
    inertia={1,2,4};subgroups=[]
    for mask in range(256):
        H={i for i in range(8) if (mask>>i)&1}
        if 0 not in H or any(a^b not in H for a in H for b in H):continue
        if not H.intersection(inertia):subgroups.append(sorted(H))
    maximum=max(map(len,subgroups));largest=[H for H in subgroups if len(H)==maximum]
    assert maximum==4 and largest==[[0,3,5,6]]
    # The remaining character is the total product: (a,u,v) -> a*u*v.
    assert all((x.bit_count()%2)==0 for x in largest[0])
    # Nonzero pair classes in J(H)[2], modeled by even branch subsets mod complement.
    pair_masks=[0b000011,0b001100,0b110000]
    canonical=lambda m:min(m,m^0b111111)
    assert len({canonical(x) for x in pair_masks})==3
    assert canonical(pair_masks[0]^pair_masks[1]^pair_masks[2])==0
    assert all((pair_masks[i]&pair_masks[j]).bit_count()%2==0 for i in range(3) for j in range(i))
    det=matrix(QQ,[q.list() for q in qs]).det();assert det
    brackets=[qs[(i+1)%3].derivative()*qs[(i+2)%3]-qs[(i+1)%3]*qs[(i+2)%3].derivative() for i in range(3)]
    assert all(q.degree()==2 and q.is_squarefree() for q in brackets)
    assert prod(brackets).is_squarefree()
    t0=QQ(inp['retained_lift']['t']);roots=list(map(QQ,inp['retained_lift']['roots']));y0=prod(roots)
    assert y0*y0==sextic(t0)
    rows=[]
    for p in (131,137):
        F=GF(p);rp=PolynomialRing(F,'t');qred=[rp(q) for q in qs]
        assert rp(sextic).degree()==6 and rp(sextic).is_squarefree()
        base_count=0;liftable=0;lift_count=0;branch_base=0;branch_liftable=0
        for x in list(F)+[None]:
            vals=[q(x) if x is not None else q[2] for q in qred]
            n=[1 if v==0 else 2 if v.is_square() else 0 for v in vals]
            product=prod(vals);ny=1 if product==0 else 2 if product.is_square() else 0
            base_count+=ny;lift_count+=prod(n)
            if product!=0:
                # On H the third squareclass follows from the first two.
                good=ny and vals[0].is_square() and vals[1].is_square()
                assert bool(good)==bool(prod(n))
            else:
                assert sum(v==0 for v in vals)==1
                # At a branch point use the two nonvanishing factors.
                good=all(v.is_square() for v in vals if v!=0)
                branch_base+=ny;branch_liftable+=int(bool(good))
            if good:liftable+=ny
        assert lift_count==4*liftable
        prior=next(row for row in verified['finite_field_character_checks'] if row['prime']==p)
        assert base_count==prior['character_quotient_points'][3] and lift_count==prior['triple_points']
        rows.append({'prime':p,'genus_two_points':int(base_count),'liftable_genus_two_points':int(liftable),
            'genus_five_points':int(lift_count),'branch_base_points':int(branch_base),'liftable_branch_points':int(branch_liftable)})
    return {'schema':'rank-jump.native-genus-two-lift-gate.v1','status':'PASS','layer':'solubility',
        'sextic_coefficients':list(map(str,sextic.list())),
        'retained_genus_two_point':{'t':str(t0),'y':str(y0)},
        'native_deck_subgroups_acting_freely':subgroups,'unique_maximal_free_native_subgroup':largest[0],
        'etale_cover_degree':4,'quotient_genus':2,'rational_isotropic_2_torsion_dimension':2,
        'quadratic_coefficient_determinant':str(det),'Richelot_brackets':[[str(x) for x in q.list()] for q in brackets],
        'finite_field_lift_checks':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,VERIFIED,Path(__file__),HERE/'retrospective.py')},
        'boundary':'Algebraic quotient and lift checks only. The note proves the Abel-Jacobi/isogeny pullback interpretation. No complete genus-two rational-point computation, genus-two rank bound, or prospective parameter selection.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=compute()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert r.read(OUTPUT)==data
    print('PASS: unique maximal free native quotient genus2; etale degree4; squareclass gate')
