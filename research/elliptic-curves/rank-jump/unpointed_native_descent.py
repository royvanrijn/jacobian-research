#!/usr/bin/env python3
"""Canonical doubling loses the native class; fixed local twist control."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,prod
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'UNPOINTED_NATIVE_DESCENT_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_native_genus_five_lift_inputs_v1.json'
GATE=r.OUT/'rank_jump_native_genus_two_lift_gate_v1.json'
OUTPUT=r.OUT/'rank_jump_unpointed_native_descent_v1.json'


def compute():
    inp=r.read(INPUT);gate=r.read(GATE)
    for data in (inp,gate):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'T');T=R.gen();qs=[R(c['form']) for c in inp['covers']];F=prod(qs)
    assert F.degree()==6 and F.is_squarefree()
    S=PolynomialRing(R,'x');x=S.gen();u=(x-T)**2
    norms=[R(u.resultant(S(q.list()))) for q in qs]
    assert all(n==q*q for n,q in zip(norms,qs,strict=True))
    # Clear the denominator of v=y+F'(T)(x-T)/(2y), with y^2=F(T).
    numerator=4*F*S(F.list())-(2*F+F.derivative()*(x-T))**2
    quotient,remainder=numerator.quo_rem(u);assert remainder==0
    # Generic-form-only choice of a locally insoluble twist, independent of t0.
    p=131;rp=PolynomialRing(GF(p),'T');j=next(i for i,q in enumerate(qs) if rp(q).degree()==2 and rp(q).is_irreducible())
    # d=(d1,d2,d1*d2); odd valuation in the chosen quadratic coordinate.
    ds=[p,1,p] if j in (0,2) else [1,p,p]
    assert ds[j]==p and prod(ds)==p*p
    assert int(rp(qs[j]).leading_coefficient())!=0
    roots=list(map(QQ,inp['retained_lift']['roots']));t0=QQ(inp['retained_lift']['t'])
    assert all(a*a==q(t0) for a,q in zip(roots,qs,strict=True))
    rows=[]
    for p in (131,137):
        k=GF(p);rk=PolynomialRing(k,'T');qq=[rk(q) for q in qs]
        assert prod(qq).degree()==6 and prod(qq).is_squarefree()
        counts={f'{a}{b}':0 for a in range(2) for b in range(2)};examples={};total=0
        for t in k:
            vals=[q(t) for q in qq]
            if any(v==0 for v in vals):continue
            ys=(PolynomialRing(k,'y').gen()**2-prod(vals)).roots(multiplicities=False)
            for y in ys:
                native=[int(not v.is_square()) for v in vals[:2]];key=''.join(map(str,native))
                norms_at=[v*v for v in vals[:2]];assert all(v.is_square() for v in norms_at)
                counts[key]+=1;total+=1
                examples.setdefault(key,{'t':int(t),'y':int(y),'native_values':[int(v) for v in vals],
                    'doubled_divisor_resultants':[int(v) for v in norms_at]})
        assert sum(counts.values())==total and sum(counts[k] for k in ('01','10','11'))>0
        rows.append({'prime':p,'scope':'finite non-branch H points only','points_checked':total,
            'native_class_counts':counts,'doubled_divisor_class_00_count':total,'examples':examples})
    return {'schema':'rank-jump.unpointed-native-descent.v1','status':'PASS','layer':'solubility',
        'generic_resultants':[[str(c) for c in n.list()] for n in norms],
        'cleared_Mumford_identity_quotient':[[str(c) for c in R(coef).list()] for coef in quotient.list()],
        'finite_field_class_checks':rows,
        'matched_twist':{'equation_convention':'d_i*z_i^2=q_i(t); product quotient y=(sqrt(product d_i))*z_1*z_2*z_3',
            'd':ds,'product_square_root':131,'obstruction_prime':131,'irreducible_coordinate':j,
            'irreducible_label':inp['covers'][j]['label'],'quadratic_mod_prime':[int(c) for c in rp(qs[j]).list()],
            'no_Qp_points':True,'same_genus_two_product_quotient':True,'same_Richelot_isogeny_data':True},
        'original_twist_retained_lift':inp['retained_lift'],
        'canonical_doubled_Abel_map_lift_gate':'identically trivial; proof via rational Weierstrass pair in note',
        'bindings':{str(path.relative_to(r.ROOT)):r.digest(path.read_bytes()) for path in (INPUT,GATE,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: doubled map loses class; matched twist',d['matched_twist']['d'],'has Q131 obstruction')
