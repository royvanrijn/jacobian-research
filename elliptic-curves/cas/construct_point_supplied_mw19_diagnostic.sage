#!/usr/bin/env sage-python
"""One fixed MW18 -> MW19 V4 incidence construction, with arithmetic costs.

The retained 531/12f61 cover is paired with the lexically first other cover
at that anchor, 0a9bf. No pair census, Selmer computation or point search.
An explicit point proves positive rank of the actual paired genus-one base;
the second cover supplies an independent section by its distinct character.
"""
import json
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path
import sys
from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, GF, matrix, prime_range

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint
import mw16_model_size as sizes

ART=ROOT/'artifacts/generated-results'
COVERS=ART/'elkies-k3-r17-extreme-anchored-mw18-covers-v1.json'
HELPER=ROOT/'elkies-k3/scripts/prepare_r17_norm12_11952_alternate_v4_rank_one_bases.sage'
GEOMETRY=ART/'elliptic-curves/mw18_generic_height_geometry_v1.json'
OUT=ART/'elliptic-curves/point_supplied_mw19_diagnostic_v1.json'


def run():
    helpers=SourceFileLoader('paired_base_helpers',str(HELPER)).load_module()
    cert=json.loads(COVERS.read_text())
    fibre=next(f for ch in cert['charts'] for f in ch['fibres'] if f['curve_id']==531)
    left=next(c for c in fibre['covers'] if c['label']=='08234-orbit-12f61')
    right=min((c for c in fibre['covers'] if c['label']!=left['label']),key=lambda c:c['label'])
    assert right['label']=='08234-orbit-0a9bf'
    R=PolynomialRing(QQ,'r');r=R.gen();K=R.fraction_field()
    q1=R(left['branch_quadratic_coefficients_low_to_high']);q2=R(right['branch_quadratic_coefficients_low_to_high'])
    assert q1.degree()==q2.degree()==2 and q1.is_squarefree() and q2.is_squarefree() and q1.gcd(q2)==1
    t0=QQ(fibre['native_parameter']);s0=QQ(left['canonical_positive_square_root']);w0=QQ(right['canonical_positive_square_root'])
    assert s0*s0==q1(t0) and w0*w0==q2(t0)
    den=1-q1[2]*r*r
    tmap=K(t0+(q1.derivative()(t0)*r*r-2*s0*r)/den)
    smap=K(s0+(tmap-t0)/r)
    assert smap*smap==q1(tmap)
    quartic=R(q2(tmap)*den**2);assert quartic.degree()==4 and quartic.is_squarefree() and quartic[0]==w0*w0
    _,_,base,opp,constants=helpers.pointed_curve(list(quartic),QQ(0),w0,'r')
    candidates=sorted((p for p in base.lift_x(opp,all=True) if p[1]),key=lambda p:tuple(p.xy()))
    if not candidates:raise ArithmeticError('fixed opposite-anchor construction supplies no base point')
    P=candidates[0]
    # Reduction is injective on rational torsion at good primes. Their group
    # order gcd is an annihilator; an exact nonzero multiple proves nontorsion.
    reductions=[];annihilator=0
    for p in prime_range(5,100):
        try:
            ep=EllipticCurve(GF(p),[GF(p)(a) for a in base.ainvs()])
            if not ep.discriminant():continue
        except (ZeroDivisionError,ArithmeticError,ValueError):continue
        order=int(ep.cardinality());annihilator=gcd(annihilator,order)
        reductions.append(dict(prime=int(p),order=order))
        if len(reductions)>=3:break
    assert len(reductions)==3 and annihilator and annihilator*P
    geometry=json.loads(GEOMETRY.read_text())
    gl=next(c for c in geometry['covers'] if c['cover_label']==left['label'])
    gr=next(c for c in geometry['covers'] if c['cover_label']==right['label'])
    G=matrix(QQ,gl['old_gram']); assert G==matrix(QQ,gr['old_gram'])
    from sage.all import vector
    v1,v2=vector(QQ,gl['trace_word']),vector(QQ,gr['trace_word'])
    Gram=(4*G).augment(matrix(QQ,17,1,list(2*G*v1))).augment(matrix(QQ,17,1,list(2*G*v2)))
    cross=v1*G*v2
    Gram=Gram.stack(matrix(QQ,1,19,list(2*G*v1)+[16,cross]))
    Gram=Gram.stack(matrix(QQ,1,19,list(2*G*v2)+[cross,16]))
    assert Gram.is_positive_definite() and Gram.det()==4**17*948*36
    # A small fixed list of base multiples measures the height of the
    # constructed t and point. These are not selected by any rank/score.
    directpath=ROOT/gl['model_source'];direct=json.loads(directpath.read_text())
    A=R(direct['weierstrass_model']['A_coefficients_low_to_high']);B=R(direct['weierstrass_model']['B_coefficients_low_to_high'])
    rows=[]
    for k in (2,3):
        image=helpers.inverse_pointed(k*P,constants)
        if image is None:
            rows.append(dict(multiple=k,status='INVERSE_MAP_EXCEPTION'));continue
        parameter,ordinate=image
        assert ordinate**2==quartic(parameter)
        if den(parameter)==0:
            rows.append(dict(multiple=k,status='BASE_PARAMETER_POLE'));continue
        t=tmap(parameter);u1=smap(parameter);u2=ordinate/den(parameter)
        assert u1*u1==q1(t) and u2*u2==q2(t)
        E=EllipticCurve(QQ,[A(t),B(t)])
        sections=[]
        for c,u in ((left,u1),(right,u2)):
            s=c['eighteenth_section']
            x=R(s['x0_coefficients_low_to_high'])(t)+R(s['x1_coefficients_low_to_high'])(t)*u
            y=R(s['y0_coefficients_low_to_high'])(t)+R(s['y1_coefficients_low_to_high'])(t)*u
            sections.append(E(x,y))
        rows.append(dict(multiple=k,status='EXACT_FIBRE_AND_TWO_POINTS',base_point=sizes.point_record(k*P),
            paired_parameter=str(parameter),t=str(t),u1=str(u1),u2=str(u2),curve=list(map(str,E.ainvs())),
            supplied_points=[sizes.point_record(p) for p in sections],
            parameter_height_bits=int(max(abs(t.numerator()),t.denominator()).nbits()),
            j_log2_height=float(E.j_invariant().global_height()/QQ(2).log()),
            raw_model_size=sizes.size(E.ainvs()),point_sizes=[sizes.size(p.xy()) for p in sections],
            specialized_independence='NOT_TESTED; generic rank 19 does not certify this particular fibre'))
    output=dict(schema='elliptic-curves.point-supplied-mw19-diagnostic.v1',
        status='EXACT_INDEPENDENT_SECTION_AND_POSITIVE_RANK_BASE',covers=[left['label'],right['label']],
        inputs={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in (COVERS,HELPER,GEOMETRY,directpath,Path(__file__))},
        quadratic_branches=[list(map(str,q)) for q in (q1,q2)],
        paired_quartic=list(map(str,quartic)),base_curve=list(map(str,base.ainvs())),
        base_point=sizes.point_record(P),base_torsion_annihilator=annihilator,
        base_good_reductions=reductions,nonzero_annihilator_multiple=sizes.point_record(annihilator*P),
        base_genus=1,base_rank_lower_bound=1,base_degree_over_old_parameter=4,
        generic_rank_lower_bound=19,gain_relative_to_retained_MW18=1,
        generic_gram=[list(map(str,row)) for row in Gram.rows()],generic_gram_determinant=str(Gram.det()),
        new_section_orthogonal_height_over_paired_base='6',
        rank_accounting='The same rank-28 anchor has 9 remaining displayed directions beyond these 19, '
        'rather than 10 beyond MW18. Its rank has not increased.',
        height_samples=rows,
        claim_boundary='One fixed standard V4 construction from retained covers, not a novel rank-19 theorem, '
        'a record search, a saturated generic lattice, or a specialized rank certificate. '
        'The actual paired genus-one base is used; the product-quartic quotient alone would not give the same rational lifts.')
    checkpoint(OUT,output);print('CONSTRUCTION',output['status'],'base_rank>=1','generic_rank>=19',flush=True)
    for row in rows:print('COST',row['multiple'],row['status'],row.get('parameter_height_bits'),row.get('raw_model_size'),flush=True)


if __name__=='__main__':run()
