#!/usr/bin/env sage-python
"""Replay recovered 245 parameters through twelve generic section directions.

The old embedding is loaded only after constructing the family and points.
This checks the actual covariant-image subgroup, not just its rational span.
It does not certify that these twelve sections are a saturated generic basis.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, prod

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from icarm_curve245_mestre import fermigier_roots, fermigier_extra_line
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import primitive_visible_points, quartic_point_to_short_jacobian
from icarm_curve245 import POINTS, GENERAL_WEIERSTRASS_COEFFICIENTS, short_coefficients

SEARCH=ROOT/'artifacts/generated-results/elliptic-curves/curve302_inverse_fermigier_245_h8_v1.json'
TRUTH=ROOT/'artifacts/generated-results/elliptic-curves/latent_lattice_calibration_truth_v1.json'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_inverse_fermigier_245_mw_control_v1.json'


def digest(p): return sha256(p.read_bytes()).hexdigest()


def build():
    search=json.loads(SEARCH.read_text())
    assert search['status']=='COMPLETE' and search['control_recovered'] is True
    E=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)))
    reconstructed=[]
    for hit in search['hits']:
        roots=list(map(QQ,hit['roots']))
        construction=SixRootMestreConstruction(tuple(Fraction(str(r)) for r in roots))
        RT=PolynomialRing(QQ,'T'); T=RT.gen()
        RX=PolynomialRing(RT,'X'); X=RX.gen()
        # Independent old rational arithmetic at seven interpolation values.
        samples=[construction.primitive_quartic_coefficients(Fraction(i)) for i in range(1,8)]
        coeff=[RT.lagrange_polynomial([(QQ(i+1),QQ(samples[i][j])) for i in range(7)]) for j in range(5)]
        quartic=sum(coeff[i]*X**i for i in range(5))
        product=prod((X-r-T)*(X-r+T) for r in roots)
        g=X**6
        for k in range(5,-1,-1): g+=(product[6+k]-(g*g)[6+k])/2*X**k
        assert g*g-product==T*T*QQ(construction.quartic_content)*quartic
        for match in hit['matches']:
            t0=QQ(match['T_canonical'])
            if t0<=0: continue
            rational_t=Fraction(str(t0))
            visible=list(primitive_visible_points(construction,rational_t))
            for origin in match['native_parameters']:
                u,v=map(Fraction,(origin['u'],origin['v']))
                native=list(map(QQ,fermigier_roots(u,v)))
                scale=t0/QQ(origin['T'])
                intercept,slope=map(QQ,fermigier_extra_line(u,v))
                # Positive scale anchors at min(root); negative at max(root).
                offset=min(native) if scale>0 else max(native)
                extra_x=scale*(intercept-offset)+slope*T
                ordinate_square=RT(quartic(extra_x))
                assert ordinate_square.is_square()
                extra_y=ordinate_square.sqrt()
                assert extra_y**2==quartic(extra_x)
                qpts=visible+[(Fraction(str(extra_x(t0))),Fraction(str(extra_y(t0))))]
                jacpts=[quartic_point_to_short_jacobian(construction,rational_t,pt) for pt in qpts]
                J=EllipticCurve(QQ,list(map(QQ,construction.primitive_jacobian_coefficients(rational_t))))
                for iso in J.isomorphisms(E):
                    pts=[iso(J(list(map(QQ,pt)))) for pt in jacpts]
                    assert sum(pts[:12],E(0))==E(0)
                    reconstructed.append({'parameters':origin,'roots':roots,'t0':t0,
                        'extra_x':extra_x,'extra_y':extra_y,'iso':iso,
                        'points':pts[:11]+[pts[12]]})
    # This is the first access to the calibration subgroup.
    truth=next(x for x in json.loads(TRUTH.read_text())['fermigier_family_controls']
               if x['label'].startswith('ICARM_245_'))
    public=[E(list(map(QQ,p))) for p in POINTS]
    # Fresh complete finite-group quotients, independent of the old height
    # relation finder. A rank20 mod2 image and E(Q)[2]=0 prove independence.
    reduction_rows=[]; reduction_primes=[]
    def key(p): return None if p.is_zero() else (int(p[0]),int(p[1]))
    for p in [11,23,29,41,47,59,61,67,73,83,97,101,113,127,139,149]:
        field=GF(p); ep=EllipticCurve(field,[field(a) for a in E.a_invariants()])
        group=ep.points(); doubles={key(2*q):2*q for q in group}
        labels={k:0 for k in doubles}; cosets=1
        for q in group:
            if key(q) not in labels:
                for h in doubles.values(): labels[key(q+h)]=cosets
                cosets+=1
        assert cosets in [1,2,4] and len(labels)==len(group)
        reduced=[ep([field(q[0]),field(q[1])]) for q in public]
        for bit in range(cosets.bit_length()-1):
            reduction_rows.append([(labels[key(q)]>>bit)&1 for q in reduced])
        reduction_primes.append(p)
    assert matrix(GF(2),reduction_rows).rank()==20
    sc=list(map(QQ,short_coefficients())); f7=GF(7)
    assert all(x**3+f7(sc[3])*x+f7(sc[4]) for x in f7)
    embedding=matrix(ZZ,truth['embedding_matrix_rows'])
    old=[sum((n*p for n,p in zip(c,public)),E(0)) for c in embedding.columns()]
    matching=[]
    for candidate in reconstructed:
        selected=[]
        for pt in candidate['points']:
            options=[(j,sign) for j,q in enumerate(old) for sign in [1,-1] if pt==sign*q]
            if len(options)!=1: break
            selected.append(options[0])
        if len(selected)==12 and len({j for j,s in selected})==12:
            matching.append((candidate,selected))
    assert matching,'reconstructed generic subgroup failed exact signed-basis comparison'
    candidate,selected=matching[0]
    new=matrix(ZZ,20,12,lambda i,j:selected[j][1]*embedding[i,selected[j][0]])
    assert new.rank()==12
    assert new.transpose().hermite_form()==embedding.transpose().hermite_form()
    diagonal=[abs(int(v)) for v in new.smith_form()[0].diagonal() if v]
    assert prod(diagonal)==2048
    return {'schema':'curve302.inverse-fermigier-mw-control.v1',
        'status':'PASS_EXACT_245_FAMILY_AND_TRANSPORTED_RANK12_SUBGROUP',
        'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),SEARCH,TRUTH,
            Path(__file__).with_name('icarm_curve245_mestre.py'),Path(__file__).with_name('nagao_1994.py'),
            Path(__file__).with_name('mestre_root_tuples.py'),Path(__file__).with_name('icarm_curve245.py')]},
        'recovered_parameters':candidate['parameters'],'T_canonical':str(candidate['t0']),
        'roots':list(map(str,candidate['roots'])),
        'extra_x_coefficients':list(map(str,candidate['extra_x'].list())),
        'extra_y_coefficients':list(map(str,candidate['extra_y'].list())),
        'to_public_u_r_s_t':list(map(str,candidate['iso'].tuple())),
        'generic_quartic_identity_checked':True,'generic_extra_section_identity_checked':True,
        'visible_section_relation_checked':True,'exact_embedding_matrix_rows':[list(map(int,r)) for r in new.rows()],
        'rank':12,'smith_factors':diagonal,'index_in_displayed_primitive_closure':2048,
        'equals_retained_generic_image_as_integer_subgroup':True,
        'public_point_independence':{'method':'complete finite elliptic groups modulo doubles',
            'primes':reduction_primes,'rows_mod2':reduction_rows,'rank':20,
            'irreducible_two_division_cubic_prime':7},
        'matching_presentations':len(matching),
        'boundary':'Retrospective equation/parameter and actual transported section-subgroup calibration, with fresh exact public-point independence. No full generic saturation, new curve rank or parent of302 is asserted.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=build()
    if args.check: assert result==json.loads(OUT.read_text())
    else: OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'])
