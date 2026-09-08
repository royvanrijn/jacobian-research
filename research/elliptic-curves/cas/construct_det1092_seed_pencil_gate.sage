#!/usr/bin/env sage-python
"""One symbolic pencil gate and eight generic-only splitting witnesses.

No point search, new parameter addresses, V3 input, Selmer or class group.
The historical u0 remains a labelled calibration, never a blind selector.
"""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_seed_pencil_gate_v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def retain(p, data):
    text = json.dumps(data, indent=2, sort_keys=True)+'\n'
    if p.exists():
        assert p.read_text() == text
    else:
        p.write_text(text)


def build():
    paths = [ART/'curve302_recovered_mw17_parent_v1.json',
             ART/'det1092_first_centre_rr_net_v1.json',
             ART/'det1092_first_witness_pencil_genus_gate_v1.json',
             ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent, net, pencil, roster = [json.loads(p.read_text()) for p in paths]
    protocol = {'classification': 'frozen equation diagnostic, not a new search policy',
        'members': 'Entire already selected pencil B+(u0+v*t)A; no v is optimized.',
        'generic_witness_rule': 'v(t)=(-B(S0(t))/A(S0(t))-u0)/t; section0 only.',
        'cases': [{'index': i, 'parameter': r['parameter']} for i,r in enumerate(roster['cases'])],
        'limits': {'wall_seconds': 25, 'symbolic_pencils': 1, 'existing_addresses': 9,
                   'new_addresses': 0, 'point_searches': 0, 'V3_inputs': 0,
                   'Selmer_runs': 0, 'class_group_runs': 0, 'pilot_changes': 0},
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'protocol.json',protocol)
    R = PolynomialRing(QQ,'t')
    t = R.gen()
    K = R.fraction_field()

    def rat(row):
        return K(R(row['numerator']))/R(row['denominator'])

    def rec(f):
        f = K(f)
        return {'numerator': list(map(str,f.numerator().list())),
                'denominator': list(map(str,f.denominator().list()))}

    E = EllipticCurve(K,[rat(r) for r in parent['a_invariants']])
    basis = [E([rat(r) for r in p]) for p in parent['basis_weierstrass_coordinates']]
    w = -vector(QQ,net['trace_word'])
    Z = sum((n*P for n,P in zip(w,basis)),E(0))
    h = R(Z[0].denominator().sqrt())
    A,B = [[R(r) for r in net[k]] for k in ['A','B']]
    u0 = QQ(pencil['u0'])
    kappa = QQ((A[2]*B[1]-A[1]*B[2])/h**3)
    assert kappa
    V = PolynomialRing(K,'v')
    v = V.gen()
    line = [V(B[j]+u0*A[j])+t*A[j]*v for j in range(3)]
    f0,f1,f2 = line
    mnum = -f1+E.a1()*f2/2
    cx = Z[0]+E.b2()/12
    cy = Z[1]+(E.a1()*Z[0]+E.a3())/2
    a = -E.c4()/48
    F = (mnum**4-6*cx*mnum*mnum*f2*f2-8*cy*mnum*f2**3
         -(3*cx*cx+4*a)*f2**4)/h**6
    assert F.degree() == 4
    assert mnum[1]*f2[0]-mnum[0]*f2[1] == kappa*t*h**3
    c0,c1,c2,c3,c4 = F.list()
    I = 12*c4*c0-3*c3*c1+c2*c2
    J = 72*c4*c2*c0+9*c3*c2*c1-27*c4*c1*c1-27*c3*c3*c0-2*c2**3
    assert I == (kappa*t)**4*E.c4()
    assert J == 2*(kappa*t)**6*E.c6()
    assert F.discriminant() == 256*(kappa*t)**12*E.discriminant()
    # Source-only witness, defined before any control evaluation.
    S0 = basis[0]
    r0 = -(B[0]+B[1]*S0[0]+B[2]*S0[1])/(A[0]+A[1]*S0[0]+A[2]*S0[1])
    vg = (r0-u0)/t
    assert vg.denominator()(0) == 0 and vg.numerator()(0) != 0
    f1g,f2g = [f(vg) for f in [f1,f2]]
    mg = -f1g+E.a1()*f2g/2
    sg = ((2*(S0[0]+E.b2()/12)+cx)*f2g*f2g-mg*mg)/h**3
    assert sg*sg == F(vg)
    q = [R(c) for c in F.list()]
    at0 = QQ(q[0](0))
    assert all(f(0) == 0 for f in q[1:]) and at0 and at0.is_square()
    rows = []
    for case in protocol['cases']:
        tau = QQ(case['parameter'])
        if tau == 0:
            row = {**case, 'status': 'CALIBRATED_CONSTANT_SQUARE_AFFINE_FIBRE',
                   'constant_square': str(at0), 'generic_rule_has_pole': True,
                   'restriction_rank': 1, 'quartic_invariants': ['0','0'],
                   'elliptic_discriminant_nonzero': bool(E.discriminant()(tau))}
        else:
            v0,s0 = vg(tau),sg(tau)
            ft = PolynomialRing(QQ,'v')([f(tau) for f in q])
            assert ft(v0) == s0*s0 and ft.discriminant()
            row = {**case, 'status': 'RATIONAL_PENCIL_POINT_IS_GENERIC_SECTION0',
                   'v': str(v0), 's': str(s0), 'restriction_rank': 2,
                   'elliptic_basis_word': [1]+[0]*16,
                   'elliptic_point': [str(S0[j](tau)) for j in range(2)],
                   'quartic_invariants': [str(I(tau)),str(J(tau))]}
        rows.append(row)
        retain(DIR/f'case-{case["index"]:02d}.json',row)
    result = {'classification': 'verified application and new deduction; independent replay required',
        'status': 'PASS_CANDIDATE_PENCIL_CALIBRATION_AND_GENERIC_CONTROL_WITNESSES',
        'kappa': str(kappa), 'u0': str(u0), 'h': list(map(str,h.list())),
        'quartic_coefficients_in_v': [rec(c) for c in F.list()],
        'generic_v_rule': rec(vg), 'generic_s_rule': rec(sg),
        'generic_v_map_degree': int(max(vg.numerator().degree(),vg.denominator().degree())),
        'invariant_identities': {'I': '(kappa*t)^4*c4(E)', 'J': '2*(kappa*t)^6*c6(E)',
                                 'quartic_discriminant': '256*(kappa*t)^12*Delta(E)',
                                 'Jacobian_short_model_scale': '6*kappa*t'},
        'cases': rows,
        'boundary': 'The fixed-cover exclusions do not extend to the pencil. At nonzero addresses its rational points parameterize the original elliptic fibre, not an independent arithmetic incidence condition. At302 the affine chart has a constant square because the pencil was fitted there; the elliptic fibre itself has good reduction.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths+[DIR/'protocol.json',Path(__file__)]}}
    retain(DIR/'construction.json',result)
    print(result['status'],'degree',result['generic_v_map_degree'],flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    DIR.mkdir(parents=True,exist_ok=True)
    build()
