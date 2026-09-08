#!/usr/bin/env sage-python
"""Independent literal elimination, invariant covariance, and marked maps.

One generic polynomial identity, nine fixed addresses, 25-second cap.
No point search, exceptional coordinates, V3 input or constructor import.
"""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_seed_pencil_gate_v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def invariants(f):
    a,b,c,d,e = [f[j] for j in range(5)]
    return (12*e*a-3*d*b+c*c,
            72*e*c*a+9*d*c*b-27*e*b*b-27*d*d*a-2*c**3)


def covariance():
    """Universal identities over QQ, not interpolation of the family data."""
    R = PolynomialRing(QQ,['cx','cy','aa','a','b','c','d'])
    cx,cy,aa,a,b,c,d = R.gens()
    S = PolynomialRing(R,'z')
    z = S.gen()
    n,den = a*z+b,c*z+d
    f = n**4-6*cx*n*n*den*den-8*cy*n*den**3-(3*cx*cx+4*aa)*den**4
    I,J = invariants(f)
    det = a*d-b*c
    assert I == -48*aa*det**4
    assert J == -1728*(cy*cy-cx**3-aa*cx)*det**6
    # Discriminant identity by a universal resultant, independent of the
    # constructor's concrete rational-function discriminant calculation.
    C = PolynomialRing(QQ,['c0','c1','c2','c3','c4'])
    Z = PolynomialRing(C,'z')
    f = Z(C.gens())
    I,J = invariants(f)
    assert f.resultant(f.derivative()) == C.gen(4)*(4*I**3-J**2)/27
    return True


def verify():
    d = json.loads((DIR/'construction.json').read_text())
    protocol = json.loads((DIR/'protocol.json').read_text())
    for collection in [d['inputs'],protocol['inputs']]:
        for path,digest in collection.items():
            assert sha(ROOT/path) == digest
    parent = json.loads((ART/'curve302_recovered_mw17_parent_v1.json').read_text())
    net = json.loads((ART/'det1092_first_centre_rr_net_v1.json').read_text())
    pencil = json.loads((ART/'det1092_first_witness_pencil_genus_gate_v1.json').read_text())
    R = PolynomialRing(QQ,'t')
    t = R.gen()
    K = R.fraction_field()

    def rat(r):
        return K(R(r['numerator']))/R(r['denominator'])

    E = EllipticCurve(K,[rat(r) for r in parent['a_invariants']])
    basis = [E([rat(r) for r in p]) for p in parent['basis_weierstrass_coordinates']]
    Z = -sum((n*P for n,P in zip(net['trace_word'],basis)),E(0))
    h = R(d['h'])
    assert h*h == Z[0].denominator()
    A,B = [[R(row) for row in net[k]] for k in ['A','B']]
    u0 = QQ(d['u0'])
    kappa = QQ(d['kappa'])
    assert u0 == QQ(pencil['u0']) and kappa
    assert A[2]*B[1]-A[1]*B[2] == kappa*h**3
    V = PolynomialRing(K,'v')
    v = V.gen()
    f0,f1,f2 = [V(B[j]+u0*A[j])+t*A[j]*v for j in range(3)]
    assert f0+f1*(-Z)[0]+f2*(-Z)[1] == 0
    X = PolynomialRing(V,'x')
    x = X.gen()
    a1,a2,a3,a4,a6 = E.a_invariants()
    line = f0+f1*x
    pol = line*line-a1*x*line*f2-a3*line*f2-f2*f2*(x**3+a2*x*x+a4*x+a6)
    residual,rem = pol.quo_rem(x-(-Z)[0])
    assert not rem and residual.degree() == 2
    F = V([rat(row) for row in d['quartic_coefficients_in_v']])
    assert residual[1]**2-4*residual[0]*residual[2] == h**6*F
    assert F.degree() == 4 and all(c.denominator().degree() == 0 for c in F)
    q = [R(c) for c in F.list()]
    # Bind to the immutable genus2 pencil, including its exact twist factor.
    scale = QQ(pencil['branch_squareclass_scale'])
    oldq = pencil['q_coefficients_t_then_v']
    for j in range(5):
        assert q[j] == scale*R([QQ(row[j]) if len(row)>j else QQ(0) for row in oldq])
    assert covariance()
    mnum = -f1+a1*f2/2
    det = mnum[1]*f2[0]-mnum[0]*f2[1]
    assert det == kappa*t*h**3
    cx = Z[0]+E.b2()/12
    cy = Z[1]+(a1*Z[0]+a3)/2
    aa,bb = -E.c4()/48,-E.c6()/864
    assert cy*cy-cx**3-aa*cx == bb
    assert (mnum**4-6*cx*mnum*mnum*f2*f2-8*cy*mnum*f2**3
            -(3*cx*cx+4*aa)*f2**4) == h**6*F
    # Covariance implies I=(kappa*t)^4*c4 and J=2*(kappa*t)^6*c6.
    # Hence the Jacobian short model is scaled from E by 6*kappa*t.
    scale_jac = 6*kappa*t
    assert -27*(kappa*t)**4*E.c4() == scale_jac**4*aa
    assert -54*(kappa*t)**6*E.c6() == scale_jac**6*bb
    vg,sg = [rat(d[k]) for k in ['generic_v_rule','generic_s_rule']]
    P = basis[0]
    ag = A[0]+A[1]*P[0]+A[2]*P[1]
    bg = B[0]+B[1]*P[0]+B[2]*P[1]
    assert ag and bg+(u0+t*vg)*ag == 0
    assert vg.numerator()(0) and vg.denominator()(0) == 0
    assert vg.denominator().derivative()(0) != 0  # simple pole at302
    degree = max(vg.numerator().degree(),vg.denominator().degree())
    assert degree == d['generic_v_map_degree'] == 7
    # The discriminant-square identity follows from the exact quadratic root.
    r0,r1,r2 = [r(vg) for r in residual.list()]
    assert r0+r1*P[0]+r2*P[0]**2 == 0
    assert sg == -(2*r2*P[0]+r1)/h**3
    results = []
    for row in d['cases']:
        tau = QQ(row['parameter'])
        assert row == json.loads((DIR/f'case-{row["index"]:02d}.json').read_text())
        assert row['parameter'] == protocol['cases'][row['index']]['parameter']
        restriction = matrix(QQ, [[B[j](tau)+u0*A[j](tau) for j in range(3)],
                                  [tau*A[j](tau) for j in range(3)]])
        assert restriction.rank() == row['restriction_rank']
        if tau == 0:
            assert E.discriminant()(tau) and restriction.rank() == 1
            assert all(c(tau) == 0 for c in q[1:])
            assert q[0](tau) == QQ(row['constant_square']) and q[0](tau).is_square()
            results.append({'index':row['index'],'status':'PASS_CALIBRATION_CHART_COLLAPSE_NOT_BAD_ELLIPTIC_FIBRE'})
            continue
        v0,s0 = map(QQ,[row['v'],row['s']])
        assert v0 == vg(tau) and s0 == sg(tau) and s0
        ft = PolynomialRing(QQ,'v')([c(tau) for c in q])
        assert ft(v0) == s0*s0
        I,J = invariants(ft)
        assert [I,J] == list(map(QQ,row['quartic_invariants']))
        assert I == (kappa*tau)**4*E.c4()(tau)
        assert J == 2*(kappa*tau)**6*E.c6()(tau)
        assert (4*I**3-J**2)/27 == 256*(kappa*tau)**12*E.discriminant()(tau) != 0
        fg = [f2[j](tau) for j in range(2)]
        f2_at = fg[0]+fg[1]*v0
        m = (mnum[0](tau)+mnum[1](tau)*v0)/f2_at
        xx = (m*m-cx(tau)+h(tau)**3*s0/f2_at**2)/2
        yy = m*(xx-cx(tau))-cy(tau)
        literal_x = xx-E.b2()(tau)/12
        literal_y = yy-(a1(tau)*literal_x+a3(tau))/2
        assert [literal_x,literal_y] == [P[0](tau),P[1](tau)]
        assert [literal_x,literal_y] == list(map(QQ,row['elliptic_point']))
        assert row['elliptic_basis_word'] == [1]+[0]*16
        results.append({'index':row['index'],'status':'PASS_GENERIC_RATIONAL_POINT_ON_WHOLE_PENCIL',
                        'elliptic_relative_class':'ZERO', 'smooth_quartic_genus':1,
                        'Jacobian_Q_isomorphic_to_marked_E':True})
    out = {'classification':'verified application and new deduction',
        'status':'PASS_INDEPENDENT_SEED_PENCIL_CALIBRATION_GATE',
        'literal_residual_elimination':True, 'universal_invariant_covariance':True,
        'discriminant_formula':'disc_v(F)=256*kappa^12*t^12*Delta(E)',
        'Jacobian_short_coordinate_scale':'6*kappa*t',
        'generic_control_witness_degree':7, 'generic_rule_simple_pole_at302':True,
        'cases':results,
        'general_anchor_deduction':'For B+(u0+v*(t-tau0))*A the same determinant equals kappa*(t-tau0)*h^3; the factor(t-tau0)^12 is a pencil-chart effect, not an elliptic rank or Selmer discriminator.',
        'boundary':'Nonsplitting of the fixed v0 member cannot be upgraded to nonexistence of rational points on the whole pencil. The eight constructed witnesses are deliberately generic, not new seeds. The single-seed construction remains calibrated and its useful parameter selector/global arithmetic class remains unknown.',
        'limits':protocol['limits'],
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [DIR/'protocol.json',DIR/'construction.json']},
        'checker_sha256':sha(Path(__file__))}
    text = json.dumps(out,indent=2,sort_keys=True)+'\n'
    p = DIR/'replay.json'
    if p.exists():
        assert p.read_text() == text
    else:
        p.write_text(text)
    print(out['status'],flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    verify()
