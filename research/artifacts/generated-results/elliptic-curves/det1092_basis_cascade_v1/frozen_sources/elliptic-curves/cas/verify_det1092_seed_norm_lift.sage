#!/usr/bin/env sage-python
"""Independent cubic reduction and factor-free specific-class obstructions.

Checks one calibrated class family. No point search or integer factorization.
The copied seed abscissa is explicitly NOT claimed to be an oracle-free input.
"""
import hashlib
import json
import math
import runpy
import signal
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, prod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_seed_norm_lift_v1'
RANK_CHECKER = ROOT/'elliptic-curves/cas/verify_det1092_single_seed_covers.sage'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify():
    d = json.loads((DIR/'construction.json').read_text())
    protocol = json.loads((DIR/'protocol.json').read_text())
    for rows in [d['inputs'],protocol['inputs']]:
        for path,digest in rows.items():
            assert sha(ROOT/path) == digest
    pp = ART/'curve302_recovered_mw17_parent_v1.json'
    parent = json.loads(pp.read_text())
    R = PolynomialRing(QQ,'t')
    K = R.fraction_field()

    def rat(row):
        return K(R(row['numerator']))/R(row['denominator'])

    E = EllipticCurve(K,[rat(r) for r in parent['a_invariants']])
    T = PolynomialRing(K,'theta')
    theta = T.gen()
    cubic = T([rat(r) for r in d['cubic_coefficients']])
    assert cubic == theta**3+5*theta**2+(16*E.a4()+8)*theta+64*E.a6()+16
    k = QQ(d['Xstar'])
    D = rat(d['D'])
    assert D == cubic(k) and str(k) == protocol['Xstar']
    alpha = T([rat(r) for r in d['alpha_coefficients']])
    assert alpha == D*(k-theta)
    norm_matrix = matrix(K,3,3,lambda i,j:((alpha*theta**j)%cubic)[i])
    assert norm_matrix.det() == D**4
    assert rat(d['norm_square_root']) == D**2
    disc = K(cubic.discriminant())
    dn = D.numerator()
    assert dn.degree() == 12 and D.denominator().degree() == 0
    assert dn.gcd(dn.derivative()) == dn.gcd(disc.numerator()) == 1
    # Rebuild (k-theta)*gamma^2 directly in the cubic algebra, independently
    # of the constructor's explicit coefficient formulas.
    Z = PolynomialRing(K,['z0','z1','z2'])
    z0,z1,z2 = Z.gens()
    S = PolynomialRing(Z,'theta')
    th = S.gen()
    f = S(cubic.list())
    gamma = z0+z1*th+z2*th**2
    residue = ((k-th)*gamma**2)%f
    assert residue[2](1,0,0) == 0 and residue[1](1,0,0) == -1
    U = PolynomialRing(K,'r')
    z = [U([rat(c) for c in row]) for row in d['conic_parametrization']]
    rr = [c(*z) for c in residue.list()]
    assert rr[2] == 0
    assert z[1] == U.gen()*z[2]
    a = cubic[2]
    assert [zz(k+a) for zz in z] == [D,0,0]
    Q = U([rat(c) for c in d['quartic_coefficients']])
    R0 = U([rat(c) for c in d['map_R0']])
    N = U([rat(c) for c in d['map_Ngamma']])
    assert Q == -D*rr[1] and R0 == rr[0]
    V = PolynomialRing(U,'theta')
    th = V.gen()
    fu = V(cubic.list())
    g = z[0]+z[1]*th+z[2]*th**2
    assert N == fu.resultant(g)  # independent of multiplication determinant
    assert (((k-th)*g*g)%fu) == V(R0)+V(rr[1])*th
    assert D*N*N == R0**3+cubic[2]*R0**2*(-rr[1])+cubic[1]*R0*(-rr[1])**2+cubic[0]*(-rr[1])**3
    # Cover map after w^2=Q. The X-map has degree4, hence C->E has degree4.
    assert Q.degree() == 4 and R0.gcd(rr[1]).degree() == 0
    assert max(R0.degree(),rr[1].degree()) == 4
    e = EllipticCurve(K,[0,cubic[2],0,cubic[1],cubic[0]])
    assert e.c4() == 16*E.c4() and e.c6() == 64*E.c6()
    c0,c1,c2,c3,c4 = Q.list()
    I = 12*c4*c0-3*c3*c1+c2*c2
    J = 72*c4*c2*c0+9*c3*c2*c1-27*c4*c1*c1-27*c3*c3*c0-2*c2**3
    assert I == D**4*e.c4() and J == 2*D**6*e.c6()
    assert (4*I**3-J**2)/27 == 256*D**12*e.discriminant()
    results = []
    for row in d['cases']:
        idx = row['index']
        assert row == json.loads((DIR/f'case-{idx:02d}.json').read_text())
        assert row['parameter'] == protocol['cases'][idx]['parameter']
        tau = QQ(row['parameter'])
        ft = PolynomialRing(QQ,'X')([v(tau) for v in cubic.list()])
        value,delta = QQ(ft(k)),QQ(ft.discriminant())
        assert value == QQ(row['D_value']) and value and delta
        assert list(map(str,ft.list())) == row['cubic']
        sf = [ZZ(2),k.denominator(),value.denominator(),abs(delta.numerator()),
              delta.denominator(),*[c.denominator() for c in ft.list()]]
        assert sf == list(map(ZZ,row['support_factors']))
        support = ZZ(prod(sf))
        rem = int(abs(value.numerator()))
        for factor in row['gcd_strips']:
            g0 = int(factor)
            assert g0 == math.gcd(rem,int(support)) and g0 > 1
            rem //= g0
        assert rem == int(row['coprime_remainder'])
        assert math.gcd(rem,int(support)) == 1
        root = math.isqrt(rem)
        assert root*root <= rem < (root+1)*(root+1)
        assert (root*root == rem) == row['remainder_square']
        if root*root != rem:
            assert row['status'] == 'NOT_SELMER_GOOD_PRIME_RAMIFICATION_OBSTRUCTION'
            gg,aa,bb = ZZ(rem).xgcd(support)
            assert gg == 1 and aa*rem+bb*support == 1
            results.append({'index':idx,'status':'PASS_EXACT_GOOD_PRIME_LOCAL_OBSTRUCTION',
                'coprime_remainder':str(rem),'floor_sqrt':str(root),
                'sqrt_remainder':str(rem-root*root),
                'bezout_with_support':[str(aa),str(bb)],
                'remainder_bits':rem.bit_length(),
                'rational_cover_points':'NONE',
                'Selmer_membership':False,'Sha_class':'NOT_APPLICABLE: not everywhere locally soluble'})
            continue
        assert row['status'] == 'RATIONAL_LIFT_AT_COPIED_ABSCISSA' and tau == 0
        rp,wp = map(QQ,row['quartic_point'])
        assert rp == k+ft[2]
        qt = PolynomialRing(QQ,'r')([v(tau) for v in Q.list()])
        assert qt(rp) == wp*wp and wp
        xx = value*QQ(R0(rp)(tau))/(wp*wp)
        yy = value*value*QQ(N(rp)(tau))/(wp**3)
        assert [xx,yy] == list(map(QQ,row['cubic_point']))
        assert xx == k and yy*yy == value
        e0 = EllipticCurve(QQ,[0,ft[2],0,ft[1],ft[0]])
        generic = []
        for p in parent['basis_weierstrass_coordinates']:
            x,y = [rat(r)(tau) for r in p]
            generic.append(e0([4*x,8*y+4*x+4]))
        proof = runpy.run_path(str(RANK_CHECKER))['finite_rank'](e0,[*generic,e0([xx,yy])])
        assert proof['rank'] == 18 and proof['generic_rank'] == 17
        assert proof['no_rational_2_torsion_prime'] == 31
        assert proof['separator'] is not None
        results.append({'index':idx,'status':'PASS_RATIONAL_SEED_CLASS_OUTSIDE_GENERIC_SPAN',
            'independence':proof,'Selmer_membership':True,'Sha_image':'ZERO',
            'all_rational_cover_images_outside_MW17_rational_span':True,
            'reason':'Image equals Pseed+2E(Q). The rank17 mod2 injection makes the generic saturation have the same mod2 image; the separating character excludes this entire coset.'})
    out = {'classification':'verified application and new deduction',
        'status':'PASS_INDEPENDENT_NORM_LIFT_AND_EIGHT_NON_SELMER_OBSTRUCTIONS',
        'norm_identity':True,'conic_parametrization_and_quartic_map':True,
        'smooth_cover_genus':1,'degree_to_elliptic_fibre':4,'two_cover_class':'D*(Xstar-theta)',
        'quartic_Jacobian_Q_isomorphic_to_E':True,
        'generic_good_fibre_ramification_degree':12,
        'generic_Q_t_rational_point':'NONE: ramified at good-fibre places',
        'cases':results,
        'calibration_boundary':'Xstar is literally copied from the known first seed. This diagnoses a specific norm-corrected continuation and does not remove the seed oracle or supply a prospective class selector.',
        'limits':protocol['limits'],
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [DIR/'protocol.json',DIR/'construction.json',RANK_CHECKER]},
        'checker_sha256':sha(Path(__file__))}
    text = json.dumps(out,indent=2,sort_keys=True)+'\n'
    path = DIR/'replay.json'
    if path.exists():
        assert path.read_text() == text
    else:
        path.write_text(text)
    print(out['status'],flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    verify()
