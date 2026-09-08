#!/usr/bin/env sage-python
"""Calibrated norm-square continuation of the historical first seed class.

This deliberately retains Xstar from the known seed. It is a diagnostic
lift, NOT an oracle-free construction/discovery. No later points or V3 data.
One class, nine frozen addresses, no factorization or point search, 25 sec.
"""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, gcd, prod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_seed_norm_lift_v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def retain(p,data):
    text = json.dumps(data,indent=2,sort_keys=True)+'\n'
    if p.exists():
        assert p.read_text() == text
    else:
        p.write_text(text)


def build():
    paths = [ART/'curve302_recovered_mw17_parent_v1.json',
             ART/'det1092_first_centre_rr_net_replay_v1.json',
             ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent,first,roster = [json.loads(p.read_text()) for p in paths]
    k = 4*QQ(first['reconstructed_point_literal302'][0])
    protocol = {'classification':'calibrated specific-class diagnostic, not prospective discovery',
        'rule':'D(t)=f_t(Xstar); alpha_t=D(t)*(Xstar-theta_t), norm alpha=D^4.',
        'Xstar':str(k),
        'oracle_boundary':'Xstar is the known historical first seed abscissa in X=4x. It is not claimed removed. Only the continuation and control arithmetic are equation-only.',
        'cases':[{'index':i,'parameter':r['parameter']} for i,r in enumerate(roster['cases'])],
        'limits':{'wall_seconds':25,'class_lifts':1,'existing_addresses':9,
                  'gcd_stripping_steps_per_case':64,'new_addresses':0,
                  'integer_factorizations':0,'class_groups':0,'point_searches':0,
                  'full_Selmer_runs':0,'V3_inputs':0,'pilot_changes':0},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'protocol.json',protocol)
    R = PolynomialRing(QQ,'t')
    K = R.fraction_field()

    def rat(r):
        return K(R(r['numerator']))/R(r['denominator'])

    def rec(f):
        f = K(f)
        return {'numerator':list(map(str,f.numerator().list())),
                'denominator':list(map(str,f.denominator().list()))}

    E = EllipticCurve(K,[rat(r) for r in parent['a_invariants']])
    assert E.a_invariants()[:3] == (K(1),K(1),K(1))
    T = PolynomialRing(K,'theta')
    theta = T.gen()
    cubic = theta**3+5*theta**2+(16*E.a4()+8)*theta+64*E.a6()+16
    D = K(cubic(k))
    assert D and D.numerator().degree() == 12
    assert D.numerator().gcd(D.numerator().derivative()) == 1
    disc = K(cubic.discriminant())
    assert D.numerator().gcd(disc.numerator()) == 1
    assert D.denominator().degree() == 0
    a,b,c = cubic[2],cubic[1],cubic[0]
    # Parametrize the conic [(k-theta)*gamma^2]_theta^2=0 through (1,0,0).
    U = PolynomialRing(K,'r')
    r = U.gen()
    z0 = (k+a)*r*r+(2*b-2*a*(k+a))*r+(k+a)*(a*a-b)-a*b+c
    z1 = 2*r*(r-k-a)
    z2 = 2*(r-k-a)
    A0 = z0*z0-2*c*z1*z2+a*c*z2*z2
    A1 = 2*z0*z1-2*b*z1*z2+(a*b-c)*z2*z2
    A2 = z1*z1+2*z0*z2-2*a*z1*z2+(a*a-b)*z2*z2
    R0,R1,R2 = k*A0+c*A2,k*A1-A0+b*A2,(k+a)*A2-A1
    assert R2 == 0
    quartic = -D*R1
    assert quartic.degree() == 4
    # A multiplication determinant records Norm(gamma), not field factorization.
    V = PolynomialRing(U,'theta')
    th = V.gen()
    fc = th**3+U(a)*th**2+U(b)*th+U(c)
    gamma = V(z0)+V(z1)*th+V(z2)*th**2
    mult = matrix(U,3,3,lambda i,j:((gamma*th**j)%fc)[i])
    ngamma = mult.det()
    assert D*ngamma**2 == R0**3+a*R0**2*(-R1)+b*R0*(-R1)**2+c*(-R1)**3
    # quartic model w^2=-D*R1: X=D*R0/w^2; Y=D^2*N(gamma)/w^3.
    rows = []
    for case in protocol['cases']:
        tau = QQ(case['parameter'])
        f = PolynomialRing(QQ,'X')([v(tau) for v in cubic.list()])
        value = QQ(f(k))
        delta = QQ(f.discriminant())
        assert value and delta
        support_factors = [ZZ(2), k.denominator(), value.denominator(),
                           abs(delta.numerator()), delta.denominator(),
                           *[v.denominator() for v in f.list()]]
        support = prod(support_factors)
        rem = abs(value.numerator())
        strips = []
        for step in range(64):
            g = gcd(rem,support)
            if g == 1:
                break
            strips.append(str(g))
            rem //= g
        complete = gcd(rem,support) == 1
        obstruction = bool(complete and not rem.is_square())
        row = {**case,'D_value':str(value),'cubic':list(map(str,f.list())),
               'support_factors':list(map(str,support_factors)),
               'gcd_strips':strips,'coprime_remainder':str(rem),
               'remainder_bits':int(rem.nbits()), 'support_separation_complete':bool(complete),
               'remainder_square':bool(rem.is_square()),
               'local_obstruction_without_factoring':obstruction}
        if obstruction:
            row['status'] = 'NOT_SELMER_GOOD_PRIME_RAMIFICATION_OBSTRUCTION'
        elif value.is_square():
            root = value.sqrt()
            rp = k+f[2]
            wp = value*root
            qp = PolynomialRing(QQ,'r')([v(tau) for v in quartic.list()])
            assert qp(rp) == wp*wp
            row.update({'status':'RATIONAL_LIFT_AT_COPIED_ABSCISSA',
                        'quartic_point':[str(rp),str(wp)],
                        'cubic_point':[str(k),str(root)]})
        else:
            row['status'] = 'UNRESOLVED_BY_FROZEN_RAMIFICATION_TEST'
        rows.append(row)
        retain(DIR/f'case-{case["index"]:02d}.json',row)
    data = {'classification':'calibrated verified application; independent replay required',
        'status':'PASS_CANDIDATE_SPECIFIC_KUMMER_LIFT_AND_RAMIFICATION_GATE',
        'cubic_coefficients':[rec(v) for v in cubic.list()],
        'Xstar':str(k),'D':rec(D),'alpha_coefficients':[rec(k*D),rec(-D),rec(0)],
        'norm_square_root':rec(D**2),
        'conic_parametrization':[[rec(v) for v in z.list()] for z in [z0,z1,z2]],
        'quartic_coefficients':[rec(v) for v in quartic.list()],
        'map_R0':[rec(v) for v in R0.list()],
        'map_Ngamma':[rec(v) for v in ngamma.list()],
        'map':'X=D*R0/w^2, Y=D^2*Ngamma/w^3; r is the quartic coordinate.',
        'generic_ramification':{'D_degree':12,'D_squarefree':True,'D_coprime_cubic_discriminant':True,
            'conclusion':'Ramified at twelve geometric good-fibre places; not a generic everywhere-locally-soluble2-cover class.'},
        'cases':rows,
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'protocol.json',Path(__file__)]}}
    retain(DIR/'construction.json',data)
    print(data['status'],[(r['index'],r['status'],r['remainder_bits']) for r in rows],flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    DIR.mkdir(parents=True,exist_ok=True)
    build()
