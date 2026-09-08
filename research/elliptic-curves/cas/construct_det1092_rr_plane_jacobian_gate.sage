#!/usr/bin/env sage-python
"""RR double-plane equation and one fixed member's Jacobian proof gate.

Exactly the source-selected B member, at most18 fixed primes through197.
Finite-field counts only; no rational point search or Selmer/class groups.
"""
import hashlib
import json
import runpy
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, lcm, gcd

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
CAS = ROOT/'elliptic-curves/cas'
OUT = ART/'det1092_rr_plane_jacobian_gate_v1.json'
CHECKPOINT = ROOT/'artifacts/local/elliptic-curves/det1092-rr-plane-jacobian-gate-v1'
PRIMES = [17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def emit_new(p, data):
    if p.exists():
        raise FileExistsError(p)
    p.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')


def run():
    if OUT.exists() or CHECKPOINT.exists():
        raise FileExistsError('Preserve the existing bounded experiment')
    pp = ART/'curve302_recovered_mw17_parent_v1.json'
    np = ART/'det1092_first_centre_rr_net_v1.json'
    loader = CAS/'load_curve302_recovered_parent.sage'
    CHECKPOINT.mkdir(parents=True)
    protocol = {'classification': 'verified application; historically calibrated centre',
                'fixed_member': 'B, i.e. u=v=0; no member substitutions',
                'primes': PRIMES,
                'stop': 'First ordinary absolutely simple Jacobian reduction, or the fixed prime list exhausted. Retain every exclusion and count. No extensions to this budget.',
                'inputs': {str(p.relative_to(ROOT)):sha(p) for p in [pp,np,loader,Path(__file__)]},
                'limits': {'members':1, 'prime_cap':18, 'prime_bound':197,
                           'point_searches':0, 'Selmer_runs':0, 'class_group_runs':0}}
    emit_new(CHECKPOINT/'protocol.json', protocol)
    net = json.loads(np.read_text())
    E, basis, unused_t0 = runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    R = E.base_ring().ring()
    centre = -sum((ZZ(n)*P for n,P in zip(net['trace_word'],basis)), E(0))
    cx = centre[0]+E.b2()/12
    cy = centre[1]+(E.a1()*centre[0]+E.a3())/2
    aa = -E.c4()/48
    h = R(centre[0].denominator().sqrt())
    nx, ny = R(cx*h*h), R(cy*h**3)
    B0 = PolynomialRing(QQ, 'r')
    r = B0.gen()
    S = PolynomialRing(B0, 'T')
    T = S.gen()
    A, B = [list(map(R, net[k])) for k in ['A','B']]
    f1, f2 = S(B[1])+r*S(A[1]), S(B[2])+r*S(A[2])
    mnum = -f1+S(R(E.a1()))*f2/2
    L, rem = f2.quo_rem(S(h))
    assert not rem
    raw = mnum**4-6*S(nx)*mnum*mnum*L*L-8*S(ny)*mnum*L**3-(3*S(nx)**2+4*S(R(aa))*S(h)**4)*L**4
    q, rem = raw.quo_rem(S(h)**6)
    assert not rem
    vals = [QQ(c) for row in q.list() for c in row.list()]
    den = lcm(c.denominator() for c in vals)
    content = gcd(ZZ(c*den) for c in vals)
    q *= QQ(den)/content
    scale = QQ(content)/den
    plane_ring = PolynomialRing(QQ, ['X','Y','Z'])
    X, Y, Z = plane_ring.gens()
    terms = []
    plane = plane_ring(0)
    for i, row in enumerate(q.list()):
        for j, c in enumerate(row.list()):
            if not c:
                continue
            assert i+j <= 6 and j <= 4
            plane += c*X**i*Y**j*Z**(6-i-j)
            terms.append({'X':i,'Y':j,'Z':6-i-j,'coefficient':str(c)})
    assert plane.is_homogeneous() and plane.degree() == 6
    # At N=[0:1:0], the tangent cone comes from the Y^4 coefficient.
    tangent = sum((c*X**i*Z**k for (i,j,k),c in plane.dict().items() if i+k == 2), plane_ring(0))
    assert tangent and all(i+k >= 2 for i,j,k in plane.dict())
    c20, c11, c02 = [tangent.monomial_coefficient(v) for v in [X*X,X*Z,Z*Z]]
    tangent_disc = c11*c11-4*c20*c02
    assert tangent_disc
    fixed = R([row(0) for row in q.list()])*scale
    assert fixed.degree() == 6 and fixed.gcd(fixed.derivative()).degree() == 0
    trials = []
    proof = None
    for p in PRIMES:
        if any(c.denominator() % p == 0 for c in fixed.list()):
            trial = {'prime':p, 'status':'EXCLUDED_NONINTEGRAL_MODEL'}
        else:
            F = GF(p)
            Rp = PolynomialRing(F, 'x')
            f = Rp(fixed.list())
            if f.degree() != 6 or not f.discriminant():
                trial = {'prime':p, 'status':'EXCLUDED_NONSmooth_DEGREE6_MODEL',
                         'coefficients':list(map(int,f.list()))}
            else:
                counts = []
                for degree in [1,2]:
                    field = F if degree == 1 else GF(p*p, 'a')
                    fK = f.change_ring(field)
                    count = sum(1 if fK(x)==0 else 2 if fK(x).is_square() else 0 for x in field)
                    count += 2 if field(f[6]).is_square() else 0
                    counts.append(int(count))
                s1, s2 = p+1-counts[0], p*p+1-counts[1]
                e2 = ZZ(s1*s1-s2)/2
                assert e2 in ZZ
                e2 = ZZ(e2)
                U = PolynomialRing(QQ, 'x')
                x = U.gen()
                frob = x**4-s1*x**3+e2*x*x-p*s1*x+p*p
                irreducible = frob.is_irreducible()
                ordinary = bool(e2 % p)
                exceptions = [s1==0, s1*s1==p+e2, s1*s1==2*e2, s1*s1==3*e2-3*p]
                absolute = ordinary and irreducible and not any(exceptions)
                trial = {'prime':p, 'status':'ABSOLUTELY_SIMPLE_REDUCTION' if absolute else 'NO_ABSOLUTE_SIMPLICITY_CERTIFICATE_AT_THIS_PRIME',
                         'coefficients':list(map(int,f.list())), 'point_counts_Fp_Fp2':counts,
                         'frobenius_coefficients':list(map(int,frob.list())),
                         'ordinary':ordinary, 'irreducible':irreducible,
                         'Howe_Zhu_exception_equalities':exceptions}
                if absolute:
                    proof = trial
        trials.append(trial)
        emit_new(CHECKPOINT/('%02d.json' % len(trials)), trial)
        print('prime',p,trial['status'], flush=True)
        if proof:
            break
    data = {'classification':'verified application and new deduction',
            'status':'PASS_EXPLICIT_RR_DOUBLE_PLANE_AND_FIXED_JACOBIAN_GATE',
            'protocol':protocol, 'protocol_sha256':sha(CHECKPOINT/'protocol.json'),
            'plane':{'equation':'W^2=scale*B(X,Y,Z) in weights(1,1,1,3)',
                     'scale':str(scale), 'sparse_B':terms,
                     'node':[0,1,0], 'tangent_cone_coefficients_X2_XZ_Z2':list(map(str,[c20,c11,c02])),
                     'tangent_discriminant':str(tangent_disc),
                     'original_fibre_lines':'X=tau*Z, through N',
                     'RR_genus2_lines':'Y=u*Z+v*X, not through N',
                     'generic_surface_map':'r=-B(t,x,y)/A(t,x,y); [X:Y:Z]=[t:r:1]'},
            'fixed_member':{'parameters':[0,0], 'q_coefficients':list(map(str,fixed.list())),
                            'Jacobian_absolute_simplicity':'PROVED' if proof else 'UNKNOWN'},
            'trials':trials, 'absolute_simplicity_witness_prime':None if proof is None else proof['prime'],
            'boundary':'This is the one fixed generic RR member, not a genus2 curve attached to302. Absolute simplicity, if proved, excludes elliptic maps from that member and the geometric generic RR member, not all special RR loci. No Selmer comparison or incidence theorem follows.'}
    emit_new(OUT, data)
    print(data['status'], 'Jacobian',data['fixed_member']['Jacobian_absolute_simplicity'], flush=True)


if __name__ == '__main__':
    run()
