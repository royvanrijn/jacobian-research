"""Finite exact transcriptions of the written all-ring/all-index proofs.

No general integration or actual-image kernel is computed. The universal
statements are proved in the two canonical notes, not extrapolated here.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import sympy as s

p = Path(__file__).resolve().parent
t, z, e = s.symbols('t z e')
sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()


def normal(f, nilpotent=False):
    f = s.expand(f)
    return s.expand(f.coeff(e, 0) + e*f.coeff(e, 1)) if nilpotent else f


def heat(f, X, covariance, sign, nilpotent=False):
    out, term, k = s.S.Zero, normal(f, nilpotent), 0
    while term != 0:
        out += sign**k*term/(2**k*s.factorial(k))
        term = normal(sum(covariance[i, j]*s.diff(term, X[i], X[j])
                          for i in range(len(X)) for j in range(len(X))), nilpotent)
        k += 1
    return normal(out, nilpotent)


def nabla(f, i, X, H, P, nilpotent=False):
    return normal(s.diff(f, X[i]) - (H*s.Matrix(X))[i]*f
                  + t*s.diff(P, X[i])*f, nilpotent)


def cancel_t(f, g, X, H, P, nilpotent=False):
    """Return original polynomial IBP witnesses; reject invalid inputs."""
    n, C = len(X), H.inv()
    N = lambda v: normal(v, nilpotent)
    if N(t*f-sum(nabla(g[i], i, X, H, P, nilpotent) for i in range(n))):
        raise ValueError('not an original t-divisibility identity')
    w = [heat(v, X, C, 1, nilpotent)
         for v in H*s.Matrix([v.subs(t, 0) for v in g])]
    if N(sum(X[i]*w[i] for i in range(n))):
        raise ValueError('not a first Koszul syzygy')
    h = s.zeros(n)
    homogeneous = []
    for v in w:
        pieces = {}
        for exponents, c in s.Poly(v, *X).terms():
            d = sum(exponents)
            pieces[d] = pieces.get(d, 0) + c*s.prod(x**a for x, a in zip(X, exponents))
        homogeneous.append(pieces)
    for i in range(n):
        for j in range(n):
            h[i, j] = heat(sum((s.diff(homogeneous[j].get(d, 0), X[i])
                               - s.diff(homogeneous[i].get(d, 0), X[j]))/(d+1)
                              for d in set(homogeneous[i]) | set(homogeneous[j])),
                              X, C, -1, nilpotent)
    h = C*h*C.T
    q = []
    for i in range(n):
        v = N(g[i]-sum(nabla(h[i, j], j, X, H, P, nilpotent) for j in range(n)))
        if N(v.subs(t, 0)):
            raise ValueError('syzygy contraction did not cancel the constant term')
        q.append(N(s.cancel(v/t)))
        if s.denom(q[-1]) != 1:
            raise ValueError('nonpolynomial time division')
    if N(f-sum(nabla(q[i], i, X, H, P, nilpotent) for i in range(n))):
        raise ValueError('incorrect cancelled original identity')
    return q, h


def multiply(P, Q):
    """Exact PBW order t^a D^i; no external algebra backend."""
    out = {}
    for (a, i), c in P.items():
        for (b, j), d in Q.items():
            for k in range(min(i, b)+1):
                v = math.comb(i, k)*s.ff(b, k)*c*d
                key = a+b-k, i+j-k
                out[key] = out.get(key, 0)+v
    return {k: s.expand(v) for k, v in out.items() if s.expand(v) != 0}


def add(P, Q):
    out = dict(P)
    for key, v in Q.items():
        out[key] = s.expand(out.get(key, 0)+v)
    return {k: v for k, v in out.items() if v != 0}


def scale(P, c):
    return {k: s.expand(v*c) for k, v in P.items() if s.expand(v*c) != 0}


def euler_A(b):
    power, op = {(0, 0): s.S.One}, {}
    for j in range(s.degree(b, z)+1):
        op = add(op, scale(power, s.expand(b).coeff(z, j)))
        power = multiply(power, {(1, 1): s.S.One})
    return add({(0, 0): s.S.One}, scale(multiply({(1, 0): s.S.One}, op), -1))


def left_point(P, j):
    v = s.Symbol('v')
    return s.expand(sum(c*(-1)**a*s.diff(v**(i+j), v, a) for (a, i), c in P.items()))


def right_point(P, j):
    v = s.Symbol('v')
    return s.expand(sum(c*v**i for (a, i), c in multiply({(0, j): s.S.One}, P).items() if a == 0))


def packet_model(b, poles):
    b = s.Poly(b, z, domain=s.QQ)
    if b.degree() < 1 or b.is_zero:
        raise ValueError('positive degree required')
    poles = tuple(s.Rational(x) for x in poles)
    if any(not x for x in poles) or len(set(poles)) != len(poles):
        raise ValueError('distinct nonzero poles required')
    roots = s.polys.polytools.ground_roots(b)
    terminating = any(x.is_Integer and x >= 0 for x in roots)
    if terminating:
        raise ValueError('packet requires nonterminating recurrence')
    r, count = b.degree(), len(poles)
    return dict(b=str(b.as_expr()), poles=[str(x) for x in poles],
                rank=1+r*count, length=1+count,
                characteristic_cycle=[1+r*count, (r+1)*count],
                boundary_delta0=1, quotient_global_semisimple=True,
                origin_minimal=True, unit_essential=True,
                relations=['D e0=0']+[f'(1-({x})*t*({b.as_expr()})(theta)) e{i}=e0'
                                       for i, x in enumerate(poles, 1)])


def verify_model(model):
    if model != packet_model(s.sympify(model['b'], locals={'z': z}), model['poles']):
        raise ValueError('damaged exact packet presentation')
    return True


def main():
    X = s.symbols('x y')
    x, y = X
    cancellation = []
    for H, nil in [(s.eye(2), False), (s.Matrix([[0, 1], [1, 0]]), False),
                   (s.Matrix([[2, 1], [1, 3]]), True)]:
        P = x*y**2+(e*x**2 if nil else x**3)
        h = s.Matrix([[0, x*x+y+e*x*y], [-x*x-y-e*x*y, 0]])
        original = [x+t*y, x*y+e*y]
        f = normal(sum(nabla(original[i], i, X, H, P, nil) for i in range(2)), nil)
        g = [normal(t*original[i]+sum(nabla(h[i, j], j, X, H, P, nil) for j in range(2)), nil)
             for i in range(2)]
        q, contracted = cancel_t(f, g, X, H, P, nil)
        cancellation.append(dict(precision=[[str(v) for v in row] for row in H.tolist()],
                                 coefficient_ring='QQ[e]/e^2' if nil else 'QQ[e]',
                                 f=str(f), g=[str(v) for v in g], cancelled=[str(v) for v in q],
                                 skew_contraction=[[str(v) for v in row] for row in contracted.tolist()]))
    # Finite original Wick jets compute localized principal parts, including nilpotents.
    pp = []
    for k in [1, 2, 4]:
        P, amplitude = e*x+x*x, 1+t*x
        series = 0
        for m in range(k+1):
            series += t**m*heat(amplitude*P**m, (x,), s.eye(1), 1, True).subs(x, 0)/s.factorial(m)
        series = normal(series, True)
        principal = s.expand(sum(series.coeff(t, j)*t**(j-k) for j in range(k)))
        derivative = sum(series.coeff(t, k-j)*(-1)**(j-1)*s.diff(1/t, t, j-1)/s.factorial(j-1)
                         for j in range(1, k+1))
        assert s.expand(principal-derivative) == 0
        pp.append(dict(denominator_power=k, amplitude=str(amplitude), original_phase=str(P),
                       jets=[str(series.coeff(t, j)) for j in range(k)], principal_part=str(principal)))
    # Genuine Gaussian source point kernel is at t=1, not the normalization origin.
    assert s.expand(nabla(1, 0, (x,), s.eye(1), x*x/2)-(t-1)*x) == 0
    assert all(heat(x**(2*j+1), (x,), s.eye(1), 1).subs(x, 0) == 0 for j in range(5))
    source_torsion = dict(phase='x^2/2', nonzero_class='[x]', annihilator='t-1',
                          derivative='[x^3/2]', actual_period='0',
                          nonzero_proof='Highest x-degree in nabla(g)=x forces g=1/(t-1), impossible in QQ[t,x].')
    point_tests, terminating = [], []
    v = s.Symbol('v')
    bs = [z+1, (2*z+1)*(2*z+3), 3*(3*z+1)*(3*z+2), z*(z+1), (z-2)*(z+1)]
    for b in bs:
        A = euler_A(b)
        L = multiply({(0, 1): s.S.One}, A)
        for j in range(6):
            lower = v**(j-1) if j else 0
            assert s.expand(left_point(A, j)-v**j-j*b.subs(z, -j-1)*lower) == 0
            assert s.expand(left_point(L, j)-v**(j+1)-j*b.subs(z, -j-1)*v**j) == 0
            assert s.expand(right_point(A, j)-v**j+j*b.subs(z, j-1)*lower) == 0
            assert s.expand(right_point(L, j)-v**(j+1)+(j+1)*b.subs(z, j)*v**j) == 0
        coefficients = [s.S.One]
        for j in range(8):
            coefficients.append(coefficients[-1]*b.subs(z, j))
        h = sum(c*t**j for j, c in enumerate(coefficients))
        applied = sum(c*t**a*s.diff(h, t, i) for (a, i), c in A.items())
        assert all(s.expand(applied-1).coeff(t, j) == 0 for j in range(9))
        roots = [j for j in range(8) if b.subs(z, j) == 0]
        if roots:
            n = min(roots)
            assert s.degree(h, t) == n and s.expand(applied) == 1
            terminating.append(dict(b=str(b), first_root=n, actual_polynomial=str(h),
                                    source_generic_rank=int(s.degree(b, z))+1, actual_rank=1,
                                    source_kernel='S_b', splitting_vector='e-h*[A_b]'))
        point_tests.append(dict(b=str(b), A=[[list(k), str(c)] for k, c in sorted(A.items())],
                                L=[[list(k), str(c)] for k, c in sorted(L.items())], tested_indices=6))
    models = [packet_model(bs[0], [1]), packet_model(bs[1], [1, 2]),
              packet_model(bs[2], [1, -1])]
    assert all(verify_model(m) for m in models)
    alpha, kappa, delta = s.symbols('alpha kappa delta')
    L = {(3, 2): kappa**2, (2, 1): 2*(alpha+1)*kappa**2,
         (1, 1): -kappa*delta, (0, 1): -1,
         (1, 0): alpha*(alpha+1)*kappa**2, (0, 0): -alpha*kappa*delta}
    for m in range(6):
        lower = v**(m-1) if m else 0
        left = -v**(m+1)+kappa*delta*(m+1-alpha)*v**m-kappa**2*m*(m-alpha)*(m+1-alpha)*lower
        right = -v**(m+1)-kappa*delta*(m+alpha)*v**m+kappa**2*m*(m+alpha)*(m+alpha-1)*lower
        assert s.expand(left_point(L, m)-left) == 0
        assert s.expand(right_point(L, m)-right) == 0
    Z, W, T = s.symbols('Z W T')
    glue = []
    for nu in [1, 2, 3]:
        h = (1-Z)**(nu+1)*T
        Q = [(1-Z)**j for j in range(1, nu+1)]
        lifts = []
        for j in range(nu-1):
            R = s.cancel((Q[j]-Q[-1])/Z)
            assert s.denom(R) == 1
            # nabla_W=-Z+t*h in the split Gaussian pair.
            assert s.expand((-Z+t*h)*(-R)-(Q[j]-Q[-1]-t*h*R)) == 0
            lifts.append(str(s.expand(h*R)))
        row = s.ones(1, nu)
        kernel = s.Matrix(nu, nu-1, lambda i, j: int(i == j)-int(i == nu-1))
        assert row.rank() == 1 and row*kernel == s.zeros(1, nu-1)
        assert kernel.rank() == nu-1
        glue.append(dict(nu=nu, rank=2*nu, boundary_dimension=nu-1,
                         principal_part_row=[1]*nu, kernel_columns=[[int(v) for v in row] for row in kernel.tolist()],
                         original_polynomial_boundary_markers=lifts,
                         origin_minimal=(nu == 1), unit_essential=(nu == 1),
                         indecomposable=True))
    powers = []
    for k, poles in [(2, [1]), (3, [1]), (5, [1]), (2, [1, -1]), (3, [1, -1])]:
        b = s.expand(k*s.prod(k*z+j for j in range(1, k)))
        folded = sorted(set(s.Rational(lam)**k for lam in poles))
        model = packet_model(b, folded)
        assert verify_model(model)
        for m in range(6):
            a = s.factorial(k*m)/s.factorial(m)
            assert s.factorial(k*(m+1))/s.factorial(m+1) == a*b.subs(z, m)
            hyper = s.rf(1, m)*s.prod(s.rf(s.Rational(j, k), m) for j in range(1, k))*k**(k*m)/s.factorial(m)
            assert hyper == a
        vandermonde = s.Matrix([[s.Rational(lam)**j for lam in poles] for j in range(len(poles))])
        assert vandermonde.det() != 0
        powers.append(dict(power=k, original_poles=poles, folded_poles=[str(x) for x in folded], model=model))
    # A selected even subsequence can vanish although its derivative-marker block survives.
    assert all(1-(-1)**(2*m) == 0 for m in range(5))
    assert all((2*m+1)*(1-(-1)**(2*m+1)) != 0 for m in range(5))
    rejected = []
    def reject(name, fn):
        try:
            fn()
        except (ValueError, TypeError, KeyError):
            rejected.append(name)
        else:
            raise AssertionError('accepted mutation: '+name)
    reject('invalid t-divisibility premise', lambda: cancel_t(1, [0, 0], X, s.eye(2), x))
    reject('constant b', lambda: packet_model(s.Integer(1), [1]))
    reject('terminating recurrence substituted', lambda: packet_model(z-2, [1]))
    reject('zero pole', lambda: packet_model(z+1, [0]))
    reject('unfolded duplicate pole', lambda: packet_model(z+1, [1, 1]))
    for field, value in [('rank', 1), ('boundary_delta0', 0), ('quotient_global_semisimple', False),
                         ('unit_essential', False), ('relations', [])]:
        damaged = json.loads(json.dumps(models[1])); damaged[field] = value
        reject('damaged '+field, lambda damaged=damaged: verify_model(damaged))
    certificate = dict(cancellation=cancellation, principal_parts=pp, genuine_away_origin_kernel=source_torsion,
                       Euler_point_tests=point_tests, terminating=terminating, global_models=models,
                       pullback_symbolic_point_indices=6, pullback_glue=glue, power_packets=powers,
                       decimation_zero_rescued_by_original_derivative_marker=True, rejections=rejected,
                       universal_claims_use_written_proofs=True, actual_image_kernel_computed=False)
    out = p/'gaussian_origin_global_lattice_certificate.json'
    out.write_text(json.dumps(certificate, indent=2)+'\n')
    names = ['TIME_ZERO_SATURATION_AND_PRINCIPAL_PART_RECONSTRUCTION.md',
             'GLOBAL_EULER_AND_GAUSSIAN_PULLBACK_LATTICES.md']
    receipt = dict(completed=True, canonical_proofs_sha256={n: sha(p/n) for n in names},
                   checker_sha256=sha(Path(__file__)), certificate_sha256=sha(out),
                   runtime_sha256=sha(p/'gaussian_origin_global_lattice_runtime.json'),
                   cancellation_controls=len(cancellation), principal_part_controls=len(pp),
                   Euler_point_models=len(point_tests), terminating_controls=len(terminating),
                   global_models=len(models), pullback_glue_controls=len(glue), power_packet_controls=len(powers),
                   rejections=len(rejected), historical_cusp_and_cone_restarted=False,
                   independent_replay=False, formal_verification=False, external_review=False,
                   assurance='Finite exact transcriptions only; universal all-ring and all-index results use the written proofs.')
    (p/'gaussian_origin_global_lattice_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: v for k, v in receipt.items() if 'sha256' not in k}))


if __name__ == '__main__':
    main()
