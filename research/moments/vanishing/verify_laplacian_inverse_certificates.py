"""Small exact inverse certificates; no high-dimensional witness producer."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def run():
    u, w, v, z, t = s.symbols('u w v z t')
    xs = (u, w, v, z)
    metric = s.Matrix([[0, 0, 1, 0], [0, 0, 0, 1],
                       [1, 0, 0, 0], [0, 1, 0, 0]])
    zero = lambda f: s.expand(f) == 0
    grad = lambda f: metric * s.Matrix([s.diff(f, x) for x in xs])
    lap = lambda f: s.expand(2*(s.diff(f, u, v)+s.diff(f, w, z)))

    def iterate(f, m):
        for _ in range(m):
            f = lap(f)
        return f

    def compose(fs, gs):
        return [s.expand(f.subs(dict(zip(xs, gs)), simultaneous=True)) for f in fs]

    rows = []
    contraction_count = 0
    multiplier_count = 0
    for label, f, g in [('homogeneous_cubic', u**2, w**3),
                        ('mixed_cubic', u**2+u, w**3+2*w**2),
                        ('nilpotent_quadratic', u, w**2/2)]:
        p = f*z+g
        potential = f*z+g.subs(w, w+t*f)
        forward = s.Matrix(xs)-t*grad(p)
        inverse = s.Matrix(xs)+t*grad(potential)
        assert all(zero(a-b) for a, b in zip(compose(forward, inverse), xs))
        assert all(zero(a-b) for a, b in zip(compose(inverse, forward), xs))
        a = metric*s.hessian(p, xs)
        assert a**4 == s.zeros(4)
        assert zero(lap(potential))
        assert zero(s.diff(potential, t)-(grad(potential).T*metric*grad(potential))[0]/2)
        at_inverse = compose([p], inverse)[0]
        gp = s.Matrix(compose(grad(p), inverse))
        assert zero(potential-at_inverse+t*(gp.T*metric*gp)[0]/2)
        nt = s.degree(potential, t)
        for m in range(1, 5):
            assert zero(iterate(p**m, m))
            contraction_count += 1
        for k in (1, 2):
            for m in range(0, k*nt+2):
                c = s.expand(potential**k).coeff(t, m)
                rhs = s.factorial(k)*iterate(p**(m+k), m)/(2**m*s.factorial(m)*s.factorial(m+k))
                assert zero(c-rhs)
                contraction_count += 1
        # These are boundary and post-boundary controls, not universal proofs.
        for q in (v, z, v*z):
            d = s.Poly(q, *xs).total_degree()
            m = (nt+1)*d+1
            assert zero(iterate(q*p**m, m))
            multiplier_count += 1
        # Reject removing the last time coefficient of the true inverse.
        max_t = max(s.degree(c, t) for c in inverse if c != 0)
        bad = s.Matrix([s.expand(c)-s.expand(c).coeff(t, max_t)*t**max_t for c in inverse])
        assert any(not zero(a-b) for a, b in zip(compose(forward, bad), xs))
        # A finite homogeneous-x recursion control, separate from a general producer.
        a0 = a.subs(dict.fromkeys(xs, 0))
        ct = sum((t**j*a0**j for j in range(4)), s.zeros(4))
        assert (s.eye(4)-t*a0)*ct == s.eye(4)
        homs = []
        for k in range(1, max(s.Poly(c, *xs).total_degree() for c in inverse)+1):
            part = s.Matrix([sum(coeff*s.prod(x**power for x, power in zip(xs, mon))
                           for mon, coeff in s.Poly(c, *xs).terms() if sum(mon) == k)
                           for c in inverse])
            if part == s.zeros(4, 1):
                continue
            assert all(s.degree(c, t) <= 7*k-4 for c in part if c != 0)
            homs.append(k)
        rows.append({'case': label, 'potential_t_degree': int(nt),
                     'inverse_t_degree': int(max_t), 'both_compositions': True,
                     'truncated_false_inverse_rejected': True,
                     'homogeneous_x_degrees_checked': homs})

    x, y, zz = s.symbols('x y zz')
    vv = 1+2*x*y
    base = s.Matrix([x-3*x**2*y-x**3*zz,
                     3*vv**2*x*zz+(9*vv**2-6*vv-2)*y,
                     vv**3*zz+(6*vv**2+2*vv)*y**2])
    assert zero(base.jacobian((x, y, zz)).det()-1)
    for point in [(0, 1, -8), (-1, 1, 4)]:
        assert list(base.subs(dict(zip((x, y, zz), point)))) == [0, 1, 0]

    # Integer checks of the tree degree estimate include binary extremality.
    tree_checks = 0
    for n in range(2, 6):
        for left in range(1, 5):
            for right in range(1, 5):
                h = lambda k: (2*n-1)*k-n
                assert n+h(left)+h(right) == h(left+right)
                tree_checks += 1
    return {'completed': True, 'python_version': platform.python_version(),
            'sympy_version': s.__version__, 'inverse_certificates': rows,
            'contraction_controls': contraction_count,
            'multiplier_tail_controls': multiplier_count,
            'tree_bound_controls': tree_checks,
            'de_bondt_base_map_determinant': '1',
            'de_bondt_base_map_collision_checked': True,
            'high_dimensional_ancillary_producer_executed': False,
            'reported_quartic_witness_independently_replayed': False,
            'general_inverse_or_closed_locus_producer_implemented': False,
            'all_degree_proof_inferred_from_finite_controls': False,
            'general_Gaussian_uniformity_proved': False,
            'nonreduced_coefficient_uniformity_proved': False,
            'ternary_Laplacian_power_frontier_changed': False,
            'independent_review': False, 'formal_verification': False,
            'novelty_claim': False}


if __name__ == '__main__':
    path = Path(__file__).resolve()
    result = run()
    result['checker_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    destination = path.with_name('laplacian_inverse_verification.json')
    destination.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
