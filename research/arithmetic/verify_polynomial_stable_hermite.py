#!/usr/bin/env python3
"""Exact CRT/norm certificates and bounded complete stable-passport controls."""
from collections import Counter, deque
from fractions import Fraction as F
from functools import reduce
from math import factorial, gcd, comb
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_phase_interpolation as j

if not __debug__:
    raise RuntimeError("Stable Hermite checker requires active assertions")
HERE = Path(__file__).resolve().parent
c, m, l, b = j.c, j.m, j.l, j.b

# Each profile lists every nontrivial cycle of one non-transposition factor.
CASES = (
    (((3,), (3,)), 1),
    (((2, 2), (3,)), 1),
    (((2, 2), (2, 2)), 2),
    (((3, 2), (3,)), 1),
    (((3, 3), (4,)), 1),
    (((3,), (3,), (3,)), 1),
    (((2, 2), (2, 2), (3,)), 3),
    (((2, 2), (3,)), 2),
    (((3, 3),), 1),
    (((2, 4),), 1),
    (((3, 2),), 1),
    (((3, 3, 3),), 2),
    (((2, 2, 2),), 2),
    (((3, 3),), 2),
    (((4, 4),), 1),
    (((4, 2),), 2),
    (((5, 2), (4, 3)), 2),
    (((2, 2), (2, 2), (2, 2)), 4),
    (((2, 2), (3,), (3,)), 1),
    (((2, 2), (2, 2), (3,)), 2),
    (((2, 2), (2, 2), (2, 2)), 3),
    (((2, 2), (3,), (3,), (3,)), 1),
    (((2, 2), (2, 2), (3,), (3,)), 2),
    (((2, 2), (2, 2), (2, 2), (3,)), 3),
    (((2, 2), (2, 2), (2, 2), (2, 2)), 4),
    (((2, 2), (3,), (3,), (3,), (3,)), 1),
    (((2, 2), (2, 2), (2, 2), (2, 2), (2, 2)), 5),
    (((2, 2), (3,), (3,), (3,), (3,), (3,)), 1),
    (((3,), (3,), (3,), (3,), (3,)), 1),
    (((3, 2), (4,), (3,), (5,)), 1),
    (((3, 3), (4,), (3,)), 1),
)


def domain(profiles, s, stable=True):
    assert profiles and s >= 1
    assert all(q and all(isinstance(x, int) and x >= 2 for x in q)
               and sum(x-1 for x in q) >= 2 for q in profiles)
    p = sum(map(len, profiles))
    I = sum(x-1 for q in profiles for x in q)
    N = I+s+1
    assert all(sum(q) <= N for q in profiles)
    if stable:
        assert s >= p-len(profiles)
    passport = [tuple(sorted((1,)*(N-sum(q))+tuple(q))) for q in profiles]
    passport += [(1,)*(N-2)+(2,)]*s
    return N, p, I, passport


def hermite_basis(profiles, points):
    assert len(points) == sum(map(len, profiles))
    assert len(set(points)) == len(points)
    units, colors, powers = [], [], []
    at = 0
    for color, profile in enumerate(profiles):
        for q in profile:
            units.append(j.power([-points[at], 1], q))
            colors.append(color)
            powers.append(q)
            at += 1
    Q = [F(1)]
    for U in units:
        Q = j.mul(Q, U)
    E = [[F(0)] for _ in profiles]
    for U, color in zip(units, colors):
        V, rest = j.divmod_poly(Q, U)
        assert rest == [0]
        indicator = j.mod(j.mul(V, j.inverse_mod(V, U)), Q)
        E[color] = j.add(E[color], indicator)
    assert reduce(j.add, E, [F(0)]) == [1]
    for i, Ei in enumerate(E):
        for color, U in zip(colors, units):
            assert j.mod(Ei, U) == [F(i == color)]
    return Q, E, powers, units, colors


def linear_solve(matrix, rhs):
    """Exact affine solution and kernel basis, or reject inconsistency."""
    n = len(matrix[0])
    A = [list(map(F, row))+[F(value)] for row, value in zip(matrix, rhs)]
    pivots = []
    for col in range(n):
        pivot = next((h for h in range(len(pivots), len(A)) if A[h][col]), None)
        if pivot is None:
            continue
        at = len(pivots)
        A[at], A[pivot] = A[pivot], A[at]
        divisor = A[at][col]
        A[at] = [x/divisor for x in A[at]]
        for h in range(len(A)):
            if h != at:
                value = A[h][col]
                A[h] = [x-value*y for x, y in zip(A[h], A[at])]
        pivots.append(col)
    assert all(any(row[:n]) or row[-1] == 0 for row in A)
    particular = [F(0)]*n
    for h, col in enumerate(pivots):
        particular[col] = A[h][-1]
    kernel = []
    for free in range(n):
        if free in pivots:
            continue
        v = [F(0)]*n
        v[free] = 1
        for h, col in enumerate(pivots):
            v[col] = -A[h][free]
        kernel.append(v)
    return particular, kernel, pivots


def matrix_rank(matrix):
    return len(linear_solve(matrix, [0]*len(matrix))[2])


def moment(numerator, denominator, ell):
    """Coefficient of z^-1 in z^ell*numerator/denominator at infinity."""
    M = len(denominator)-1
    inv = [F(1)]
    for at in range(1, ell+1):
        inv.append(-sum(denominator[M-h]*inv[at-h]
                        for h in range(1, at+1) if h <= M))
    return sum(numerator[M-1-ell+h]*inv[h] for h in range(ell+1)
               if 0 <= M-1-ell+h < len(numerator))


def moment_matrix(E, Q, h):
    return [[moment(Ei, Q, ell) for Ei in E] for ell in range(h+1)]


def chart(profiles, s, points, values, residual=None):
    N, p, I, _ = domain(profiles, s)
    Q, E, powers, units, colors = hermite_basis(profiles, points)
    M = len(Q)-1
    values = list(map(F, values))
    assert len(values) == len(profiles)
    d = N-M
    if d < 0:
        h = M-N-1
        leading = moment_matrix(E, Q, h)
        rhs = [F(0)]*h+[F(1)]
        particular, kernel, pivots = linear_solve(leading, rhs)
        assert len(pivots) == h+1
        free = [at for at in range(len(E)) if at not in pivots]
        values = [particular[at]+sum(values[t]*v[at] for t, v in zip(free, kernel))
                  for at in range(len(E))]
        assert [sum(x*v for x, v in zip(row, values)) for row in leading] == rhs
        assert all(sum(row) == 0 for row in leading)
        assert residual is None
    else:
        assert d >= 0 and residual is not None
        residual = j.trim(residual)
        assert len(residual) == d+1 and residual[-1] == 1
        leading = None
    Phi = reduce(j.add, (j.scale(Ei, v) for Ei, v in zip(E, values)), [F(0)])
    f = Phi if d < 0 else j.add(Phi, j.mul(Q, residual))
    assert len(f) == N+1 and f[-1] == 1
    for color, U in zip(colors, units):
        assert j.mod(j.sub(f, [values[color]]), U) == [0]
    Z = [F(1)]
    for u, q in zip(points, powers):
        Z = j.mul(Z, j.power([-u, 1], q-1))
    A, rest = j.divmod_poly(j.derivative(f), j.scale(Z, N))
    assert rest == [0] and len(A) == s+1 and A[-1] == 1
    return f, A, values, leading


def value_norm(A, f):
    """Characteristic polynomial by Newton traces, checked by determinants."""
    d = len(A)-1
    columns = [j.mod([0]*k+f, A) for k in range(d)]
    M = [[columns[h][i] if i < len(columns[h]) else F(0)
          for h in range(d)] for i in range(d)]
    power = [[F(i == h) for h in range(d)] for i in range(d)]
    traces = [F(0)]
    for k in range(1, d+1):
        power = j.matmul(power, M)
        traces.append(sum(power[i][i] for i in range(d)))
    coeff = [F(1)]
    for k in range(1, d+1):
        coeff.append(-sum(coeff[k-h]*traces[h] for h in range(1, k+1))/k)
    C = coeff[::-1]
    # Both sides are monic of degree d. These d+1 exact values certify equality.
    for y in range(d+1):
        assert j.evaluate(C, y) == j.determinant([
            [F(y)*(i == h)-M[i][h] for h in range(d)] for i in range(d)])
    return C


def good(points, values, A, C):
    P = j.from_roots(points)
    return (len(set(values)) == len(values)
            and j.gcd_poly(A, j.derivative(A)) == [1]
            and j.gcd_poly(A, P) == [1]
            and j.gcd_poly(C, j.derivative(C)) == [1]
            and all(j.evaluate(C, v) for v in values))


def residue_at(numerator, denominator, point, order):
    divisor = j.power([-point, 1], order)
    local, rest = j.divmod_poly(denominator, divisor)
    assert rest == [0]
    local = j.compose(local, [point, 1])
    assert local[0]
    inv = [1/local[0]]
    for t in range(1, order):
        inv.append(-sum(local[h]*inv[t-h] for h in range(1, t+1)
                        if h < len(local))/local[0])
    numer = j.compose(numerator, [point, 1])
    return sum(numer[h]*inv[order-1-h] for h in range(order) if h < len(numer))


def residue_matrix(Q, points, powers, colors, k, h):
    return [[sum(residue_at([0]*ell+[1], Q, u, q)
                 for u, q, color in zip(points, powers, colors) if color == i)
             for i in range(k)] for ell in range(h+1)]


def jacobian_control(profiles, s, points, deficient=False):
    """Differentiate moving-contour residues, independently of CRT derivatives."""
    N, p, _, _ = domain(profiles, s, stable=False)
    Q, E, powers, _, colors = hermite_basis(profiles, points)
    h = len(Q)-N-2
    assert h >= 1
    matrix = moment_matrix(E, Q, h)
    assert matrix == residue_matrix(Q, points, powers, colors, len(E), h)
    particular, kernel, pivots = linear_solve(matrix, [0]*h+[1])
    r = len(pivots); delta = h+1-r
    assert bool(delta) == deficient
    # Span of all w with Fw in C*e, with the zero scalar included.
    S = [particular]+kernel
    P = j.from_roots(points)
    Z = [F(1)]
    for u, q in zip(points, powers):
        Z = j.mul(Z, j.power([-u, 1], q-1))
    Bs = []
    derivative_count = 0
    for w in S:
        g = reduce(j.add, (j.scale(Ei, v) for Ei, v in zip(E, w)), [F(0)])
        assert len(g) <= N+1
        Bg, rest = j.divmod_poly(j.derivative(g), Z)
        assert rest == [0] and len(Bg) <= s+1
        Bs.append(Bg)
        for at, (u, q) in enumerate(zip(points, powers)):
            denominator = j.mul(Q, [-u, 1])
            pole_orders = [v+(at == t) for t, v in enumerate(powers)]
            for ell in range(h+1):
                actual = q*sum(w[color]*residue_at([0]*ell+[1], denominator, x, order)
                               for x, order, color in zip(points, pole_orders, colors))
                expected = -u**ell*j.evaluate(Bg, u)/j.evaluate(j.derivative(P), u)
                assert actual == expected, (profiles, ell, at, actual, expected)
                derivative_count += 1
    width = max(map(len, Bs))
    Vdim = matrix_rank([list(Bg)+[0]*(width-len(Bg)) for Bg in Bs])
    assert Vdim == len(E)-r
    _, beta_basis, _ = linear_solve(list(map(list, zip(*matrix))), [0]*len(E))
    assert len(beta_basis) == delta and all(beta[-1] == 0 for beta in beta_basis)
    products = [j.mul(j.trim(beta), Bg) for beta in beta_basis for Bg in Bs]
    if products:
        width = max(map(len, products))
        product_dim = matrix_rank([list(v)+[0]*(width-len(v)) for v in products])
        evaluation_dim = matrix_rank([[j.evaluate(v, u) for u in points] for v in products])
        assert product_dim == evaluation_dim >= delta+Vdim-1
    else:
        product_dim = 0
    # For an actual good monic point this Vandermonde Jacobian is surjective.
    g = reduce(j.add, (j.scale(Ei, v) for Ei, v in zip(E, particular)), [F(0)])
    Bg, _ = j.divmod_poly(j.derivative(g), Z)
    smooth = matrix_rank([[-u**ell*j.evaluate(Bg, u)/j.evaluate(j.derivative(P), u)
                           for u in points] for ell in range(h+1)])
    return {"profiles": profiles, "s": s, "degree": N, "points": list(map(str, points)),
            "h": h, "moment_rank": r, "deficiency": delta,
            "residue_derivative_equalities": derivative_count,
            "derivative_polynomial_space_dimension": Vdim,
            "conormal_product_dimension": product_dim,
            "monicity_jacobian_rank": smooth}


def differential_controls(charts):
    rows = []
    for h in range(1, 5):
        row = next(x for x in charts if x["moment_rows"] == h+1)
        rows.append(jacobian_control(row["profiles"], row["s"],
                                    list(map(F, row["points"]))))
    # The top row vanishes when each color consists of an antipodal pair.
    bad = jacobian_control(((2, 2),)*3, 3, list(map(F, (-3, 3, -2, 2, -1, 1))),
                           deficient=True)
    assert bad["deficiency"] == 1 and bad["conormal_product_dimension"] == 2
    rows.append(bad)
    # At one level below the new bound the same conormal only has dimension one.
    boundary = jacobian_control(((2, 2),)*2, 1, list(map(F, (-2, 2, -1, 1))),
                                deficient=True)
    assert boundary["conormal_product_dimension"] == boundary["deficiency"] == 1
    rows.append(boundary)
    return rows


def chart_controls():
    rows = []
    for profiles, s in CASES:
        N, p, I, _ = domain(profiles, s)
        d = s+1-p
        for trial in range(1, 21):
            points = [F(h*h+trial*h) for h in range(p)]
            values = [F((h+1)*(h+trial+1), 101+2*trial) for h in range(len(profiles))]
            R = None if d < 0 else j.from_roots([-(p+trial+3*h) for h in range(d)])
            try:
                f, A, values, leading = chart(profiles, s, points, values, R)
            except AssertionError:
                continue
            C = value_norm(A, f)
            if good(points, values, A, C):
                break
        else:
            raise AssertionError(("no exact chart witness in specified attempts", profiles, s))
        rows.append({
            "profiles": profiles, "s": s, "degree": N, "p": p,
            "regime": ("monicity_hyperplane" if d == -1 else "moment_matrix"
                       if d < -1 else "monic_residual"),
            "moment_rows": -d if d < 0 else 0,
            "points": list(map(str, points)), "branch_values": list(map(str, values)),
            "residual": None if R is None else list(map(str, R)),
            "monicity_linear_coefficients": None if leading is None else list(map(str, leading)),
            "polynomial": list(map(str, f)), "critical_quotient": list(map(str, A)),
            "remaining_critical_value_polynomial": list(map(str, C)),
            "exact_open_conditions": True,
        })
    assert any(x["regime"] == "monicity_hyperplane" for x in rows)
    assert any(x["regime"] == "monic_residual" for x in rows)
    assert max(x["moment_rows"] for x in rows) >= 5
    return rows


def tree_seed(profiles, s):
    """Find one finite cycle-incidence tree; no classification is assumed."""
    N, _, _, _ = domain(profiles, s, stable=False)
    colored = [list(q) for q in profiles]+[[2] for _ in range(s)]
    initial = tuple(tuple() for _ in colored)

    def extend(remaining, used, edges, next_label):
        if not any(remaining):
            return edges
        # Same-color cycles must be disjoint; every new cycle attaches once.
        options = sorted(range(len(colored)), key=lambda i: -len(remaining[i]))
        for color in options:
            if not remaining[color]:
                continue
            q = remaining[color][0]
            parents = [0] if next_label == 0 else [
                h for h in range(next_label) if h not in used[color]]
            for parent in parents:
                if next_label == 0:
                    cy = tuple(range(q)); new_label = q
                else:
                    cy = (parent,)+tuple(range(next_label, next_label+q-1))
                    new_label = next_label+q-1
                rem = list(remaining); rem[color] = remaining[color][1:]
                us = list(used); us[color] = used[color] | frozenset(cy)
                ed = list(edges); ed[color] = edges[color]+(cy,)
                answer = extend(tuple(rem), tuple(us), tuple(ed), new_label)
                if answer is not None:
                    return answer
        return None

    tree = extend(tuple(map(tuple, colored)), tuple(frozenset() for _ in colored),
                  initial, 0)
    assert tree is not None
    raw = tuple(m.perm(N, cycles) for cycles in tree)
    seed = c.v.one.fixed_product(raw)
    c.validate(seed)
    return seed


def group_control(seed, expected_base, expected_blocks):
    N = len(seed[0])
    n, D, edges = l.group_from_transposition(seed)
    assert (n, D) == (expected_base, expected_blocks)
    assert edges == D*comb(n, 2)
    q = tuple(part[0] for part in l.to_wreath(seed, D))
    quotient_order = c.h.group_closure(q)
    assert quotient_order == D
    order = factorial(n)**D*D
    exact = c.h.group_closure(seed) if order <= 70000 else None
    if exact is not None:
        assert exact == order
    blocks = sorted(set(b.divisors(D)) | {N})
    assert c.block_divisors(seed) == blocks
    deck = 2 if n == 2 else 1
    assert c.stabilizer_size(seed) == deck
    return {"normal_transposition_base": n, "canonical_block_count": D,
            "normal_transposition_edges": edges, "quotient_order": quotient_order,
            "actual_group_order": order, "direct_group_closure": exact,
            "all_block_count_divisors": blocks, "deck": deck}


def complete_orbits():
    rows = []
    for profiles, s in CASES:
        N, p, _, passport = domain(profiles, s)
        total = c.character_total(passport)*c.type_placements(passport)
        assert total == c.comparison_total(passport)*c.type_placements(passport)
        # Delimited finite controls, not a search claimed to cover all passports.
        if total > 15000:
            continue
        seed = tree_seed(profiles, s)
        D = reduce(gcd, passport[0]) if len(profiles) == 1 else 1
        control = group_control(seed, N//D, D)
        orbit, _ = c.closure(seed)
        assert len(orbit) == total
        assert all(sorted(tuple(m.cycle_lengths(g)) for g in state) == sorted(passport)
                   for state in orbit)
        patterns = Counter(tuple(tuple(m.cycle_lengths(g)) for g in state) for state in orbit)
        assert len(patterns) == c.type_placements(passport)
        assert set(patterns.values()) == {c.character_total(passport)}
        rows.append(dict(control, profiles=profiles, s=s, p=p, degree=N,
                         literal_seed=seed, literal_orbit=total,
                         absolute_orbit=total*control["deck"]//N,
                         orbit_sha256=c.digest(orbit)))
        print(json.dumps({"degree": N, "profiles": profiles, "s": s,
                          "complete_literal_orbit": total}), flush=True)
    assert any(x["s"] == x["p"]-2 for x in rows)
    assert any(x["deck"] == 2 for x in rows)
    return rows


def colored_braid_controls():
    """Complete fixed-pattern absolute closures, using the colored braid group."""
    rows = []
    for profiles, s in ((((2, 2), (3,), (3,)), 1),
                        (((2, 2), (3,), (3,), (3,)), 1)):
        N, p, _, passport = domain(profiles, s)
        seed = tree_seed(profiles, s)
        group = group_control(seed, N, 1)
        pattern = tuple(tuple(m.cycle_lengths(g)) for g in seed)
        r = len(seed); words = set()
        for i in range(r):
            for end in range(i+1, r):
                path = tuple(range(end, i+1, -1))
                back = tuple(-x for x in path[::-1])
                word = path+(i+1, i+1)+back
                words.add(word); words.add(tuple(-x for x in word[::-1]))
        for part in set(pattern):
            positions = [i for i, q in enumerate(pattern) if q == part]
            for i, end in zip(positions, positions[1:]):
                path = tuple(range(end, i+1, -1))
                word = path+(i+1,)+tuple(-x for x in path[::-1])
                words.add(word); words.add(tuple(-x for x in word[::-1]))
        for word in words:
            assert tuple(tuple(m.cycle_lengths(g)) for g in m.braid(seed, word)) == pattern
        start = c.cyclic_label(seed)
        seen = {start}; queue = deque([start])
        while queue:
            state = queue.popleft()
            for word in words:
                new = c.cyclic_label(m.braid(state, word))
                if new not in seen:
                    seen.add(new); queue.append(new)
        ordered = c.character_total(passport)
        assert ordered == c.comparison_total(passport) == N*len(seen)
        rows.append(dict(group, profiles=profiles, s=s, p=p, degree=N,
                         literal_seed=seed, fixed_pattern_absolute_orbit=len(seen),
                         fixed_pattern_literal_domain=ordered,
                         full_type_placements=c.type_placements(passport),
                         complete_full_literal_domain=ordered*c.type_placements(passport),
                         colored_generators=len(words), orbit_sha256=c.digest(seen)))
        print(json.dumps({"degree": N, "h": p-s-2, "fixed_pattern_absolute_orbit": len(seen),
                          "complete_full_literal_domain": rows[-1]["complete_full_literal_domain"]}),
              flush=True)
    return rows


def sharpness():
    profiles, s = ((2, 2), (2, 2)), 1
    N, p, _, passport = domain(profiles, s, stable=False)
    union = set(); rows = []
    for descriptor in b.component_descriptors(6, 2, 2, 1):
        seed = b.component_seed(descriptor, 1)
        orbit, _ = c.closure(seed)
        assert not union & orbit
        union.update(orbit)
        n, D, edges = l.group_from_transposition(seed)
        q = tuple(part[0] for part in l.to_wreath(seed, D))
        qo = c.h.group_closure(q)
        order = factorial(n)**D*qo
        assert c.h.group_closure(seed) == order
        blocks = c.block_divisors(seed); deck = c.stabilizer_size(seed)
        assert len(orbit) == descriptor["literal_component_size"]
        rows.append({"literal_seed": seed, "literal_orbit": len(orbit),
                     "absolute_orbit": len(orbit)*deck//N,
                     "group_order": order, "normal_base": n, "quotient_order": qo,
                     "blocks": blocks, "deck": deck, "orbit_sha256": c.digest(orbit)})
    assert len(union) == c.character_total(passport)*c.type_placements(passport) == 243
    assert sorted((x["literal_orbit"], x["group_order"], x["deck"]) for x in rows) == [
        (27, 48, 2), (216, 720, 1)]
    return {"degree": N, "p": p, "s": s, "complete_literal_domain": 243,
            "components": rows}


def negatives():
    rejected = 0
    bad = ((((2,),), 1), (((2, 2),), 0), (((2, 2), (2, 2)), 1),
           (((3, 3, 3),), 1))
    for profiles, s in bad:
        try:
            domain(profiles, s)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError("invalid domain accepted")
    try:
        hermite_basis(((2, 2), (3,)), [0, 0, 1])
    except AssertionError:
        rejected += 1
    else:
        raise AssertionError("coincident prescribed critical points accepted")
    # A repeated remaining critical point, a branch-value collision, a repeated
    # remaining value, and a remaining value equal to a prescribed value.
    assert not good([0], [0], [1, -2, 1], [2, 3, 1])
    assert not good([0], [1, 1], [1, 1], [2, 1])
    assert not good([0], [0], [2, 3, 1], [1, -2, 1])
    assert not good([0], [1], [1, 1], [-1, 1])
    # A prescribed critical point with excess multiplicity.
    assert not good([0], [1], [0, 1], [2, 1])
    return {"invalid_domains_or_supports_rejected": rejected,
            "five_distinct_open_condition_failures_rejected": True}


def main():
    inputs = [
        "POLYNOMIAL_STABLE_HERMITE_CONNECTIVITY.md",
        "POLYNOMIAL_PASSPORT_COMPONENTS.md",
        "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
        "POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md",
        "verify_polynomial_phase_interpolation.py",
        "verify_polynomial_involution_components.py",
        "verify_polynomial_lift_transport.py",
        "verify_polynomial_passport_components.py",
        "sources/goulden-jackson-cacti-1992.pdf",
        "sources/mullane-second-kind-hurwitz-2022.pdf",
        "sources/jones-zvonkin-braid-cacti-2002.pdf",
    ]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    charts = chart_controls()
    differential = differential_controls(charts)
    orbits = complete_orbits()
    colored = colored_braid_controls()
    receipt = {
        "schema_version": 1, "kind": "stable_Hermite_arbitrary_passports",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": sha(Path(__file__).resolve()),
        "input_sha256": {p: sha(HERE/p) for p in dict.fromkeys(inputs)},
        "exact_chart_witnesses": charts, "complete_stable_orbits": orbits,
        "moment_differential_controls": differential,
        "colored_braid_complete_orbits": colored,
        "complete_stable_literal_states": sum(x["literal_orbit"] for x in orbits),
        "complete_all_literal_domains": (sum(x["literal_orbit"] for x in orbits)
                                        +sum(x["complete_full_literal_domain"] for x in colored)),
        "complete_passport_count": len(orbits)+len(colored),
        "sharpness_counterexample": sharpness(), "negative_controls": negatives(),
        "scope": "Exact rational CRT/norm examples in both regimes and bounded complete "
                 "permutation orbits. Uniform connectivity and actual groups use K1-K3; "
                 "existence/counts use the corrected primary cactus formula. "
                 "Passports below s=p-k remain outside this theorem.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    (HERE/"polynomial_stable_hermite_verification.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print("PASS: "+str(len(charts))+" exact charts; "+str(len(differential))+
          " moment differential controls; "+str(receipt["complete_passport_count"])+
          " complete stable passports representing "+str(receipt["complete_all_literal_domains"])+
          " literals; 243-tuple two-component boundary counterexample")


if __name__ == "__main__":
    main()
