#!/usr/bin/env python3
"""Exact first-deficit monodromy, determinant and bounded component controls."""
from collections import deque
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_stable_hermite as h

if not __debug__:
    raise RuntimeError("First-deficit checker requires active assertions")
HERE = Path(__file__).resolve().parent
c, m, l, j = h.c, h.m, h.l, h.j


def first_deficit(profiles, s):
    N, p, I, passport = h.domain(profiles, s, stable=False)
    k = len(profiles)
    assert k >= 2 and s == p-k-1
    return N, p, I, passport


def colored_generators(pattern):
    r = len(pattern)
    words = set()
    for i in range(r):
        for end in range(i+1, r):
            path = tuple(range(end, i+1, -1))
            word = path+(i+1, i+1)+tuple(-x for x in path[::-1])
            words.add(word)
            words.add(tuple(-x for x in word[::-1]))
    for part in set(pattern):
        positions = [i for i, q in enumerate(pattern) if q == part]
        for i, end in zip(positions, positions[1:]):
            path = tuple(range(end, i+1, -1))
            word = path+(i+1,)+tuple(-x for x in path[::-1])
            words.add(word)
            words.add(tuple(-x for x in word[::-1]))
    return sorted(words)


def absolute_colored_closure(seed):
    pattern = tuple(tuple(m.cycle_lengths(g)) for g in seed)
    words = colored_generators(pattern)
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
    return seen, len(words)


def group_control(seed, base, blocks, quotient_order):
    N = len(seed[0])
    B, d, edges = l.group_from_transposition(seed)
    assert (B, d) == (base, blocks)
    assert edges == d*B*(B-1)//2
    q = tuple(part[0] for part in l.to_wreath(seed, d))
    assert c.h.group_closure(q) == quotient_order
    block_counts = [1, N] if d == 1 else [1, d, N]
    assert c.block_divisors(seed) == block_counts
    deck = 2 if B == 2 else 1
    assert c.stabilizer_size(seed) == deck
    return {"normal_transposition_base": B, "canonical_block_count": d,
            "normal_transposition_support_edges": edges,
            "quotient_group_order": quotient_order,
            "actual_group_order": factorial(B)**d*quotient_order,
            "all_block_count_divisors": block_counts, "deck": deck}


def imprimitive_seed(a, b):
    d = a+b-1
    A = m.perm(d, [tuple(range(a))])
    B = m.perm(d, [(0,)+tuple(range(a, d))])
    outer = c.v.one.fixed_product((A, B))
    seed = l.free_seed(outer, 2)
    c.validate(seed)
    return seed


def family_controls():
    rows = []
    for a in range(2, 6):
        for b in range(a, 6):
            profiles = ((a, a), (b, b))
            N, p, _, passport = first_deficit(profiles, 1)
            d = a+b-1
            W = c.type_placements(passport)
            assert W == 3*(2-(a == b))
            D = (2*a-1)*(2*b-1)
            ordered_total = c.character_total(passport)
            assert ordered_total == c.comparison_total(passport) == d*d*D
            prim = h.tree_seed(profiles, 1)
            prim_group = group_control(prim, N, 1, 1)
            prim_orbit, generators = absolute_colored_closure(prim)
            assert len(prim_orbit)*N == d*d*(D-1)
            imp = imprimitive_seed(a, b)
            qo = factorial(d)//(2 if a % 2 and b % 2 else 1)
            imp_group = group_control(imp, 2, d, qo)
            imp_orbit, _ = absolute_colored_closure(imp)
            assert len(imp_orbit) == d and not prim_orbit & imp_orbit
            assert N*len(prim_orbit)+(N//2)*len(imp_orbit) == ordered_total
            components = []
            for seed, orbit, control in ((prim, prim_orbit, prim_group),
                                          (imp, imp_orbit, imp_group)):
                components.append(dict(
                    control, literal_seed=seed, fixed_pattern_absolute_orbit=len(orbit),
                    fixed_pattern_literal_orbit=len(orbit)*N//control["deck"],
                    full_absolute_orbit=W*len(orbit),
                    full_literal_orbit=W*len(orbit)*N//control["deck"],
                    orbit_sha256=c.digest(orbit)))
            rows.append({"a": a, "b": b, "degree": N, "p": p, "k": 2, "s": 1,
                         "ordered_character_total": ordered_total,
                         "full_type_placements": W,
                         "complete_full_literal_domain": W*ordered_total,
                         "colored_generators": generators, "components": components})
            print(json.dumps({"a": a, "b": b, "fixed_pattern_absolute_components":
                              [len(prim_orbit), len(imp_orbit)]}), flush=True)
    return rows


def mixed_controls():
    rows = []
    for profiles, s in ((((2, 2), (2, 2), (3,)), 1),
                        (((2, 2), (2, 2), (2, 2)), 2)):
        N, p, _, passport = first_deficit(profiles, s)
        seed = h.tree_seed(profiles, s)
        control = group_control(seed, N, 1, 1)
        orbit, generators = absolute_colored_closure(seed)
        ordered = c.character_total(passport)
        assert ordered == c.comparison_total(passport) == N*len(orbit)
        W = c.type_placements(passport)
        rows.append(dict(control, profiles=profiles, degree=N, p=p, k=len(profiles),
                         s=s, literal_seed=seed, fixed_pattern_absolute_orbit=len(orbit),
                         fixed_pattern_literal_domain=ordered, full_type_placements=W,
                         complete_full_literal_domain=W*ordered,
                         full_absolute_orbit=W*len(orbit), colored_generators=generators,
                         orbit_sha256=c.digest(orbit)))
        print(json.dumps({"degree": N, "profiles": profiles, "s": s,
                          "fixed_pattern_absolute_orbit": len(orbit),
                          "full_literal_domain": W*ordered}), flush=True)
    assert [x["complete_full_literal_domain"] for x in rows] == [38400, 1771470]
    return rows


def degree6_determinant():
    # Each cleared numerator has degree at most four in x and in y.
    # The 5-by-5 grid therefore proves the displayed polynomial identity.
    checked = 0
    for x in range(2, 7):
        for y in range(-5, 0):
            Q = j.power(j.from_roots([0, 1, x, y]), 2)
            residue = h.residue_at([1], Q, F(0), 2)+h.residue_at([1], Q, F(1), 2)
            U, V = F(x+y-1), F(x*y)
            P = -3*V*V+U*(2*U+3)*V-(U+1)*U*U
            assert residue*V**3*(V-U)**3 == 2*U*P
            checked += 1
    Un = [1, -3, 3]; den = [-1, 2]
    Vn = j.mul([0, 1], Un); Sn = j.add(Un, den)
    identity = j.add(j.add(j.scale(j.mul(j.power(Vn, 2), den), -3),
                           j.mul(j.mul(Un, j.add(j.scale(Un, 2), j.scale(den, 3))), Vn)),
                     j.scale(j.mul(Sn, j.power(Un, 2)), -1))
    assert identity == [0]
    disc_numerator = j.sub(j.power(Sn, 2), j.scale(j.mul(Vn, den), 4))
    quartic = j.scale(j.mul(j.mul([0, 1], [-1, 1]), [4, -15, 15]), -1)
    assert disc_numerator == quartic
    assert j.gcd_poly(quartic, j.derivative(quartic)) == [1]
    points = list(map(F, (0, 1, 2, -1)))
    Q, E, _, _, _ = h.hermite_basis(((2, 2), (2, 2)), points)
    matrix = h.moment_matrix(E, Q, 1)
    values, kernel, pivots = h.linear_solve(matrix, [0, 1])
    assert len(pivots) == 1
    f = j.add(j.scale(E[0], values[0]), j.scale(E[1], values[1]))
    Z = j.from_roots(points)
    A, rest = j.divmod_poly(j.derivative(f), j.scale(Z, 6))
    assert rest == [0]
    C = h.value_norm(A, f)
    assert h.good(points, values, A, C)
    assert j.compose(f, [1, -1]) == f
    return {"cleared_identity_degree_bounds": [4, 4], "identity_grid_size": checked,
            "rational_parameter_U_numerator": Un, "rational_parameter_denominator": den,
            "ordered_root_discriminant_numerator": list(map(str, quartic)),
            "squarefree_discriminant_numerator": True,
            "geometric_irreducibility_argument": "Written U4: nonsquare rational-function discriminant over C(w).",
            "imprimitive_exact_chart": {
                "critical_points": list(map(str, points)), "branch_values": list(map(str, values)),
                "polynomial": list(map(str, f)), "critical_quotient": list(map(str, A)),
                "critical_value_polynomial": list(map(str, C)),
                "source_deck_involution": "z -> 1-z"}}


def negative_controls():
    rejected = 0
    for profiles, s in ((((2, 2), (2, 2)), 2), (((2, 2), (2, 2)), 0),
                        (((2,), (2, 2)), 1), (((2, 2, 2),), 1)):
        try:
            first_deficit(profiles, s)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError("invalid first-deficit input accepted")
    Q = j.power(j.from_roots([0, 1, 2, 3]), 2)
    residue = h.residue_at([1], Q, F(0), 2)+h.residue_at([1], Q, F(1), 2)
    U, V = F(4), F(6)
    P = -3*V*V+U*(2*U+3)*V-(U+1)*U*U
    assert residue*V**3*(V-U)**3 != 2*U*(P+1)
    assert rejected == 4
    return {"invalid_hypothesis_domains_rejected": rejected,
            "corrupted_determinant_factor_identity_rejected": True,
            "geometric_factorization_scope": "No rational factorization is used to certify geometric factorization; U4 uses its nonsquare function-field proof."}


def main():
    inputs = [
        "POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md",
        "POLYNOMIAL_STABLE_HERMITE_CONNECTIVITY.md",
        "POLYNOMIAL_PASSPORT_COMPONENTS.md", "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
        "verify_polynomial_stable_hermite.py", "verify_polynomial_phase_interpolation.py",
        "verify_polynomial_involution_components.py", "verify_polynomial_lift_transport.py",
        "verify_polynomial_passport_components.py",
        "sources/jones-zvonkin-braid-cacti-2002.pdf",
        "sources/goulden-jackson-cacti-1992.pdf",
        "sources/kochetkov-zolotarev-trees-2012.pdf",
    ]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    family = family_controls()
    mixed = mixed_controls()
    receipt = {
        "schema_version": 1, "kind": "first_deficit_monodromy_and_component_reduction",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": sha(Path(__file__).resolve()),
        "input_sha256": {p: sha(HERE/p) for p in dict.fromkeys(inputs)},
        "complete_family_controls": family, "complete_mixed_controls": mixed,
        "complete_passport_count": len(family)+len(mixed),
        "complete_full_literal_domains": sum(x["complete_full_literal_domain"]
                                           for x in family+mixed),
        "degree6_determinant": degree6_determinant(), "negative_controls": negative_controls(),
        "scope": "Uniform U1-U3 monodromy/imprimitive-component/determinant reduction are written proofs, not finite extrapolations. U4's geometric irreducibility is a written quadratic-cover argument with exact algebraic controls. U5 is complete only for the ten listed2<=a<=b<=5 cases and two additional mixed passports. General determinant factor orbits and separate primitive degrees remain open.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    (HERE/"polynomial_first_deficit_verification.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print("PASS: "+str(receipt["complete_passport_count"])+" complete first-deficit passports, "+
          str(receipt["complete_full_literal_domains"])+" literal tuples represented; "
          "degree6 determinant/geometric-factor certificates")


if __name__ == "__main__":
    main()
