#!/usr/bin/env python3
"""Reduced residue equations, discriminant certificates and full orbit controls."""
from collections import deque
from fractions import Fraction as F
from functools import reduce
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_first_deficit as u

if not __debug__:
    raise RuntimeError("Duplicated-three-cycle checker requires active assertions")
HERE = Path(__file__).resolve().parent
h, c, m, l, j = u.h, u.c, u.m, u.l, u.j


# Polynomials in V, with coefficients in Q[t].
def vtrim(p):
    p = [j.trim(q) for q in p]
    while len(p) > 1 and p[-1] == [0]:
        p.pop()
    return p or [[F(0)]]


def vadd(p, q):
    return vtrim([j.add(p[i] if i < len(p) else [0],
                        q[i] if i < len(q) else [0])
                  for i in range(max(len(p), len(q)))])


def vscale(p, scalar):
    return vtrim([j.mul(q, scalar) for q in p])


def vmul(p, q):
    out = [[F(0)] for _ in range(len(p)+len(q)-1)]
    for i, x in enumerate(p):
        for k, y in enumerate(q):
            out[i+k] = j.add(out[i+k], j.mul(x, y))
    return vtrim(out)


def vpower(p, n):
    out = [[F(1)]]
    for _ in range(n):
        out = vmul(out, p)
    return out


def residue_terms(a, b):
    assert isinstance(a, int) and isinstance(b, int) and a >= 2 and b >= 2
    terms = {}
    for ell in range(a):
        for k in range(ell//2+1):
            coefficient = ((-1)**k*comb(2*a-ell-2, a-ell-1)
                           *comb(b+ell-k-1, ell-k)*comb(ell-k, k))
            exponents = (ell-2*k, a-ell+k-1)
            terms[exponents] = terms.get(exponents, F(0))+coefficient
    return {key: value for key, value in terms.items() if value}


def substitute_terms(terms, S, V):
    out = [[F(0)]]
    for (r, s), coefficient in terms.items():
        out = vadd(out, vscale(vmul(vpower(S, r), vpower(V, s)), [coefficient]))
    return out


def reduced_curve(a, b):
    terms = residue_terms(a, b)
    V = [[F(0)], [F(1)]]
    plus = [[F(1)], [F(1), F(-1)]]
    minus = [[F(1)], [F(-1), F(1)]]
    tV = [[F(0)], [F(0), F(1)]]
    R = vadd(vscale(substitute_terms(terms, plus, V), [0]*(a+b-1)+[1]),
             vscale(substitute_terms(terms, minus, tV), [-1]))
    result = []
    for coefficient in R:
        q, rest = j.divmod_poly(coefficient, [-1, 1])
        assert rest == [0]
        result.append(q)
    result = vtrim(result)
    assert len(result) == a
    assert reduce(j.gcd_poly, result) == [1]
    C = comb(a+b-2, a-1)
    assert [j.evaluate(coef, 0) for coef in result] == [
        F(C*(-1)**i*comb(a-1, i)) for i in range(a)]
    return result


def scalar_H(a, b, S, V):
    return sum(value*S**r*V**s for (r, s), value in residue_terms(a, b).items())


def evaluate_curve(curve, V, t):
    return sum(j.evaluate(coefficient, t)*V**i for i, coefficient in enumerate(curve))


def general_residue_controls():
    rows = []
    for a in range(2, 8):
        for b in range(2, 10):
            curve = reduced_curve(a, b)
            n = a-1; degree = a+b-1
            C, Cp = comb(b+n-1, n), comb(b+n-2, n-1)
            tm = [0]*degree+[1]
            R0 = j.scale(j.sub(tm, [1]), C)
            R1 = j.add(j.scale(j.mul([1, -1], j.add(tm, [1])), C*n),
                        j.scale(j.sub(tm, [0, 1]), 2*Cp))
            for i, expected in enumerate((R0, R1)):
                assert j.mul(curve[i], [-1, 1]) == expected
            for S, V in ((F(0), F(1)), (F(1), F(2)),
                         (F(-1), F(3)), (F(2), F(5))):
                W = 1-S+V; t = W/V
                B = [V, -S, 1]
                Q = j.mul(j.power([0, -1, 1], a), j.power(B, b))
                actual = h.residue_at([1], Q, F(0), a)+h.residue_at([1], Q, F(1), a)
                formula = (-1)**a*(scalar_H(a, b, S, V)/V**degree
                                  -scalar_H(a, b, 2-S, W)/W**degree)
                assert actual == formula
                assert actual*W**degree == (-1)**a*(t-1)*evaluate_curve(curve, V, t)
            rows.append({"a": a, "b": b, "degree_in_V": a-1,
                         "content": "1", "independent_residue_points": 4,
                         "first_two_coefficients_and_t0_identity": True})
    return rows


def quadratic_formula(b):
    assert isinstance(b, int) and b >= 3
    C = F(b*(b+1), 2); one_minus_t = [1, -1]; tm = [0]*(b+2)+[1]
    r0 = j.scale(j.sub(tm, [1]), C)
    r1 = j.scale(j.sub(j.mul([b+3, -b-1], tm), [-b-1, b+3]), b)
    A = j.add(j.add([6], j.scale(one_minus_t, 3*b)),
              j.scale(j.power(one_minus_t, 2), C))
    B = j.add(j.add([0, 0, 6], j.scale(j.mul([0, 1], one_minus_t), -3*b)),
              j.scale(j.power(one_minus_t, 2), C))
    r2 = j.sub(j.mul(A, tm), B)
    result = []
    for coefficient in (r0, r1, r2):
        quotient, rest = j.divmod_poly(coefficient, [-1, 1])
        assert rest == [0]
        result.append(quotient)
    return result


def verify_discriminant_jet(curve, b):
    disc = j.sub(j.power(curve[1], 2), j.scale(j.mul(curve[0], curve[2]), 4))
    assert disc[0] == 0 and disc[1] == 2*b*b*(b+1)
    return disc


def quadratic_controls():
    rows = []
    for b in list(range(3, 21))+[32, 64]:
        curve = quadratic_formula(b)
        assert curve == reduced_curve(3, b)
        disc = verify_discriminant_jet(curve, b)
        assert j.gcd_poly(disc, [0, 0, 1]) == [0, 1]
        rows.append({"b": b, "quadratic_coefficients_Qt": [list(map(str, q)) for q in curve],
                     "discriminant_degree": len(disc)-1,
                     "discriminant_t0_order": 1, "discriminant_first_coefficient": str(disc[1])})
    # Independent exact identity in Q[b], not sampling the parameter b.
    C = [F(0), F(1, 2), F(1, 2)]
    r0, r1 = j.scale(C, -1), j.scale(C, 2)
    derivative_r1, derivative_r2 = [0, -3, -1], [0, 4, 1]
    jet = j.sub(j.scale(j.mul(r1, derivative_r1), 2),
                j.scale(j.mul(r0, derivative_r2), 4))
    assert jet == [0, 0, 2, 2]
    return {"rows": rows, "symbolic_first_jet_over_Qb": list(map(str, jet)),
            "symbolic_identity": "2*b^2*(b+1)",
            "limits": "Odd valuation is an all-b written proof; the displayed finite polynomials are exact controls."}


def alternating_support_certificate(quotient):
    degree = len(quotient[0])
    three = next(g for g in quotient
                 if m.cycle_lengths(g) == [1]*(degree-3)+[3])
    start = tuple(sorted(next(cy for cy in c.cycles(three, False) if len(cy) == 3)))
    seen = {start}; queue = deque([start])
    generators = quotient+tuple(m.inverse(g) for g in quotient)
    while queue:
        support = queue.popleft()
        for g in generators:
            new = tuple(sorted(g[x] for x in support))
            if new not in seen:
                seen.add(new); queue.append(new)
    assert len(seen) == comb(degree, 3)
    even = all(c.index(g) % 2 == 0 for g in quotient)
    return factorial(degree)//(2 if even else 1), len(seen), even


def group_control(seed, primitive, b):
    N = len(seed[0]); base, blocks, edges = l.group_from_transposition(seed)
    assert edges == blocks*base*(base-1)//2
    if primitive:
        assert (base, blocks) == (N, 1)
        order = factorial(N); all_blocks = [1, N]; deck = 1
        triples = None; quotient_even = None
    else:
        assert base == 2 and blocks == b+2
        quotient = tuple(entry[0] for entry in l.to_wreath(seed, blocks))
        quotient_order, triples, quotient_even = alternating_support_certificate(quotient)
        assert quotient_even == bool(b % 2)
        order = 2**blocks*quotient_order
        all_blocks = [1, blocks, N]; deck = 2
    assert c.block_divisors(seed) == all_blocks and c.stabilizer_size(seed) == deck
    return {"actual_group_order": order, "normal_transposition_base": base,
            "canonical_block_count": blocks, "normal_transposition_edges": edges,
            "quotient_all_three_subsets": triples, "quotient_even": quotient_even,
            "all_block_count_divisors": all_blocks, "deck": deck}


def complete_components():
    rows = []
    for b in range(3, 13):
        N, _, _, passport = u.first_deficit(((3, 3), (b, b)), 1)
        degree = b+2; W = 3*(2-(b == 3))
        total = c.character_total(passport)
        assert total == c.comparison_total(passport) == 5*degree**2*(2*b-1)
        components = []; union = set(); weighted_total = 0
        for primitive, seed in ((True, h.tree_seed(((3, 3), (b, b)), 1)),
                                (False, u.imprimitive_seed(3, b))):
            control = group_control(seed, primitive, b)
            orbit, generators = u.absolute_colored_closure(seed)
            assert not union & orbit; union.update(orbit)
            expected = degree*(5*b-3) if primitive else degree
            assert len(orbit) == expected
            literal = N*len(orbit)//control["deck"]
            weighted_total += literal
            components.append(dict(control, component="primitive" if primitive else "imprimitive",
                                   literal_seed=seed, fixed_pattern_absolute_orbit=len(orbit),
                                   fixed_pattern_literal_orbit=literal,
                                   full_literal_orbit=W*literal, full_absolute_orbit=W*len(orbit),
                                   colored_generators=generators, orbit_sha256=c.digest(orbit)))
        assert weighted_total == total and W == c.type_placements(passport)
        rows.append({"b": b, "degree": N, "ordered_character_total": total,
                     "full_type_placements": W, "complete_full_literal_domain": W*total,
                     "components": components})
        print(json.dumps({"b": b, "degree": N, "fixed_pattern_absolute_components":
                          [x["fixed_pattern_absolute_orbit"] for x in components]}), flush=True)
    return rows


def negative_controls():
    rejected = 0
    for fn, args in ((reduced_curve, (1, 3)), (reduced_curve, (3, 1)),
                     (quadratic_formula, (2,))):
        try:
            fn(*args)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError("invalid parameter domain accepted")
    corrupted = quadratic_formula(3)
    corrupted[0] = j.add(corrupted[0], [1])
    try:
        verify_discriminant_jet(corrupted, 3)
    except AssertionError:
        rejected += 1
    else:
        raise AssertionError("corrupted discriminant certificate accepted")
    # The branch point used for valuation is a collision configuration.
    assert j.gcd_poly([0, -1, 1], [1, -2, 1]) != [1]
    assert rejected == 4
    return {"invalid_domains_or_corrupted_jets_rejected": rejected,
            "t0_not_claimed_admissible": True}


def main():
    inputs = [
        "POLYNOMIAL_DUPLICATED_THREE_CYCLE_COMPONENTS.md",
        "POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md",
        "POLYNOMIAL_DT_DUPLICATED_CYCLE_COMPONENTS.md",
        "POLYNOMIAL_STABLE_HERMITE_CONNECTIVITY.md",
        "POLYNOMIAL_PASSPORT_COMPONENTS.md", "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
        "verify_polynomial_first_deficit.py", "verify_polynomial_stable_hermite.py",
        "verify_polynomial_phase_interpolation.py", "verify_polynomial_involution_components.py",
        "verify_polynomial_lift_transport.py", "verify_polynomial_passport_components.py",
        "sources/goulden-jackson-cacti-1992.pdf",
        "sources/bogatyrev-gendron-pell-abel-2025.pdf",
        "sources/wajnryb-two-special-fibers-1995.pdf",
    ]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    general = general_residue_controls(); quadratic = quadratic_controls()
    components = complete_components()
    receipt = {
        "schema_version": 1, "kind": "duplicated_three_cycle_uniform_components",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": sha(Path(__file__).resolve()),
        "input_sha256": {p: sha(HERE/p) for p in dict.fromkeys(inputs)},
        "general_reduced_curve_controls": general, "quadratic_controls": quadratic,
        "complete_component_controls": components, "complete_passport_count": len(components),
        "complete_literal_domains": sum(x["complete_full_literal_domain"] for x in components),
        "negative_controls": negative_controls(),
        "scope": "Uniform geometric irreducibility at a3 and complete all-b classification use W1-W3. Exact controls:48 general residue models,192 independent local-residue points,20 quadratic discriminants throughb64 and a symbolic Q[b] first-jet identity;10 complete two-component passports b3..12 via quotient-braid closures. No general a>=4 irreducibility, rationality or genus assertion is inferred.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    (HERE/"polynomial_duplicated_three_cycle_verification.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print("PASS:48 reduced residue models;20 quadratic/discriminant controls;"+
          str(len(components))+" complete two-component passports, "+
          str(receipt["complete_literal_domains"])+" literal tuples represented")


if __name__ == "__main__":
    main()
