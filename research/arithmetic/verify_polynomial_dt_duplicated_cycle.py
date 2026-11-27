#!/usr/bin/env python3
"""Exact rational two-component charts and bounded full orbit certificates."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_first_deficit as u

if not __debug__:
    raise RuntimeError("DT duplicated-cycle checker requires active assertions")
HERE = Path(__file__).resolve().parent
h, c, m, l, j = u.h, u.c, u.m, u.l, u.j


def parameter_polynomials(b):
    assert isinstance(b, int) and b >= 2
    A = [F((at+1)*(b-at)) for at in range(b)]
    B = [F(1)]*(b+1)
    D = j.add(j.mul([b+2, -b], [0]*(b+1)+[1]), [b, -(b+2)])
    E = j.scale(j.sub([1], [0]*(b+1)+[1]), b)
    assert D == j.scale(j.mul(j.power([-1, 1], 3), A), -1)
    assert E == j.scale(j.mul([-1, 1], B), -b)
    assert j.gcd_poly(D, E) == [-1, 1]
    assert A == A[::-1] and B == B[::-1]
    return A, B, D, E


def primitive_parameters(b, t):
    A, B, _, _ = parameter_polynomials(b)
    t = F(t)
    assert t and t != 1 and j.evaluate(A, t)
    V = b*j.evaluate(B, t)/((t-1)**2*j.evaluate(A, t))
    assert V
    U = V*(1-t)
    return 1+U, V


def chart(b, S, V):
    assert isinstance(b, int) and b >= 2
    S, V = F(S), F(V)
    N = 2*b+2; P = [F(0), F(-1), F(1)]; B = [V, -S, F(1)]
    assert V and 1-S+V and S*S-4*V
    assert j.gcd_poly(P, B) == [1]
    P2, Bb = j.power(P, 2), j.power(B, b)
    E = j.mul(P2, j.inverse_mod(P2, Bb))
    assert j.mod(E, P2) == [0] and j.mod(j.sub(E, [1]), Bb) == [0]
    assert len(E) == N+1
    L = E[-1]
    f = j.scale(E, 1/L); v = 1/L
    Z = j.mul(P, j.power(B, b-1))
    A, rest = j.divmod_poly(j.derivative(f), j.scale(Z, N))
    assert rest == [0] and len(A) == 2 and A[-1] == 1
    alpha = -A[0]; w = j.evaluate(f, alpha)
    assert j.gcd_poly(A, j.mul(P, B)) == [1] and w and w != v
    assert j.mod(f, P2) == [0] and j.mod(j.sub(f, [v]), Bb) == [0]
    C = h.value_norm(A, f)
    assert C == [-w, 1]
    return f, v, A, w


def exact_charts():
    rows = []
    for b in list(range(2, 21))+[32]:
        S, V = primitive_parameters(b, 2)
        expected = F(b*(2**(b+1)-1), (b-2)*2**(b+1)+b+4)
        assert V == expected and S == 1-V and 1 < V < 3
        assert S*S-4*V < 0
        for sector, ss, vv in (("primitive", S, V), ("imprimitive", F(1), F(2))):
            f, value, A, w = chart(b, ss, vv)
            alpha = -A[0]
            assert 0 < alpha < 1 and w > 0
            if sector == "imprimitive":
                assert alpha == F(1, 2) and j.compose(f, [1, -1]) == f
            rows.append({"b": b, "degree": 2*b+2, "component": sector,
                         "S": str(ss), "V": str(vv), "B_value": str(value),
                         "polynomial": list(map(str, f)),
                         "critical_quotient": list(map(str, A)),
                         "T_value": str(w), "all_exact_open_conditions": True})
    return rows


def residue_and_action_controls():
    rows = []
    for b in range(2, 13):
        A, B, D, E = parameter_polynomials(b)
        for t in (F(2), F(3), F(4, 3)):
            S, V = primitive_parameters(b, t)
            f, v, _, _ = chart(b, S, V)
            Q = j.mul(j.power([0, -1, 1], 2), j.power([V, -S, 1], b))
            delta = h.residue_at([1], Q, F(0), 2)+h.residue_at([1], Q, F(1), 2)
            assert delta == 0
            U = S-1; W = V-U
            residue_formula = (2*V+b*S)/V**(b+1)-(2*V-(b+2)*U+b)/W**(b+1)
            assert residue_formula == delta
            assert j.evaluate(D, t)*V == j.evaluate(E, t)
            S2, V2 = primitive_parameters(b, 1/t)
            assert (S2, V2) == (2-S, W)
            f2, v2, _, _ = chart(b, S2, V2)
            assert f2 == j.compose(f, [1, -1]) and v2 == v
            q = t+1/t
            if b == 2:
                L0 = S-S*S+2*V
                qp = L0*L0/(V*W)-2
                assert qp == 1/(q+2)-2
                assert (qp+2)+1/(qp+2) == (q+2)+1/(q+2)
            rows.append({"b": b, "t": str(t), "primitive_quotient_q": str(q),
                         "residue_zero": True, "reflected_polynomial_equal": True})
    for V in (F(2), F(3), F(4, 3)):
        V2 = V/(4*V-1)
        assert V2/(4*V2-1) == V
        assert 4*V2*V2/(4*V2-1) == 4*V*V/(4*V-1)
        chart(2, 1, V); chart(2, 1, V2)
    return {"primitive_residue_and_reflection_rows": rows,
            "equal_profile_imprimitive_color_involutions": 3,
            "equal_profile_primitive_color_involutions": 3}


def group_control(seed, primitive):
    N = len(seed[0]); B, d, edges = l.group_from_transposition(seed)
    assert edges == d*B*(B-1)//2
    if primitive:
        assert (B, d) == (N, 1)
        order = factorial(N); blocks = [1, N]; deck = 1
        quotient_edges = None
    else:
        assert B == 2 and N == 2*d
        quotient = tuple(entry[0] for entry in l.to_wreath(seed, d))
        # Its first heavy quotient is a transposition. Its normal graph
        # is connected, which certifies S_d without enumerating d! elements.
        qb, qd, quotient_edges = l.group_from_transposition(quotient)
        assert (qb, qd, quotient_edges) == (d, 1, d*(d-1)//2)
        order = 2**d*factorial(d); blocks = [1, d, N]; deck = 2
    assert c.block_divisors(seed) == blocks and c.stabilizer_size(seed) == deck
    return {"actual_group_order": order, "normal_transposition_base": B,
            "canonical_block_count": d, "normal_transposition_edges": edges,
            "quotient_normal_transposition_edges": quotient_edges,
            "all_block_count_divisors": blocks, "deck": deck}


def complete_components():
    rows = []
    for b in range(2, 13):
        N, _, _, passport = u.first_deficit(((2, 2), (b, b)), 1)
        d = b+1; W = 3*(2-(b == 2))
        total = c.character_total(passport)
        assert total == c.comparison_total(passport) == 3*d*d*(2*b-1)
        components = []; union = set(); weighted_total = 0
        for primitive, seed in ((True, h.tree_seed(((2, 2), (b, b)), 1)),
                                (False, u.imprimitive_seed(2, b))):
            control = group_control(seed, primitive)
            orbit, generators = u.absolute_colored_closure(seed)
            assert not union & orbit; union.update(orbit)
            expected_absolute = d*(3*b-2) if primitive else d
            assert len(orbit) == expected_absolute
            literal = len(orbit)*N//control["deck"]
            weighted_total += literal
            components.append(dict(control, component="primitive" if primitive else "imprimitive",
                                   literal_seed=seed,
                                   fixed_pattern_absolute_orbit=len(orbit),
                                   fixed_pattern_literal_orbit=literal,
                                   full_literal_orbit=W*literal,
                                   full_absolute_orbit=W*len(orbit),
                                   colored_generators=generators,
                                   orbit_sha256=c.digest(orbit)))
        assert weighted_total == total and c.type_placements(passport) == W
        rows.append({"b": b, "degree": N, "full_type_placements": W,
                     "ordered_character_total": total,
                     "complete_full_literal_domain": W*total, "components": components})
        print(json.dumps({"b": b, "degree": N, "fixed_pattern_absolute_components":
                          [x["fixed_pattern_absolute_orbit"] for x in components]}), flush=True)
    return rows


def negative_controls():
    rejected = 0
    for fn, args in ((primitive_parameters, (1, F(2))),
                     (primitive_parameters, (3, F(0))),
                     (primitive_parameters, (3, F(1))),
                     (primitive_parameters, (2, F(-1))),
                     (chart, (3, F(1), F(0))),
                     (chart, (3, F(2), F(1))),
                     (chart, (3, F(4), F(5)))):
        try:
            fn(*args)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError(("invalid coefficient/parameter certificate accepted", args))
    assert rejected == 7
    A, _, D, _ = parameter_polynomials(3)
    assert j.add(D, [1]) != j.scale(j.mul(j.power([-1, 1], 3), A), -1)
    return {"invalid_parameters_or_charts_rejected": rejected,
            "corrupted_denominator_identity_rejected": True}


def main():
    inputs = [
        "POLYNOMIAL_DT_DUPLICATED_CYCLE_COMPONENTS.md",
        "POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md",
        "POLYNOMIAL_STABLE_HERMITE_CONNECTIVITY.md",
        "POLYNOMIAL_PASSPORT_COMPONENTS.md", "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
        "verify_polynomial_first_deficit.py", "verify_polynomial_stable_hermite.py",
        "verify_polynomial_phase_interpolation.py", "verify_polynomial_involution_components.py",
        "verify_polynomial_lift_transport.py", "verify_polynomial_passport_components.py",
        "sources/bogatyrev-gendron-pell-abel-2025.pdf",
        "sources/wajnryb-two-special-fibers-1995.pdf",
    ]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    charts = exact_charts(); actions = residue_and_action_controls()
    components = complete_components()
    receipt = {
        "schema_version": 1, "kind": "DT_duplicated_cycle_uniform_two_components",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": sha(Path(__file__).resolve()),
        "input_sha256": {p: sha(HERE/p) for p in dict.fromkeys(inputs)},
        "exact_rational_charts": charts, "residue_and_label_actions": actions,
        "complete_component_controls": components,
        "complete_passport_count": len(components),
        "complete_literal_domains": sum(x["complete_full_literal_domain"] for x in components),
        "negative_controls": negative_controls(),
        "scope": "Uniform all-b classification, rationality and nonemptiness use written V1-V4. Exact bounded controls:40 coefficient charts at b2..20 and32;33 residue/reflection points, equal-profile color actions; complete quotient-braid coverage for b2..12. No arbitrary a,b>=3 or general first-deficit classification is inferred.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    (HERE/"polynomial_dt_duplicated_cycle_verification.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print("PASS:40 rational charts;33 residue/reflection controls;"+
          str(len(components))+" complete two-component passports, "+
          str(receipt["complete_literal_domains"])+" literal tuples represented")


if __name__ == "__main__":
    main()
