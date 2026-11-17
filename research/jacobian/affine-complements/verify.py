"""Exact controls for PF2RULED1; geometric theorems remain written/imported.

No finite quartic map or compatible completion is supplied by this checker.
"""
from pathlib import Path
from functools import reduce
from itertools import combinations
import hashlib
import json
import platform
import time
import sympy as S

HERE = Path(__file__).resolve().parent
s, t, a, b, e, m, L, r = S.symbols("s t a b e m L r")
CHECKS = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def zero(name, expression):
    check(name, S.expand(expression) == 0)


def homogeneous_gcd(forms):
    polys = [S.Poly(f, s, t, domain=S.QQ) for f in forms if f != 0]
    return reduce(S.gcd, polys).monic() if polys else S.Poly(0, s, t)


def support_degree(poly):
    if poly.is_zero:
        raise ValueError("identical columns: coincidence is not finite")
    repeated = reduce(S.gcd, [poly, poly.diff(s), poly.diff(t)])
    return poly.exquo(repeated).total_degree()


def recognize_columns(splitting, degrees, columns):
    """Return criterion data for specified split-bundle homogeneous columns."""
    if len(columns) != 2 or len(degrees) != 2 or len(splitting) < 2:
        raise ValueError("need two columns in rank at least two")
    for d, column in zip(degrees, columns):
        if len(column) != len(splitting):
            raise ValueError("column rank")
        for aj, value in zip(splitting, column):
            if value == 0:
                continue
            poly = S.Poly(value, s, t, domain=S.QQ)
            if not poly.is_homogeneous or poly.total_degree() != aj+d:
                raise ValueError("incorrect homogeneous entry degree")
        if homogeneous_gcd(column).total_degree() != 0 or all(v == 0 for v in column):
            raise ValueError("column has a projective zero")
    v, w = columns
    minors = [S.expand(v[i]*w[j]-v[j]*w[i])
              for i, j in combinations(range(len(splitting)), 2)]
    g = homogeneous_gcd(minors)
    points = support_degree(g)
    adjacent = abs(degrees[0]-degrees[1]) == 1
    return dict(splitting=splitting, degrees=degrees,
                columns=[[str(f) for f in c] for c in columns],
                minors=list(map(str, minors)), homogeneous_gcd=str(g.as_expr()),
                coincidence_length=g.total_degree(), geometric_support_points=points,
                adjacent_classes=adjacent,
                affine_space_by_imported_criterion=adjacent and points == 1)


def must_reject(name, splitting, degrees, columns):
    try:
        recognize_columns(splitting, degrees, columns)
    except ValueError as exc:
        check(name, True)
        return str(exc)
    raise AssertionError(name)


def intersect(x, y):
    """Classes are (coefficient of f, coefficient of h) on F_a."""
    return S.expand(x[0]*y[1]+x[1]*y[0]+a*x[1]*y[1])


def main():
    start = time.monotonic()
    manifest = json.loads((HERE/"sources/manifest.json").read_text())
    for source in manifest["sources"]:
        raw = (HERE/"sources"/source["file"]).read_bytes()
        check("source integrity "+source["file"],
              hashlib.sha256(raw).hexdigest() == source["sha256"]
              and len(raw) == source["bytes"])

    controls = {}
    controls["F1_contact_four"] = recognize_columns(
        [0, 1], [1, 2], [[s, t**2], [s*t, s**3+t**3]])
    controls["product_surface_contact_three"] = recognize_columns(
        [0, 0], [1, 2], [[s, t], [s**2+t**2, s*t]])
    controls["rank_three_contact_two"] = recognize_columns(
        [0, 0, 0], [1, 2], [[s, t, 0], [s**2, s*t, t**2]])
    controls["two_support_points"] = recognize_columns(
        [0, 0], [1, 2], [[s, t], [s**2, s*t+t**2]])
    controls["no_coincidence_adjacent"] = recognize_columns(
        [0, 0, 0], [0, 1], [[1, 0, 0], [0, s, t]])
    controls["one_point_nonadjacent"] = recognize_columns(
        [0, 0], [0, 2], [[1, 0], [s**2, t**2]])
    expected = {
        "F1_contact_four": (4, 1, True),
        "product_surface_contact_three": (3, 1, True),
        "rank_three_contact_two": (2, 1, True),
        "two_support_points": (3, 2, False),
        "no_coincidence_adjacent": (0, 0, False),
        "one_point_nonadjacent": (2, 1, False),
    }
    for name, (length, points, affine) in expected.items():
        record = controls[name]
        check(name, (record["coincidence_length"], record["geometric_support_points"],
                     record["affine_space_by_imported_criterion"]) == (length, points, affine))
    # The first control's support is s=0, missed by the chart s=1.
    check("homogeneous test retains the point outside s=1 chart",
          controls["F1_contact_four"]["homogeneous_gcd"] == "s**4")
    rejected = {
        "same_divisor": must_reject("identical divisor rejected", [0, 0], [0, 0],
                                    [[1, 0], [2, 0]]),
        "vertical_component": must_reject("nonsaturated column rejected", [0, 0], [1, 2],
                                          [[s, 0], [s**2, t**2]]),
        "wrong_degree": must_reject("wrong twist rejected", [0, 1], [1, 2],
                                    [[s, t], [s*t, s**3+t**3]]),
    }
    # Explicit chart orders at the contact point s=0,t=1.
    # On q1!=0, the sections have fibre coordinate q2/q1=-s and -s/(s^3+1).
    difference = S.cancel(-s+s/(s**3+1))
    check("contact is exactly four in a regular local chart",
          S.cancel(difference/s**4).subs(s, 0) == -1)

    B, E, K = (b, 1), (e, 1), (a-2, -2)
    zero("B square", intersect(B, B)-(a+2*b))
    zero("E square", intersect(E, E)-(a+2*e))
    zero("B E", intersect(B, E)-(a+b+e))
    zero("adjunction B", intersect(K, B)+intersect(B, B)+2)
    zero("adjunction E", intersect(K, E)+intersect(E, E)+2)
    zero("canonical class after quotient by B", (a-2)-b*(-2)-(a+2*b-2))
    delta = S.symbols("delta")
    k = delta+2
    zero("canonical and different imply B E=2(delta+1)", k+delta-2*(delta+1))
    check("negative adjacency forces empty coincidence", (k+delta).subs(delta, -1) == 0)
    check("positive adjacency gives contact four", (k+delta).subs(delta, 1) == 4)
    zero("F1 canonical identity K=E-3B, f coefficient",
         (a-2-(e-3*b)).subs({a: 1, b: 1, e: 2}))
    check("F1 canonical identity h coefficient", -2 == 1-3)

    # a+2b=3 with a,b>=0 implies 0<=a<=3 and 0<=b<=1.
    candidates = [(aa, bb) for aa in range(4) for bb in range(2) if aa+2*bb == 3]
    check("complete normalized-section candidates", candidates == [(1, 1), (3, 0)])
    surviving = [(aa, bb) for aa, bb in candidates if aa == 0 or bb > 0]
    check("affineness removes the disjoint negative section", surviving == [(1, 1)])
    check("two omitted primes plus infinity exceed boundary rank", 2+1 > 2)

    def projective_class(d):
        return sum(L**i for i in range(d+1)) if d >= 0 else 0
    for rank in range(2, 9):
        inclusion_exclusion = ((1+L)*projective_class(rank-1)
                              -2*(1+L)*projective_class(rank-2)
                              +(1+L)*projective_class(rank-3)+r*L**(rank-2))
        zero("stratified class identity rank "+str(rank),
             inclusion_exclusion-(L**rank+(r-1)*L**(rank-2)))

    # Independent formal differentiation verifies the leading Jacobian coefficient.
    z, w = S.symbols("z w")
    p, q = S.Function("p")(w), S.Function("q")(w)
    wedge_controls = []
    for mp in range(1, 5):
        for mq in range(1, 5):
            P, Q = z**(-mp)*p, z**(-mq)*q
            wedge = S.diff(P, z)*S.diff(Q, w)-S.diff(P, w)*S.diff(Q, z)
            zero(f"leading differential coefficient {mp},{mq}",
                 z**(mp+mq+1)*wedge-(mq*S.diff(p, w)*q-mp*p*S.diff(q, w)))
            if mp+mq+1 == 3:
                wedge_controls.append((mp, mq))
    check("noncancelled pole three forces both pole orders one",
          wedge_controls == [(1, 1)])
    check("two pole-one coordinate divisors cannot give quartic degree", 3 < 4)
    check("no rational basepoint-free extension degree four",
          S.Poly(3*m*m-4, m, domain=S.QQ).ground_roots() == {})

    # At m=2 the two necessary budgets classify arithmetic patterns only.
    patterns = set()
    for twos in range(3):
        ones = 8-4*twos
        for marked_twos in range(twos+1):
            for marked_ones in range(ones+1):
                if 2*marked_twos+marked_ones == 6:
                    patterns.add((twos, ones, marked_twos, marked_ones))
    check("complete m2 arithmetic patterns",
          patterns == {(0, 8, 0, 6), (1, 4, 1, 4)})

    for epsilon in (0, 1):
        check(f"proper maximal base multiplicity gives forbidden form pole {epsilon}",
              epsilon-2 < 0)
    geometric_patterns = {x for x in patterns if x[0] == 0}
    check("m2 geometric refinement leaves only eight simple centers",
          geometric_patterns == {(0, 8, 0, 6)})

    report = dict(
        as_of="2026-11-14", outcome="PASS", check_count=len(CHECKS), checks=CHECKS,
        controls=controls, rejected_inputs=rejected,
        ruled_candidates_before_affineness=candidates,
        ruled_candidates_after_affineness=surviving,
        m2_basepoint_patterns=[
            dict(multiplicity_two_count=x[0], multiplicity_one_count=x[1],
                 on_strict_B_multiplicity_two=x[2], on_strict_B_multiplicity_one=x[3])
            for x in sorted(patterns)],
        m2_patterns_after_proper_multiplicity_obstruction=[
            dict(multiplicity_two_count=x[0], multiplicity_one_count=x[1],
                 on_strict_B_multiplicity_two=x[2], on_strict_B_multiplicity_one=x[3])
            for x in sorted(geometric_patterns)],
        assurance=dict(recognition_theorem="external import, not proved by checker",
                       geometry="written argument",
                       hypothetical_finite_map_constructed=False,
                       actual_normalization_ruled_compatibility_proved=False,
                       quartic_packet_exclusion=False,
                       formal_verification=False, independent_full_reproof=False),
        software=dict(python=platform.python_version(), sympy=S.__version__),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=round(time.monotonic()-start, 6))
    (HERE/"receipt.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"PASS {len(CHECKS)} exact controls; conditional F1 reduction; no quartic exclusion")


if __name__ == "__main__":
    main()

