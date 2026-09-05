"""Exact, retrospective observability of known points in pointed search boxes.

Consumes the existing PointedQuarticSearch.chart_record/search transcript.
No point enumeration, Sage, class group, numerical height or rank inference.
Keep oracle points out of prospective centre/coordinate selection. A negative
answer concerns this point in this chart, not its entire Mordell--Weil coset.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, gcd, isqrt, lcm
import json


def rational(value):
    if isinstance(value, (float, bool)):
        raise ValueError("exact rational input required")
    return Q(value)


def point(value):
    if value is None:
        return None
    values = (value["x"], value["y"]) if isinstance(value, dict) else value
    if len(values) != 2:
        raise ValueError("expected an affine point or null for O")
    return tuple(map(rational, values))


def multiply(m, n):
    a, b, c, d = m
    e, f, g, h = n
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def matrix(values):
    if len(values) != 4:
        raise ValueError("four matrix entries required")
    a, b, c, d = map(rational, values)
    if not a*d-b*c:
        raise ValueError("singular horizontal matrix")
    return a, b, c, d


def transform(coefficients, m):
    a, b, c, d = m
    out = [Q(0)]*5
    for i, f in enumerate(coefficients):
        for j in range(i+1):
            for k in range(5-i):
                out[j+k] += f*comb(i, j)*a**j*b**(i-j)*comb(4-i, k)*c**k*d**(4-i-k)
    return tuple(out)


def primitive(n, d):
    n, d = rational(n), rational(d)
    if not n and not d:
        raise ValueError("zero projective coordinate")
    scale = lcm(n.denominator, d.denominator)
    n, d = int(n*scale), int(d*scale)
    g = gcd(n, d)
    n, d = n//g, d//g
    return (-n, -d) if d < 0 or (d == 0 and n < 0) else (n, d)


def prepare_chart(record):
    """Recompute the original-short-model transport and quartic identity."""
    curve = tuple(map(rational, record["input"]["curve"]))
    if len(curve) == 2:
        curve = (Q(0), Q(0), Q(0), *curve)
    if len(curve) != 5:
        raise ValueError("two or five Weierstrass coefficients required")
    a1, a2, a3, a4, a6 = curve
    c2, c4, c6 = a2+a1*a1/4, a4+a1*a3/2, a6+a3*a3/4
    shift = c2/3
    A, B = c4-c2*c2/3, c6-c2*c4/3+2*c2**3/27
    if not 4*A**3+27*B**2:
        raise ValueError("singular input curve")
    if tuple(map(rational, record["short_model"])) != (0, 0, 0, A, B) or rational(record["short_model_x_shift"]) != shift:
        raise ArithmeticError("wrong short-model transport")
    centre = point(record["base_point"])
    if centre is None:
        raise ValueError("finite pointed centre required")
    x0, y0 = centre
    xq, yq = x0+shift, y0+(a1*x0+a3)/2
    if yq*yq != xq**3+A*xq+B:
        raise ArithmeticError("centre is off the curve")
    raw = (-3*xq*xq-4*A, -8*yq, -6*xq, Q(0), Q(1))
    saved = record["pointed_chart"]
    scale = rational(saved["curve_coordinate_scale"])
    den = rational(saved["point_denominator_root"])
    k = rational(saved["shift_mod_denominator_squared"])
    ordinate = rational(record["ordinate_scale"])
    if not scale or not ordinate or den <= 0 or den.denominator != 1:
        raise ValueError("invalid chart scales")
    base_matrix = matrix(saved["unimodular_horizontal_matrix"])
    horizontal = matrix(record["horizontal_matrix"])
    # Original short-model slope t = ((den/scale)*upper+k/(den*scale)*lower)/lower.
    raw_matrix = multiply((den/scale, k/(den*scale), Q(0), Q(1)), multiply(base_matrix, horizontal))
    coefficients = tuple(map(rational, record["coefficients"]))
    if len(coefficients) != 5 or any(f.denominator != 1 for f in coefficients):
        raise ValueError("five integral final quartic coefficients required")
    factor = (den/(scale*scale*ordinate))**2
    if transform(raw, raw_matrix) != tuple(factor*f for f in coefficients):
        raise ArithmeticError("horizontal/ordinate transport does not represent the declared curve")
    return curve, (A, B), centre, (xq, yq), shift, raw_matrix, coefficients


def point_visibility(record, known_point):
    """Locate a supplied point, then separate geometry, coverage and output.

    'VISIBLE_NOT_RECORDED' means the supplied point is inside declared completed
    coverage but missing from finite_curve_points. It is an audit discrepancy,
    not automatically a proven bug: worker completeness still has its original
    trust boundary. Missing transcripts return OBSERVABLE_WITHOUT_TRANSCRIPT.
    """
    curve, (A, B), centre, (xq, yq), shift, m, coefficients = prepare_chart(record)
    p = point(known_point)
    a1, _, a3, _, _ = curve
    if p is not None:
        x, y = p[0]+shift, p[1]+(a1*p[0]+a3)/2
        if y*y != x**3+A*x+B:
            raise ValueError("oracle point is off the input curve")
    if p is None or p == centre:
        return {"status": "KNOWN_POINTED_ENDPOINT", "minimum_affine_height": None,
                "claim_boundary": "O and the supplied centre are not new-point discoveries."}
    if x == xq:
        if not yq:  # Here -Q=Q, already handled above.
            raise ArithmeticError("unexpected point above a two-torsion x-coordinate")
        t = -(3*xq*xq+A)/(2*yq)  # Removable 0/0 at R=-Q: tangent at -Q.
    else:
        t = (y+yq)/(x-xq)
    a, b, c, d = m
    n, denominator = primitive(d*t-b, -c*t+a)
    f = sum(coefficient*n**i*denominator**(4-i) for i, coefficient in enumerate(coefficients))
    if f.denominator != 1 or f < 0 or isqrt(f.numerator)**2 != f:
        raise ArithmeticError("rational oracle does not lift to an integral square in the declared chart")
    result = {"coordinate": [str(n), str(denominator)],
              "minimum_affine_height": max(abs(n), denominator) if denominator else None,
              "at_parameter_infinity": not denominator,
              "square_root_absolute": str(isqrt(f.numerator)),
              "status": "OBSERVABLE_WITHOUT_TRANSCRIPT",
              "claim_boundary": "One known point, one chart; no rank, coset-completeness or absence claim."}
    if "height_bound" not in record:
        return result
    h = record["height_bound"]
    first, last = record.get("denominator_start", 1), record.get("denominator_end", h)
    if any(type(v) is not int for v in (h, first, last)) or not 1 <= first <= last <= h:
        raise ValueError("invalid declared height interval")
    requested = (abs(n) <= h and first <= denominator <= last) if denominator else True
    done = record.get("completed_denominator")
    if done is not None and (type(done) is not int or not first-1 <= done <= last):
        raise ValueError("invalid completed interval")
    completed = (requested and done is not None and denominator <= done) if denominator else bool(record.get("infinity_checked", False))
    returned = any(point(q) == p for q in record.get("finite_curve_points", ()))
    hit = any(tuple(map(int, q[:2])) == (n, denominator) for q in record.get("primitive_square_hits", ()))
    if returned and not completed:
        raise ArithmeticError("recorded oracle point lies outside declared completed coverage")
    status = ("VISIBLE_AND_RECORDED" if returned else "VISIBLE_NOT_RECORDED") if completed else (
        "UNSEARCHED_INTERVAL" if requested else "OUTSIDE_BOX")
    result.update(status=status, in_requested_box=requested, in_completed_box=completed,
                  recorded_square_hit=hit, recorded_curve_point=returned)
    return result


def masked_control(curve, points, gram, withheld_indices):
    """Split an independently certified subgroup into detector input and oracle.

    The principal Gram block is valid on the retained subgroup. Independence
    is inherited from the input certificate, NOT established by this helper.
    Do not include the oracle object in a search request or use its points to
    choose centres. Recovering withheld generic points tests exposure, not
    exceptional-rank incidence. No per-curve target rank enters the input.
    """
    points = [point(p) for p in points]
    if not points or any(p is None for p in points) or len(set(points)) != len(points):
        raise ValueError("distinct finite basis points required")
    model = tuple(map(rational, curve))
    if len(model) == 2:
        model = (Q(0), Q(0), Q(0), *model)
    if len(model) != 5:
        raise ValueError("two or five curve coefficients required")
    a1, a2, a3, a4, a6 = model
    if any(y*y+a1*x*y+a3*y != x**3+a2*x*x+a4*x+a6 for x, y in points):
        raise ValueError("a supplied control point is off the curve")
    n = len(points)
    if len(gram) != n or any(len(row) != n for row in gram):
        raise ValueError("Gram/basis dimensions differ")
    g = tuple(tuple(map(rational, row)) for row in gram)
    if any(g[i][j] != g[j][i] for i in range(n) for j in range(n)):
        raise ValueError("Gram must be symmetric")
    # Exact Schur pivots reject indefinite/degenerate scoring Grams. This
    # validates the metric, not that it is the curve's canonical height form.
    reduced = [list(row) for row in g]
    for k in range(n):
        pivot = reduced[k][k]
        if pivot <= 0:
            raise ValueError("positive definite Gram required")
        for i in range(k+1, n):
            for j in range(k+1, n):
                reduced[i][j] -= reduced[i][k]*reduced[k][j]/pivot
    indices = list(withheld_indices)
    if any(type(i) is not int or not 0 <= i < n for i in indices) or len(set(indices)) != len(indices) or not 0 < len(indices) < n:
        raise ValueError("withhold a nonempty proper set of basis indices")
    hidden = set(indices)
    kept = [i for i in range(n) if i not in hidden]
    serialize = lambda p: dict(zip(("x", "y"), map(str, p)))
    search_input = {"curve": list(map(str, map(rational, curve))),
                    "points": [serialize(points[i]) for i in kept],
                    "metric_gram": [[str(g[i][j]) for j in kept] for i in kept]}
    binding = sha256(json.dumps(search_input, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    oracle = {"input_sha256": binding, "withheld_points": [serialize(points[i]) for i in sorted(hidden)],
              "endpoint": "WITHHELD_KNOWN_DIRECTIONS_NOT_NEW_RANK",
              "assumption": "Original points have a separate independence certificate."}
    return search_input, oracle
