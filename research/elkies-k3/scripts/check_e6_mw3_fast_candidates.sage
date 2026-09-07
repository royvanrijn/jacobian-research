from sage.all import *
from pathlib import Path
import argparse


ap = argparse.ArgumentParser(
    description=(
        "Replay E6FASTCORE records, enforce exact P1/component opens, and "
        "exhaust the frame-glue P2/P3 search with exact height checks."
    )
)
ap.add_argument("--input", type=Path, nargs="+", required=True)
ap.add_argument("--max-print", type=int, default=20)
args = ap.parse_args()

def parse_records(paths):
    records = []
    seen = set()
    observed_primes = set()
    for path in paths:
        for line in path.read_text().splitlines():
            if line.startswith("E6FASTSUMMARY|"):
                summary = dict(
                    item.split("=", 1) for item in line.split("|")[1:]
                )
                if "p" in summary:
                    observed_primes.add(int(summary["p"]))
            if not line.startswith("E6FASTCORE|"):
                continue
            fields = {}
            for item in line.split("|")[1:]:
                key, value = item.split("=", 1)
                fields[key] = value
            fields["p"] = int(fields.get("p", 31))
            observed_primes.add(fields["p"])
            for key in ("r0", "s0", "x1", "a1", "a2", "a4", "s1"):
                fields[key] = int(fields[key])
            for key in ("roots", "A", "B", "X1", "Y1"):
                fields[key] = tuple(int(value) for value in fields[key].split(","))
            identity = (
                fields["p"], fields["A"], fields["B"],
                fields["X1"], fields["Y1"],
            )
            if identity not in seen:
                seen.add(identity)
                records.append(fields)
    return records, observed_primes


records, primes = parse_records(args.input)
if len(primes) != 1:
    raise SystemExit(f"inputs must contain exactly one field, got {sorted(primes)}")
p = primes.pop()
if not is_prime(p) or p in (2, 3, 79):
    raise SystemExit(f"unsupported field GF({p})")
K = GF(p)
T = PolynomialRing(K, "t")
t = T.gen()
F = FractionField(T)


def valuation_at(poly, point):
    factor = t-point
    valuation = 0
    while poly != 0 and poly(point) == 0:
        poly //= factor
        valuation += 1
    return valuation


def node_at(A, B, point):
    matches = [
        node for node in K
        if A(point) == -3*node**2
        and B(point) == 2*node**3
        and B.derivative(t)(point)+node*A.derivative(t)(point) == 0
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one node at {point}, got {matches}")
    return matches[0]


def meets_node(section, point, node):
    if section.is_zero():
        return False
    x_coordinate = F(section[0])
    y_coordinate = F(section[1])
    if (
        T(x_coordinate.denominator())(point) == 0
        or T(y_coordinate.denominator())(point) == 0
    ):
        return False
    return x_coordinate(point) == node and y_coordinate(point) == 0


def meets_ivstar_singular(section):
    if section.is_zero():
        return False
    x_coordinate = F(section[0])
    y_coordinate = F(section[1])
    x_order = (
        T(x_coordinate.numerator()).degree()
        - T(x_coordinate.denominator()).degree()
    )
    y_order = (
        T(y_coordinate.numerator()).degree()
        - T(y_coordinate.denominator()).degree()
    )
    return x_order < 4 and y_order < 6


def section_O_intersection(section):
    """Intersection with O on the degree-(4,6) K3 Weierstrass model."""
    if section.is_zero():
        raise ValueError("O.O is not used in the height calculation")
    x_coordinate = F(section[0])
    numerator = T(x_coordinate.numerator())
    denominator = T(x_coordinate.denominator())
    finite_order = denominator.degree()
    assert finite_order % 2 == 0
    infinity_order = numerator.degree()-denominator.degree()-4
    if infinity_order < 0:
        infinity_order = 0
    assert infinity_order % 2 == 0
    return ZZ((finite_order+infinity_order)//2)


def component_add(left, right):
    moduli = (3, 4, 4, 2, 2)
    return tuple((a+b) % modulus for a, b, modulus in zip(left, right, moduli))


def component_self_contribution(label):
    e6, i4a, i4b, i2a, i2b = label
    contribution = QQ(0) if e6 == 0 else QQ(4)/3
    for modulus, value in ((4, i4a), (4, i4b), (2, i2a), (2, i2b)):
        if value:
            contribution += QQ(value*(modulus-value))/modulus
    return contribution


def section_height(section, label):
    if section.is_zero():
        return QQ(0)
    return QQ(4) + 2*section_O_intersection(section) - component_self_contribution(label)


def height_gram(sections, labels):
    diagonal = [section_height(section, label) for section, label in zip(sections, labels)]
    result = matrix(QQ, len(sections), len(sections))
    for i in range(len(sections)):
        result[i, i] = diagonal[i]
        for j in range(i):
            sum_height = section_height(
                sections[i]+sections[j],
                component_add(labels[i], labels[j]),
            )
            result[i, j] = result[j, i] = (sum_height-diagonal[i]-diagonal[j])/2
    return result


def square_roots_by_coefficients(poly):
    if poly == 0:
        return [T(0)]
    degree = poly.degree()
    if degree % 2:
        return []
    root_degree = degree//2
    roots = []
    for leading in K(poly[degree]).sqrt(all=True):
        coefficients = [K(0)]*(root_degree+1)
        coefficients[root_degree] = leading
        for target_degree in range(degree-1, root_degree-1, -1):
            index = target_degree-root_degree
            known = K(0)
            for i in range(index+1, root_degree+1):
                j = target_degree-i
                if 0 <= j <= root_degree and j != index:
                    known += coefficients[i]*coefficients[j]
            coefficients[index] = (K(poly[target_degree])-known)/(2*leading)
        candidate = T(coefficients)
        if candidate**2 == poly:
            roots.append(candidate)
    return roots


def interpolate_x(points, nodes, pole):
    result = T(0)
    for index, (point, node) in enumerate(zip(points, nodes)):
        basis = T(1)
        denominator = K(1)
        for other_index, other in enumerate(points):
            if index == other_index:
                continue
            basis *= t-other
            denominator *= point-other
        result += node*(point-pole)**2*basis/denominator
    return result


def search_p2(A, B, curve, P1, points, nodes):
    fiber_product = prod(t-point for point in points)
    tested = square_hits = finite_hits = 0
    hits = []
    for pole in K:
        if pole in points:
            continue
        constant = interpolate_x(points, nodes, pole)
        z = t-pole
        for q0 in K:
            tested += 1
            X2 = constant+q0*fiber_product
            if X2(pole) == 0:
                continue
            H = X2**3+A*X2*z**4+B*z**6
            for Y2 in square_roots_by_coefficients(H):
                if any(Y2(point) != 0 for point in points):
                    continue
                square_hits += 1
                P2 = curve(F(X2)/F(z)**2, F(Y2)/F(z)**3)
                twice = 2*P2
                if not meets_node(twice, points[0], nodes[0]):
                    continue
                if meets_node(twice, points[1], nodes[1]):
                    continue
                finite_hits += 1
                if not meets_ivstar_singular(P2):
                    continue
                # P1 and target P2 have the same nonzero IV* component label.
                if meets_ivstar_singular(P1-P2):
                    continue
                hits.append((pole, q0, X2, Y2, P2))
    return tested, square_hits, finite_hits, hits


def search_p3(A, B, curve, P1, points, nodes):
    tested = square_hits = finite_hits = 0
    hits = []
    for x1 in K:
        for x2 in K:
            tested += 1
            X3 = nodes[0]+x1*t+x2*t**2
            H = X3**3+A*X3+B
            for Y3 in square_roots_by_coefficients(H):
                if Y3(0) != 0 or X3.degree() != 2 or Y3.degree() != 4:
                    continue
                square_hits += 1
                P3 = curve(F(X3), F(Y3))
                if any(
                    meets_node(P3, points[index], nodes[index])
                    for index in range(1, 4)
                ):
                    continue
                if not meets_node(2*P3, points[0], nodes[0]):
                    continue
                finite_hits += 1
                if not meets_ivstar_singular(P3):
                    continue
                # P1 and target P3 have opposite nonzero IV* component labels.
                if meets_ivstar_singular(P1+P3):
                    continue
                hits.append((x1, x2, X3, Y3, P3))
    return tested, square_hits, finite_hits, hits


def rank3_certificate(A, B, Delta, pole, sections, bad_parameters):
    remaining_masks = set(range(8))
    trace = []
    odd_fiber = None
    excluded = set(bad_parameters) | {pole}
    for parameter in K:
        if parameter in excluded or Delta(parameter) == 0:
            continue
        fiber = EllipticCurve(K, [A(parameter), B(parameter)])
        specialized = tuple(
            fiber(F(point[0])(parameter), F(point[1])(parameter))
            for point in sections
        )
        twice = {2*point for point in fiber}
        remaining_masks = {
            mask for mask in remaining_masks
            if sum(
                (specialized[index] for index in range(3) if (mask >> index) & 1),
                fiber(0),
            ) in twice
        }
        trace.append((int(parameter), fiber.cardinality(), sorted(remaining_masks)))
        if odd_fiber is None and fiber.cardinality() % 2:
            odd_fiber = (int(parameter), fiber.cardinality())
        if remaining_masks == {0} and odd_fiber is not None:
            break
    return remaining_masks == {0} and odd_fiber is not None, trace, odd_fiber


target_gram = matrix(QQ, [[23, -10, -8], [-10, 23, 1], [-8, 1, 23]])/12
label_p1 = (1, 0, 1, 0, 0)
label_p2 = (1, 1, 2, 1, 1)
label_p3_polynomial = (2, 3, 0, 0, 0)

valid_p1 = residual_rejections = p2_surfaces = p2_hits_total = 0
height2_rejections = relative_rejections = height3_rejections = 0
rank3_hits = 0
for index, record in enumerate(records, start=1):
    A, B = T(record["A"]), T(record["B"])
    X1, Y1 = T(record["X1"]), T(record["Y1"])
    lam, mu = map(K, record["roots"])
    Delta = -16*(4*A**3+27*B**2)
    fixed_discriminant = t**4*(t-1)**4*(t-lam)**2*(t-mu)**2
    if Delta.degree() != 16 or Delta % fixed_discriminant:
        raise RuntimeError(f"record {index}: fast multiplicity replay failed")
    residual = Delta//fixed_discriminant
    if (
        residual.degree() != 4
        or gcd(residual, residual.derivative()).degree() != 0
        or any(residual(point) == 0 for point in (K(0), K(1), lam, mu))
    ):
        residual_rejections += 1
        continue
    points = [K(0), K(1), lam, mu]
    nodes = [node_at(A, B, point) for point in points]
    curve = EllipticCurve(F, [F(A), F(B)])
    P1 = curve(F(X1), F(Y1))
    if not meets_ivstar_singular(P1):
        continue
    if meets_node(P1, points[0], nodes[0]):
        continue
    if not meets_node(P1, points[1], nodes[1]):
        continue
    if not meets_node(2*P1, points[1], nodes[1]):
        continue
    if any(meets_node(P1, points[i], nodes[i]) for i in (2, 3)):
        continue
    valid_p1 += 1

    p2_tested, p2_squares, p2_finite, p2_hits = search_p2(
        A, B, curve, P1, points, nodes
    )
    if p2_hits:
        p2_surfaces += 1
        p2_hits_total += len(p2_hits)
    print(
        f"E6CHECK|record={index}|r0={record['r0']}|s0={record['s0']}"
        f"|x1={record['x1']}|roots={record['roots'][0]},{record['roots'][1]}"
        f"|p2_tested={p2_tested}|p2_squares={p2_squares}"
        f"|p2_finite={p2_finite}|p2_hits={len(p2_hits)}",
        flush=True,
    )
    if not p2_hits:
        continue

    target_p2_hits = []
    for p2_index, (pole, q0, X2, Y2, P2) in enumerate(p2_hits, start=1):
        if p2_index <= args.max_print:
            print(
                f"E6CHECKP2_HIT|record={index}|r={int(pole)}|q0={int(q0)}"
                f"|X={','.join(str(int(c)) for c in X2.list())}"
                f"|Y={','.join(str(int(c)) for c in Y2.list())}",
                flush=True,
            )
        pair_gram = height_gram((P1, P2), (label_p1, label_p2))
        sum_O = section_O_intersection(P1+P2)
        assert (pair_gram == target_gram[:2, :2]) == (sum_O == 1)
        if pair_gram != target_gram[:2, :2]:
            height2_rejections += 1
            print(
                f"E6CHECKHEIGHT2_REJECT|record={index}|p2={p2_index}"
                f"|P1plusP2O={sum_O}"
                f"|gram={','.join(str(value) for value in pair_gram.list())}",
                flush=True,
            )
            continue
        target_p2_hits.append((p2_index, pole, q0, X2, Y2, P2))

    if not target_p2_hits:
        continue
    p3_tested, p3_squares, p3_finite, p3_hits = search_p3(
        A, B, curve, P1, points, nodes
    )
    print(
        f"E6CHECKP3|record={index}|tested={p3_tested}|squares={p3_squares}"
        f"|finite={p3_finite}|hits={len(p3_hits)}",
        flush=True,
    )
    for p2_index, pole, q0, X2, Y2, P2 in target_p2_hits:
        for p3_index, (x1, x2, X3, Y3, P3) in enumerate(p3_hits, start=1):
            # Target labels at I4(0) are P2=1 and P3=3 up to simultaneous
            # orientation, so their sum must lie on the identity component.
            if meets_node(P2+P3, points[0], nodes[0]):
                relative_rejections += 1
                continue
            gram = height_gram(
                (P1, P2, P3),
                (label_p1, label_p2, label_p3_polynomial),
            )
            if gram != target_gram:
                height3_rejections += 1
                continue
            # Specialization helper needs simplified affine coordinates.
            independent, trace, odd_fiber = rank3_certificate(
                A, B, Delta, pole,
                (P1, P2, P3),
                points,
            )
            if not independent:
                continue
            rank3_hits += 1
            print(
                f"E6CHECKRANK3_HIT|record={index}|p2={p2_index}|p3={p3_index}"
                f"|P3X={','.join(str(int(c)) for c in X3.list())}"
                f"|P3Y={','.join(str(int(c)) for c in Y3.list())}"
                f"|odd_fiber={odd_fiber[0]},{odd_fiber[1]}",
                flush=True,
            )

print(
    f"E6CHECKSUMMARY|p={p}|records={len(records)}"
    f"|residual_rejections={residual_rejections}"
    f"|valid_p1={valid_p1}|p2_surfaces={p2_surfaces}"
    f"|p2_hits={p2_hits_total}|height2_rejections={height2_rejections}"
    f"|relative_rejections={relative_rejections}"
    f"|height3_rejections={height3_rejections}"
    f"|rank3_hits={rank3_hits}",
    flush=True,
)
