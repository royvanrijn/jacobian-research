from sage.all import *
from pathlib import Path
import argparse
import json
import re


ap = argparse.ArgumentParser(
    description="Verify P1/P2/P3 equations, orientations, and exact Shioda height Gram over a finite field."
)
ap.add_argument("--targets", required=True)
ap.add_argument("--p3-dir", required=True)
args = ap.parse_args()

records = {}
for record in map(json.loads, Path(args.targets).read_text().splitlines()):
    records.setdefault((record["seed"], record["hit"]), []).append(record)
primes = {record["p"] for surface_records in records.values() for record in surface_records}
if len(primes) != 1:
    raise SystemExit(f"target file must contain exactly one prime, got {sorted(primes)}")
p = primes.pop()
p3_re = re.compile(
    r"MW3A10P3_HIT\|r=(?P<r>\d+)\|q3=(?P<q3>\d+)"
    r"\|q0=(?P<q0>\d+)\|q1=(?P<q1>\d+)\|q2=(?P<q2>\d+)"
    r"\|X=(?P<X>[\d,]+)\|Y=(?P<Y>[\d,]+)"
)

K = GF(p)
Kt = PolynomialRing(K, "t")
t = Kt.gen()
F = FractionField(Kt)
fiber_orders = [11, 3, 2, 2]
profiles = [
    [2, 1, 0, 1],
    [6, 2, 1, 1],
    [10, 2, 0, 1],
]
target_scaled = matrix(ZZ, [[79, 17, -1], [17, 106, 19], [-1, 19, 259]])


def parse_coefficients(text):
    return [int(value) for value in text.split(",")]


def intersection_with_zero(point):
    x = F(point[0])
    numerator = x.numerator()
    denominator = x.denominator()
    finite_poles = denominator.degree()
    infinity_pole = max(0, numerator.degree() - denominator.degree())
    value = QQ(finite_poles + infinity_pole - 4) / 2
    if value not in ZZ or value < 0:
        raise RuntimeError(f"invalid zero-section intersection {value} for x={x}")
    return ZZ(value)


def component_contribution(profile):
    return sum(QQ(i * (n - i)) / n for i, n in zip(profile, fiber_orders))


def section_height(point, profile):
    return QQ(4) + 2 * intersection_with_zero(point) - component_contribution(profile)


def profile_sum(left, right):
    return [(a + b) % n for a, b, n in zip(left, right, fiber_orders)]


def infinity_distance(point):
    """Return min(c,11-c) from the infinity degree of scaled y."""
    y = F(point[1])
    degree_at_infinity = y.numerator().degree() - y.denominator().degree()
    distance = 6 - degree_at_infinity
    if not 0 <= distance <= 5:
        raise RuntimeError(f"invalid I11 component distance {distance} for y={y}")
    return ZZ(distance)


def raw_odd_components(point, P1, kind):
    """Recover oriented I11/I3 labels using the calibrated P1=(2,1)."""
    translated = point + P1

    # At I3, c+1 is zero exactly when c=2.  The identity component maps to a
    # smooth point (or O), while the nonidentity components map to the node.
    x = F(translated[0])
    y = F(translated[1])
    if x.denominator()(0) == 0 or y.denominator()(0) == 0:
        c3 = 2
    else:
        at_node = K(x(0)) == K(3) and K(y(0)) == K(0)
        c3 = 1 if at_node else 2

    distance = infinity_distance(translated)
    if kind == 2:
        # c11=5 gives c11+2=7 (distance 4); c11=6 gives 8 (distance 3).
        if distance not in (3, 4):
            raise RuntimeError(f"unexpected translated P2 I11 distance {distance}")
        c11 = 6 if distance == 3 else 5
    elif kind == 3:
        # c11=1 gives 3 (distance 3); c11=10 gives 1 (distance 1).
        if distance not in (1, 3):
            raise RuntimeError(f"unexpected translated P3 I11 distance {distance}")
        c11 = 10 if distance == 1 else 1
    else:
        raise RuntimeError(f"unknown section kind {kind}")
    return c11, c3


def orient_to_target(point, P1, kind):
    desired = (6, 2) if kind == 2 else (10, 2)
    raw = raw_odd_components(point, P1, kind)
    if raw == desired:
        return point, raw, 1
    negative = ((11 - raw[0]) % 11, (3 - raw[1]) % 3)
    if negative == desired:
        return -point, raw, -1
    return None, raw, 0


verified = []
evaluated = 0
for log_path in sorted(Path(args.p3_dir).glob("seed*-hit*.p3.log")):
    name_match = re.fullmatch(r"seed(\d+)-hit(\d+)\.p3\.log", log_path.name)
    if not name_match:
        continue
    key = (int(name_match.group(1)), int(name_match.group(2)))
    surface_records = records.get(key)
    if surface_records is None:
        continue
    hits = [p3_re.fullmatch(line) for line in log_path.read_text().splitlines()]
    hits = [match for match in hits if match]
    if not hits:
        continue

    for p2_candidate, record in enumerate(surface_records, 1):
        A = Kt(record["A"])
        B = Kt(record["B"])
        E = EllipticCurve(F, [0, 0, 0, F(A), F(B)])
        P1 = E(F(Kt(record["X1"])), F(Kt(record["Y1"])))
        pole2 = K(record["P2"]["r"])
        z2 = t - pole2
        P2_raw = E(
            F(Kt(record["P2"]["X2"])) / z2**2,
            F(Kt(record["P2"]["Y2"])) / z2**3,
        )
        P2, P2_raw_components, P2_sign = orient_to_target(P2_raw, P1, 2)
        if P2 is None:
            print(
                "MW3A10TRIPLE_SURFACE|" + json.dumps(
                    {
                        "seed": key[0], "hit": key[1],
                        "P2_candidate": p2_candidate,
                        "target_P2_profile": False,
                        "P2_raw_components": [int(value) for value in P2_raw_components],
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
            continue

        for candidate_index, match in enumerate(hits, 1):
            evaluated += 1
            pole3 = K(int(match.group("r")))
            z3 = t - pole3
            P3_raw = E(
                F(Kt(parse_coefficients(match.group("X")))) / z3**2,
                F(Kt(parse_coefficients(match.group("Y")))) / z3**3,
            )
            P3, P3_raw_components, P3_sign = orient_to_target(P3_raw, P1, 3)
            if P3 is None:
                output = {
                    "seed": key[0], "hit": key[1],
                    "P2_candidate": p2_candidate, "candidate": candidate_index,
                    "P3_r": int(pole3), "target_P3_profile": False,
                    "P2_raw_components": [int(value) for value in P2_raw_components],
                    "P2_sign": int(P2_sign),
                    "P3_raw_components": [int(value) for value in P3_raw_components],
                }
                print("MW3A10TRIPLE|" + json.dumps(output, sort_keys=True), flush=True)
                continue
            points = [P1, P2, P3]

            zero_intersections = [intersection_with_zero(point) for point in points]
            self_heights = [section_height(point, profile) for point, profile in zip(points, profiles)]
            gram = matrix(QQ, 3)
            for i in range(3):
                gram[i, i] = self_heights[i]
                for j in range(i + 1, 3):
                    sum_point = points[i] + points[j]
                    sum_profile = profile_sum(profiles[i], profiles[j])
                    pairing = (
                        section_height(sum_point, sum_profile)
                        - self_heights[i]
                        - self_heights[j]
                    ) / 2
                    gram[i, j] = gram[j, i] = pairing

            scaled = (66 * gram).change_ring(ZZ)
            exact_target = scaled == target_scaled
            signed_target = any(
                diagonal_matrix(ZZ, signs) * scaled * diagonal_matrix(ZZ, signs) == target_scaled
                for signs in cartesian_product_iterator([[1, -1]] * 3)
            )
            determinant_ok = gram.det() == QQ(79) / 11
            output = {
                "seed": key[0],
                "hit": key[1],
                "P2_candidate": p2_candidate,
                "candidate": candidate_index,
                "P3_r": int(pole3),
                "P2_raw_components": [int(value) for value in P2_raw_components],
                "P2_sign": int(P2_sign),
                "P3_raw_components": [int(value) for value in P3_raw_components],
                "P3_sign": int(P3_sign),
                "zero_intersections": [int(value) for value in zero_intersections],
                "scaled_gram": [[int(value) for value in row] for row in scaled.rows()],
                "determinant": str(gram.det()),
                "exact_target": bool(exact_target),
                "signed_target": bool(signed_target),
                "determinant_ok": bool(determinant_ok),
            }
            print("MW3A10TRIPLE|" + json.dumps(output, sort_keys=True), flush=True)
            if signed_target and determinant_ok:
                verified.append(output)

print(
    f"MW3A10TRIPLE|done=1|verified={len(verified)}|evaluated={evaluated}|p3_candidates="
    + str(sum(1 for path in Path(args.p3_dir).glob("*.p3.log") for line in path.read_text().splitlines() if line.startswith("MW3A10P3_HIT|"))),
    flush=True,
)
