
#!/usr/bin/env sage-python

from sage.all import *
from pathlib import Path
import argparse
import re
import subprocess
import sys
import time

sys.path.insert(0, "elliptic-curves/cas")

from anonymous_rank_candidate import (
    GENERAL_WEIERSTRASS_COEFFICIENTS,
)

from elkies_klagsbrun_rank29 import (
    COEFFICIENT_B as E29_B,
    PUBLISHED_POINTS as E29_POINTS,
)

PROTOCOL = "ANONPOINT"

A1, A2, A3, A4, A6 = [
    QQ(v) for v in GENERAL_WEIERSTRASS_COEFFICIENTS
]

E = EllipticCurve(QQ, GENERAL_WEIERSTRASS_COEFFICIENTS)


def gp_int(x):
    return str(ZZ(x))


def gp_rat(x):
    x = QQ(x)
    if x.denominator() == 1:
        return str(x.numerator())
    return f"({x.numerator()}/{x.denominator()})"


def point_ok(x, y):
    return (
        y*y + x*y
        == x**3 + A4*x + A6
    )


def x_cubic_value(x):
    # z = 2y+x
    #
    # z^2 = 4x^3 + x^2 + 4A4*x + 4A6
    return 4*x**3 + x**2 + 4*A4*x + 4*A6


def transformed_polynomial(center, scale):
    """
    Return coefficients of

        F(center + scale*u)

    where F(x)=(2y+x)^2.
    """
    R = PolynomialRing(QQ, "u")
    u = R.gen()
    x = QQ(center) + QQ(scale)*u
    f = 4*x**3 + x**2 + 4*A4*x + 4*A6
    return R(f)


def polynomial_gp(f):
    terms = []

    for degree in range(f.degree() + 1):
        c = QQ(f[degree])
        if c == 0:
            continue
        terms.append(f"({gp_rat(c)})*x^{degree}")

    return "+".join(terms) if terms else "0"


POINT_RE = re.compile(r"^PT\|([^|]+)\|([^|]+)$", re.MULTILINE)


def run_chart(chart_id, center, scale, height, timeout):
    f = transformed_polynomial(center, scale)

    program = f"""
default(parisizemax, 32G);
P={polynomial_gp(f)};
gettime();
R=hyperellratpoints(P,{height});
ms=gettime();
print("TIME|",ms);
for(i=1,#R,print("PT|",R[i][1],"|",R[i][2]));
quit;
"""

    started = time.monotonic()

    try:
        proc = subprocess.run(
            ["gp", "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
        )

    except subprocess.TimeoutExpired:
        print(
            f"{PROTOCOL}|stage=chart"
            f"|id={chart_id}"
            f"|height={height}"
            f"|status=TIMEOUT"
            f"|seconds={timeout}",
            flush=True,
        )
        return []

    wall = time.monotonic() - started

    if proc.returncode != 0:
        err = proc.stderr.strip().replace("\n", " ")[:500]
        print(
            f"{PROTOCOL}|stage=chart"
            f"|id={chart_id}"
            f"|status=ERROR"
            f"|stderr={err}",
            flush=True,
        )
        return []

    raw = POINT_RE.findall(proc.stdout)

    points = []

    for u_text, z_text in raw:
        u = QQ(u_text)
        z = QQ(z_text)

        x = QQ(center) + QQ(scale)*u

        for zsign in (z, -z):
            y = (-x + zsign) / 2

            if not point_ok(x, y):
                raise AssertionError(
                    f"off-curve point from chart {chart_id}"
                )

            points.append((x, y))

    print(
        f"{PROTOCOL}|stage=chart"
        f"|id={chart_id}"
        f"|center={center}"
        f"|scale={scale}"
        f"|height={height}"
        f"|raw={len(raw)}"
        f"|images={len(points)}"
        f"|wall={wall:.3f}",
        flush=True,
    )

    return points


def characteristic_x_scale():
    RF = RealField(256)
    return ZZ(round(RF(abs(A6)) ** (RF(1)/3)))


def real_root_centers():
    RF = RealField(256)
    R = PolynomialRing(RF, "x")
    x = R.gen()

    f = 4*x**3 + x**2 + 4*RF(A4)*x + 4*RF(A6)

    roots = f.roots()

    result = []

    for root, _ in roots:
        if abs(root.imag()) < RF("1e-50"):
            result.append(ZZ(round(root.real())))

    return sorted(set(result))


def e29_scaled_centers(limit):
    """
    Purely heuristic.

    Match characteristic B^(1/3) scales between E29 and this curve and
    transport a subset of the known E29 x-coordinate locations.
    """

    RF = RealField(256)

    scale = (
        RF(abs(A6))
        / RF(abs(E29_B))
    ) ** (RF(1)/3)

    candidates = []

    for index, point in enumerate(E29_POINTS):
        x29 = QQ(point[0])

        center = ZZ(round(RF(x29) * scale))

        candidates.append(
            (
                index + 1,
                center,
                abs(center),
            )
        )

    # Prefer moderate absolute centers first.
    candidates.sort(key=lambda x: (x[2], x[0]))

    return [
        (index, center)
        for index, center, _ in candidates[:limit]
    ]


def canonical_key(P):
    return (QQ(P[0]), QQ(P[1]))


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--heights",
        default="1000,10000,100000",
        help="comma-separated hyperellratpoints heights",
    )

    ap.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="seconds per individual PARI chart",
    )

    ap.add_argument(
        "--e29-centers",
        type=int,
        default=12,
        help="number of rescaled E29 centers",
    )

    ap.add_argument(
        "--scales",
        default="1,1000000,1000000000000,1000000000000000000",
        help="affine x scales",
    )

    ap.add_argument(
        "--out",
        default="artifacts/local/elliptic-curves/anonymous-initial-points.txt",
    )

    args = ap.parse_args()

    heights = [
        int(x)
        for x in args.heights.split(",")
        if x.strip()
    ]

    scales = [
        ZZ(x)
        for x in args.scales.split(",")
        if x.strip()
    ]

    print(
        f"{PROTOCOL}|stage=start"
        f"|heights={heights}"
        f"|timeout={args.timeout}"
        f"|scales={scales}",
        flush=True,
    )

    xscale = characteristic_x_scale()

    print(
        f"{PROTOCOL}|stage=scale"
        f"|B_cuberoot={xscale}",
        flush=True,
    )

    charts = []

    # Direct chart: useful for unexpectedly small rational points.
    charts.append(("direct", ZZ(0), ZZ(1)))

    # Real roots of the x-cubic are geometrically distinguished locations.
    for i, center in enumerate(real_root_centers()):
        for scale in scales:
            charts.append(
                (
                    f"realroot{i+1}_s{scale}",
                    center,
                    scale,
                )
            )

    # Characteristic +/- B^(1/3) centers.
    for sign in (-1, 1):
        center = sign * xscale

        for scale in scales:
            charts.append(
                (
                    f"bscale_{sign:+d}_s{scale}",
                    center,
                    scale,
                )
            )

    # Rescaled E29 generator locations.
    for index, center in e29_scaled_centers(args.e29_centers):
        for scale in scales:
            charts.append(
                (
                    f"e29p{index:02d}_s{scale}",
                    center,
                    scale,
                )
            )

    # Deduplicate equal (center, scale) charts.
    unique = []
    seen_chart = set()

    for chart in charts:
        key = (chart[1], chart[2])

        if key in seen_chart:
            continue

        seen_chart.add(key)
        unique.append(chart)

    charts = unique

    print(
        f"{PROTOCOL}|stage=charts"
        f"|count={len(charts)}",
        flush=True,
    )

    all_points = {}

    for height in heights:
        for chart_id, center, scale in charts:
            points = run_chart(
                chart_id,
                center,
                scale,
                height,
                args.timeout,
            )

            for P in points:
                key = canonical_key(P)

                # +/- have same x but are genuinely distinct points.
                if key not in all_points:
                    all_points[key] = {
                        "point": P,
                        "chart": chart_id,
                        "height": height,
                    }

                    print(
                        f"{PROTOCOL}|stage=NEW"
                        f"|chart={chart_id}"
                        f"|height={height}"
                        f"|x={P[0]}"
                        f"|y={P[1]}",
                        flush=True,
                    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as h:
        for item in all_points.values():
            P = item["point"]

            h.write(
                f"{P[0]}\t{P[1]}"
                f"\t{item['chart']}"
                f"\t{item['height']}\n"
            )

    print(
        f"{PROTOCOL}|stage=done"
        f"|distinct_points={len(all_points)}"
        f"|out={out}",
        flush=True,
    )


if __name__ == "__main__":
    main()
