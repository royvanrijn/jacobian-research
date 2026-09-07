#!/usr/bin/env sage -python
"""Certify the rational non-CM H3 anchor and the H21 entrance cubic.

The level-474 normalization identifies the published rational point

    (x,y) = (13/7, 12048/343)

with a rational point of the Elkies--Kumar H92 chart.  This checker then
matches that K3, through the four weighted Igusa invariants, with a rational
presentation of the H21 E7+E8 model.  The weighted scale is a rational
square, and the five short-Weierstrass coefficients are related by a
rational base scaling and rational quadratic twist.

The oriented Hilbert-cover coordinates from ``21/21.txt`` and ``92/92.txt``
are evaluated separately.  Both have square class ``-52203427``, so the two
orientations use the same quadratic field rather than a biquadratic one.
This is arithmetic data for the oriented Hilbert covers; it neither proves
nor obstructs rationality of the K3 Mordell--Weil section without an explicit
comparison of the markings.

Finally, the checker replays the entrance neighbor from ``21/21.txt``.  Its
generic fiber is a plane cubic over QQ(u); the displayed rational non-flex
point proves that the cubic is an elliptic curve over QQ(u).  This is the
``A2+A6+E8 -> E7+E8`` H21 entrance, not yet the desired
``E7+E8 -> E8+E6`` q=6 pencil with divisor ``O+(-P1)-F``.  Identifying and
executing that marked chord is the next equation-level step.
"""

from sage.all import *

import argparse
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_NORMALIZATION = (
    ROOT / "artifacts/generated-results/elkies-k3-h21-h92-level474-normalization.json"
)
DEFAULT_FACTOR = (
    ROOT / "artifacts/generated-results/elkies-k3-h21-h92-level474-factor-qq.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-h3-noncm-q6-source-anchor.json"
)
H92_SHA256 = "427559ecd4c2c19d0a4ed7df1019c8a351ed34f454691e9ef1080a8834e74ea1"
H21_ENTRANCE_SHA256 = "e0d0ea5ae18502fc0b51cf1999ab4e8b5755a40bffe2b7418ace6891d40a71a6"
H92_ENTRANCE_SHA256 = "5dde58046d9770fa78b514ae48509a238090ad4de7057b41e43ea308047101c2"
EK_SOURCE_URL = "https://export.arxiv.org/e-print/1209.3527"

PUBLISHED_X = QQ(13) / 7
PUBLISHED_Y = QQ(12048) / 343
EXPECTED_H92 = (QQ(-3621005) / 690947, QQ(158286) / 143585)
EXPECTED_H21 = (
    QQ(271621954946208883) / 51863976688786080,
    QQ(33855626850015548642165023481946578233811)
    / 37279290235251051531892004373889536000,
)
EXPECTED_CUBIC_POINT_X = -QQ(
    3233275090230032482678539458877352074438903840035914649417374578109347444492011378944062315708712660085121
) / QQ(
    150062187952117779106160007243308171967247482935454785199197675557527907483162593795348889600000000
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path):
    """Use a repository-relative path for repository-owned files."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def stage(name, **values):
    payload = "|".join(f"{key}={value}" for key, value in values.items())
    print(f"H3NONCMQ6|stage={name}" + (f"|{payload}" if payload else ""), flush=True)


def parse_h92(path):
    if digest(path) != H92_SHA256:
        raise ValueError(f"unexpected H92 ancillary hash for {path}")
    text = path.read_text()
    ring = PolynomialRing(QQ, names=("r", "s"))
    r, s = ring.gens()
    environment = {"r": r, "s": s}
    values = []
    for name in ("A1", "A", "B1", "B", "B2"):
        match = re.search(rf"\b{name}\s*=\s*(.*?);", text, flags=re.S)
        if match is None:
            raise ValueError(f"H92 source has no {name} assignment")
        expression = re.sub(
            r"\s+", " ", match.group(1).replace("^", "**")
        ).strip()
        values.append(ring(sage_eval(expression, locals=environment)))
    return ring, tuple(values)


def rational_function(expression, field):
    return field(expression.replace("^", "**"))


def normalization_point(data):
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    field = ring.fraction_field()

    published_x = rational_function(data["published_level_474"]["x"], field)
    roots = (published_x - PUBLISHED_X).numerator().roots(QQ)
    assert len(roots) == 1 and roots[0][1] == 1
    x0 = roots[0][0]

    t = rational_function(data["parameter"]["t"], field)(x0)
    a = rational_function(data["parameter"]["a"], field)(x0)
    y_formula = data["published_level_474"]["y"]
    assert y_formula.endswith("*Y")
    multiplier = rational_function(y_formula[:-2], field)(x0)
    Y = PUBLISHED_Y / multiplier

    r = (a + Y) / 2
    inverse_s = (Y - a) / 2
    s = 1 / inverse_s
    assert r / s == t
    assert (r, s) == EXPECTED_H92
    return x0, t, a, Y, r, s


def verify_level474_factor(data, r, s):
    value = sum(
        QQ(term["coefficient"]) * r ** term["r"] * s ** term["s"]
        for term in data["factor"]["coefficients"]
    )
    assert value == 0


def h21_coefficients(r, s):
    A1 = QQ(-1)
    A = -(s-r)**2 * (
        s**2 + (24*r**3+30*r**2-26*r-30)*s
        - 15*r**4-30*r**3+7*r**2+30*r+9
    ) / 3
    B1 = (
        r**4 - 2*s*r**3 + (s**2-QQ(5)/3)*r**2
        + QQ(10)/3*s*r - QQ(2)/3*s**2 + 2*s + 1
    )
    B = 2*(s-r)**4 * (
        s**2 + (-72*r**3-63*r**2+70*r+63)*s
        - 27*r**6-189*r**5-63*r**4+441*r**3
        + 280*r**2-252*r-189
    ) / 27
    B2 = (r-1)**6*(r+1)**4*(s-r)**6
    return A1, A, B1, B, B2


def igusa(coefficients):
    A1, A, B1, B, B2 = coefficients
    return (
        -24*B1/A1,
        -12*A,
        96*(A/A1)*B1-36*B,
        -4*A1*B2,
    )


def model_isomorphism(h21, h92):
    A1, A, B1, B, B2 = h21
    a1, a, b1, b, b2 = h92
    mu = (a/a1)/(A/A1)
    d2 = a1/(A1*mu**3)
    d3 = b1/(B1*mu**5)
    assert d2.is_square() and d3**2 == d2**3
    d = d2.sqrt()
    if d**3 != d3:
        d = -d
    # A rational coefficient twist would not by itself identify the two
    # marked elliptic K3 surfaces over QQ.  Here the twist parameter is
    # another rational square, so the coordinate scaling is already over QQ.
    assert d.is_square()
    assert (
        a1 == d**2*A1*mu**3
        and a == d**2*A*mu**4
        and b1 == d**3*B1*mu**5
        and b == d**3*B*mu**6
        and b2 == d**3*B2*mu**7
    )
    return mu, d


def h21_oriented_cover_value(r, s):
    return (
        16*s**4 - 8*r*(27*r**2-23)*s**3
        + (621*r**4-954*r**2+349)*s**2
        - 18*(r**3-r)*(33*r**2-29)*s
        + (r**2-1)*(189*(r**4-r**2)+16)
    )


def h92_oriented_cover_value(r, s):
    extra_i2 = (
        27*r**2*s**5 - r**4*s**4 + 42*r**3*s**4 - 30*r**2*s**4
        - 46*r*s**4 - r**5*s**3 + 17*r**4*s**3 - 69*r**3*s**3
        - 77*r**2*s**3 + 30*r*s**3 + 27*s**3 + 3*r**5*s**2
        - 28*r**4*s**2 - 4*r**3*s**2 + 69*r**2*s**2 + 42*r*s**2
        - 3*r**5*s + 9*r**4*s + 28*r**3*s + 17*r**2*s + r*s
        + r**5 + 3*r**4 + 3*r**3 + r**2
    )
    return -(s + r**2 + r) * (s**2 + r*s - r) * extra_i2


def h21_entrance_cubic(r, s):
    base = FunctionField(QQ, "u")
    u = base.gen()
    ring = PolynomialRing(base, names=("x", "t"))
    x, t = ring.gens()

    e = (r-s)**2*(r**2-1)
    univariate = PolynomialRing(QQ, "tt")
    tt = univariate.gen()
    c = e**2*(1-tt) + e**2*r**2*tt
    b = e*(1-tt) + e*r*s*tt
    modulus = tt**3-tt**2
    a = (b**2*c.inverse_mod(modulus)).mod(modulus)

    old_rhs = (
        x**3 + ring(a(t))*x**2
        + 2*ring(b(t))*(t**3-t**2)*x
        + ring(c(t))*(t**3-t**2)**2
    )
    linear_x = (
        ((2*r-2)*s-r**2+1)*t**2
        + (-2*r*s+r**2+1)*t - 2
    )
    constant = (
        (r**2-1)*(s-r)**2*(t-1)*t**2
        * ((r-1)**2*t**2+(1-r**2)*t-2)
    )
    y = ((t-1)**2*t**4*u - (linear_x*x+constant))/2
    raw = (old_rhs-y**2)(x=x*(t**3-t**2))
    denominator = (t**3-t**2)**3
    cubic = ring(raw/denominator)
    assert cubic.degree(x) == cubic.degree(t) == 3
    assert len(cubic.dict()) == 10

    projective_ring = PolynomialRing(base, names=("X", "T", "Z"))
    X, T, Z = projective_ring.gens()
    homogeneous = projective_ring(cubic(X/Z, T/Z)*Z**3)
    point = (base(EXPECTED_CUBIC_POINT_X), base(0), base(1))
    assert homogeneous(*point) == 0
    hessian = matrix([
        [homogeneous.derivative(v1, v2) for v1 in (X, T, Z)]
        for v2 in (X, T, Z)
    ]).det()
    hessian_value = hessian(*point)
    assert hessian_value != 0
    return cubic, hessian_value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h92", required=True, type=Path)
    parser.add_argument("--h92-entrance", required=True, type=Path)
    parser.add_argument("--h21-entrance", required=True, type=Path)
    parser.add_argument("--normalization", type=Path, default=DEFAULT_NORMALIZATION)
    parser.add_argument("--factor", type=Path, default=DEFAULT_FACTOR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    if digest(arguments.h21_entrance) != H21_ENTRANCE_SHA256:
        raise ValueError(
            f"unexpected H21 entrance ancillary hash for {arguments.h21_entrance}"
        )
    if digest(arguments.h92_entrance) != H92_ENTRANCE_SHA256:
        raise ValueError(
            f"unexpected H92 entrance ancillary hash for {arguments.h92_entrance}"
        )

    normalization = json.loads(arguments.normalization.read_text())
    factor = json.loads(arguments.factor.read_text())
    assert normalization["status"] == "PASS_LEVEL474_NORMALIZATION"
    assert factor["status"] == "PASS_CHARACTERISTIC_ZERO_FACTOR"

    x0, t, a, Y, r92, s92 = normalization_point(normalization)
    verify_level474_factor(factor, r92, s92)
    stage("h92_point", r=r92, s=s92, factor=0)

    h92_ring, h92_formulas = parse_h92(arguments.h92)
    unused_r, unused_s = h92_ring.gens()
    h92 = tuple(QQ(value(r92, s92)) for value in h92_formulas)
    r21, s21 = EXPECTED_H21
    h21 = h21_coefficients(r21, s21)
    I92 = igusa(h92)
    I21 = igusa(h21)
    scale = I21[0]/I92[0]
    weights = (1, 2, 3, 5)
    assert all(I21[index] == scale**weights[index]*I92[index] for index in range(4))
    assert scale.is_square()
    stage("h21_match", r=r21, s=s21, scale_square=1)

    mu, twist = model_isomorphism(h21, h92)
    stage("model_isomorphism", base_scale=mu, twist_square=1)

    h21_oriented_cover = h21_oriented_cover_value(r21, s21)
    h92_oriented_cover = h92_oriented_cover_value(r92, s92)
    oriented_square_class = ZZ(-52203427)
    h21_oriented_ratio = h21_oriented_cover / oriented_square_class
    h92_oriented_ratio = h92_oriented_cover / oriented_square_class
    assert not h21_oriented_cover.is_square() and h21_oriented_ratio.is_square()
    assert not h92_oriented_cover.is_square() and h92_oriented_ratio.is_square()
    oriented_root_ratio = h21_oriented_ratio.sqrt() / h92_oriented_ratio.sqrt()
    assert oriented_root_ratio in QQ
    stage(
        "oriented_covers",
        square_class=oriented_square_class,
        shared_field=1,
    )

    cubic, hessian_value = h21_entrance_cubic(r21, s21)
    stage(
        "h21_entrance_cubic",
        bidegree="3,3",
        terms=len(cubic.dict()),
        rational_point=1,
        nonflex=int(hessian_value != 0),
    )

    output = {
        "schema": "elkies-k3.h3-noncm-q6-source-anchor.v4",
        "status": "PASS_H3_NONCM_Q6_SOURCE_ANCHOR",
        "proof_boundary": (
            "This certifies the rational H92 and H21 presentations, their "
            "short-Weierstrass isomorphism over QQ, and a rational non-flex "
            "point on the H21 entrance plane cubic over QQ(u).  This entrance "
            "cubic is not yet identified with the E7+E8 to E8+E6 q=6 chord; "
            "the marked divisor-function step O+(-P1)-F remains open.  The "
            "matching nonsquare H21/H92 oriented Hilbert-cover values are "
            "recorded separately and are not interpreted as either a proof "
            "or an obstruction to rationality of the K3 section."
        ),
        "inputs": {
            "normalization": {
                "path": display_path(arguments.normalization),
                "sha256": digest(arguments.normalization),
            },
            "factor": {
                "path": display_path(arguments.factor),
                "sha256": digest(arguments.factor),
            },
            "h92": {
                "sha256": digest(arguments.h92),
                "source_url": EK_SOURCE_URL,
                "archive_member": "92/igusa92.txt",
            },
            "h21_entrance": {
                "sha256": digest(arguments.h21_entrance),
                "source_url": EK_SOURCE_URL,
                "archive_member": "21/21.txt",
            },
            "h92_entrance": {
                "sha256": digest(arguments.h92_entrance),
                "source_url": EK_SOURCE_URL,
                "archive_member": "92/92.txt",
            },
        },
        "published_point": {"x": str(PUBLISHED_X), "y": str(PUBLISHED_Y)},
        "normalization_point": {
            "x0": str(x0), "t": str(t), "a": str(a), "Y": str(Y),
        },
        "h92": {
            "r": str(r92), "s": str(s92),
            "coefficients": [str(value) for value in h92],
            "igusa": [str(value) for value in I92],
        },
        "h21": {
            "r": str(r21), "s": str(s21),
            "coefficients": [str(value) for value in h21],
            "igusa": [str(value) for value in I21],
        },
        "weighted_igusa_scale": {
            "scale": str(scale), "square_root": str(scale.sqrt()),
            "weights": list(weights),
        },
        "short_weierstrass_isomorphism": {
            "base_scale": str(mu),
            "quadratic_twist_scale": str(twist),
            "coordinate_scale": str(twist.sqrt()),
            "identity": (
                "H92=(d^2*(A1*mu^3,A*mu^4), "
                "d^3*(B1*mu^5,B*mu^6,B2*mu^7))"
            ),
        },
        "oriented_hilbert_covers": {
            "square_class": str(oriented_square_class),
            "shared_quadratic_field": "QQ(sqrt(-52203427))",
            "h21": {
                "value": str(h21_oriented_cover),
                "square_ratio_root": str(h21_oriented_ratio.sqrt()),
                "is_rational_square": False,
            },
            "h92": {
                "value": str(h92_oriented_cover),
                "square_ratio_root": str(h92_oriented_ratio.sqrt()),
                "is_rational_square": False,
            },
            "h21_root_over_h92_root": str(oriented_root_ratio),
            "interpretation": (
                "The two oriented Hilbert-cover coordinates use the same "
                "quadratic field.  This is not, by itself, the field of "
                "definition of the K3 MW section."
            ),
        },
        "h21_entrance_cubic": {
            "base": "QQ(u)", "bidegree_x_t": [3, 3],
            "terms": len(cubic.dict()),
            "rational_point": [str(EXPECTED_CUBIC_POINT_X), "0", "1"],
            "hessian_at_point": str(hessian_value),
            "nonflex": True,
        },
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    stage("artifact", path=display_path(arguments.output), sha256=digest(arguments.output))
    stage("complete", status=output["status"])


if __name__ == "__main__":
    main()
