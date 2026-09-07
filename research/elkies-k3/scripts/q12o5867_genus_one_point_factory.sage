#!/usr/bin/env sage-python
"""Exact birational point maps for the q12/orbit5867 neighbour.

The parent is the P1229-pointed ``4I2+16I1`` / MW13 model and the child is
the published compact rootless R17 model.  On the common open chart this
script exposes both maps

    (parent base s, parent point P) <-> (published base t, R17 point Q).

The forward map first evaluates the exact chord pencil, then uses the
nonbranch point at s=0 on the resulting binary quartic, and finally applies
the certified raw-q12 to published-R17 coordinate change.  The inverse map
uses the literal inverse formulas.  The default ``controls`` mode transports
the public complement bases at the four rank-25--28 fibres and t=3/8 back to
the parent and verifies exact round trips.

This is a rational map on an affine open set.  Points at the zero section,
the pointed quartic origin, and zeros of displayed denominators require the
corresponding projective chart and are rejected explicitly here.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, gcd
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json"
PARENT_MODEL = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
HORIZONTAL = LOCAL / "q12o5867-target-horizontal-qq.json"
Q12_MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
MATCH = ROOT / "artifacts/generated-results/elkies-k3-h3-q12o5867-elkies-2026-coordinate-match.json"
PUBLISHED_MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
CONTROL_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_specialization_controls_v1.json"
)
VISIBILITY_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_visibility_record_curves_v1.json"
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def qtext(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def point_record(point) -> list[str]:
    return [qtext(point[0]), qtext(point[1])]


def primitive_pair(value) -> tuple[int, int]:
    value = QQ(value)
    numerator = int(value.numerator())
    denominator = int(value.denominator())
    common = gcd(abs(numerator), denominator)
    return numerator // common, denominator // common


def evaluate_polynomial(poly, value):
    answer = value.parent().zero() if hasattr(value, "parent") else QQ.zero()
    for coefficient in reversed(poly.list()):
        answer = answer * value + coefficient
    return answer


class Q12O5867PointFactory:
    """Load and compose the exact q12/orbit5867 point maps."""

    def __init__(self):
        self.parent_record = json.loads(PARENT_MODEL.read_text())
        self.horizontal_record = json.loads(HORIZONTAL.read_text())
        self.q12_record = json.loads(Q12_MODEL.read_text())
        self.match_record = json.loads(MATCH.read_text())
        self.published_record = json.loads(PUBLISHED_MODEL.read_text())
        if self.parent_record["status"] != "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO":
            raise ArithmeticError("the exact pointed 4A1 parent is not pinned")
        if self.horizontal_record["status"] != "PASS_EXACT_QQ_Q12O5867_TARGET_HORIZONTAL_SECTION":
            raise ArithmeticError("the exact q12 horizontal is not pinned")
        if self.q12_record["status"] != "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN":
            raise ArithmeticError("the exact q12 child is not pinned")
        if self.match_record["status"] != "PASS_EXACT_Q12O5867_IS_ELKIES_2026_PUBLISHED_MODEL":
            raise ArithmeticError("the compact published chart is not pinned")

        self.Rs = PolynomialRing(QQ, "s")
        self.svar = self.Rs.gen()
        self.Ks = self.Rs.fraction_field()
        child = self.parent_record["child"]
        self.parent_A = self.Rs(child["minimal_A_coefficients_low_to_high"])
        self.parent_B = self.Rs(child["minimal_B_coefficients_low_to_high"])
        reducible = child["finite_reducible_fibres"]
        if [row["factor"] for row in reducible][-1:] != ["u"]:
            raise ArithmeticError("the parent I2 support ordering changed")
        nonzero_supports = [
            QQ(row["factor"].split("u - ", 1)[1])
            for row in reducible
            if row["factor"].startswith("u - ")
        ]
        if len(nonzero_supports) != 2:
            raise ArithmeticError("the parent no longer has two finite nonzero I2 supports")
        self.parent_i2_r1, self.parent_i2_r2 = nonzero_supports

        target = self.horizontal_record["target"]
        self.HX = self.Rs(target["x_numerator_coefficients_low_to_high"])
        dx = self.Rs(target["x_denominator_coefficients_low_to_high"])
        self.HY = self.Rs(target["y_numerator_coefficients_low_to_high"])
        dy = self.Rs(target["y_denominator_coefficients_low_to_high"])
        self.HZ = self.Rs.one()
        for factor, exponent in dx.factor():
            if int(exponent) % 2:
                raise ArithmeticError("the horizontal x denominator is not a square")
            self.HZ *= factor.monic() ** (int(exponent) // 2)
        if self.HZ**2 != dx or self.HZ**3 != dy:
            raise ArithmeticError("the horizontal projective denominator changed")
        if self.HY**2 != self.HX**3 + self.parent_A * self.HX * self.HZ**4 + self.parent_B * self.HZ**6:
            raise ArithmeticError("the q12 horizontal misses the parent surface")

        pairs = []
        for item in self.q12_record["smooth_RR"]["basis_pairs"]:
            pairs.append(
                (
                    self.Rs(item["AA_coefficients_low_to_high"]),
                    self.Rs(item["BB_coefficients_low_to_high"]),
                )
            )
        if len(pairs) != 2:
            raise ArithmeticError("the q12 pencil no longer has a two-plane")
        (self.AA0, self.BB0), (self.AA1, self.BB1) = pairs

        self.Ru = PolynomialRing(QQ, "u")
        self.uvar = self.Ru.gen()
        self.Ku = self.Ru.fraction_field()
        self.Rus = PolynomialRing(self.Ru, "s")
        self.quartic = self.Rus(
            [
                self.Ru(row)
                for row in self.q12_record["binary_quartic"][
                    "coefficients_in_old_v_low_to_high"
                ]
            ]
        )
        if self.quartic.degree() != 4:
            raise ArithmeticError("the q12 residual is not a binary quartic")
        self.square_factor = self._reconstruct_square_factor()
        self._prepare_pointed_quartic()
        self._prepare_compact_chart()

    def _lift_parent(self, poly):
        return self.Rus([self.Ru(value) for value in self.Rs(poly).list()])

    def _reconstruct_square_factor(self):
        aa = self._lift_parent(self.AA0) + self.uvar * self._lift_parent(self.AA1)
        bb = self._lift_parent(self.BB0) + self.uvar * self._lift_parent(self.BB1)
        X, Y, Z, A = map(
            self._lift_parent, (self.HX, self.HY, self.HZ, self.parent_A)
        )
        raw = (
            aa**4
            - 6 * X * aa**2 * bb**2
            + 8 * Y * aa * bb**3
            - 3 * X**2 * bb**4
            - 4 * A * bb**4 * Z**4
        )
        after_collision, remainder = raw.quo_rem(Z**4)
        if remainder:
            raise ArithmeticError("the q12 chord radicand lost Z^4 divisibility")
        square_quotient, remainder = after_collision.quo_rem(self.quartic)
        if remainder:
            raise ArithmeticError("the stored q12 quartic no longer divides the radicand")
        factorization = square_quotient.factor()
        if any(int(exponent) % 2 for _, exponent in factorization):
            raise ArithmeticError("the q12 residual quotient is not a square")
        unit = self.Ru(factorization.unit())
        if not unit.is_square():
            raise ArithmeticError("the q12 residual square has a nonsquare unit")
        square = self.Rus(unit.sqrt())
        for factor, exponent in factorization:
            square *= factor ** (int(exponent) // 2)
        if after_collision != self.quartic * square**2:
            raise ArithmeticError("the reconstructed q12 square factor failed replay")
        return square

    def _prepare_pointed_quartic(self):
        e, d, c, b, a = [self.Ku(value) for value in self.quartic.list()]
        if not self.Ru(e).is_square():
            raise ArithmeticError("the stored s=0 quartic point is not rational")
        q0 = self.Ku(self.Ru(e).sqrt())
        self.q0 = q0
        self.qd = d
        self.qc = c
        self.qa1 = d / q0
        self.qa2 = c - d**2 / (4 * q0**2)
        self.qa3 = 2 * q0 * b
        self.qb2 = self.qa1**2 + 4 * self.qa2
        qa4 = -4 * q0**2 * a
        qa6 = self.qa2 * qa4
        qb4 = 2 * qa4 + self.qa1 * self.qa3
        qb6 = self.qa3**2 + 4 * qa6
        c4 = self.qb2**2 - 24 * qb4
        c6 = -self.qb2**3 + 36 * self.qb2 * qb4 - 216 * qb6
        A_pointed = -c4 / 48
        B_pointed = -c6 / 864
        A_raw = self.Ru(self.q12_record["child"]["minimal_A_coefficients_low_to_high"])
        B_raw = self.Ru(self.q12_record["child"]["minimal_B_coefficients_low_to_high"])
        if 81 * A_pointed != self.Ku(A_raw) or 729 * B_pointed != self.Ku(B_raw):
            raise ArithmeticError("the pointed quartic scaling to the raw child changed")
        self.raw_A = A_raw
        self.raw_B = B_raw

    def _prepare_compact_chart(self):
        base = self.match_record["base_change"]
        self.mobius_a = QQ(base["a"])
        self.mobius_b = QQ(base["b"])
        self.mobius_c = QQ(base["c"])
        self.mobius_d = QQ(base["d"])
        self.compact_scale = QQ(self.match_record["weierstrass_isomorphism"]["s"])
        if self.mobius_a * self.mobius_d - self.mobius_b * self.mobius_c == 0:
            raise ArithmeticError("the raw-to-published base map is singular")
        self.published_A = self.Ru(self.published_record["A_coefficients_low_to_high"])
        self.published_B = self.Ru(self.published_record["B_coefficients_low_to_high"])

    def _eval_parent(self, poly, s):
        return evaluate_polynomial(self.Rs(poly), self.Ks(s))

    def _eval_bivariate(self, poly, s, u):
        s = u.parent()(s)
        answer = u.parent().zero()
        for coefficient in reversed(self.Rus(poly).list()):
            answer = answer * s + evaluate_polynomial(self.Ru(coefficient), u)
        return answer

    def parent_curve(self, s):
        s = QQ(s)
        return EllipticCurve(
            QQ,
            [
                0,
                0,
                0,
                QQ(self.parent_A(s)),
                QQ(self.parent_B(s)),
            ],
        )

    def parent_crossratio_base(self, s):
        """Map the finite I2 supports (0,r1,r2) to (0,1,infinity)."""

        s = QQ(s)
        r1, r2 = self.parent_i2_r1, self.parent_i2_r2
        if s == r2:
            raise ZeroDivisionError("the normalized parent base is infinite")
        return s * (r1 - r2) / (r1 * (s - r2))

    def raw_to_published_base(self, u):
        u = QQ(u)
        denominator = self.mobius_c * u + self.mobius_d
        if denominator == 0:
            raise ZeroDivisionError("the compact published base is infinite")
        return (self.mobius_a * u + self.mobius_b) / denominator

    def published_to_raw_base(self, t):
        t = QQ(t)
        denominator = self.mobius_c * t - self.mobius_a
        if denominator == 0:
            raise ZeroDivisionError("the raw q12 base is infinite")
        return (self.mobius_b - self.mobius_d * t) / denominator

    def raw_to_published_point(self, u, point):
        u = QQ(u)
        x, y = map(QQ, point)
        denominator = self.mobius_c * u + self.mobius_d
        scale_x = self.compact_scale**2 * denominator**4
        scale_y = self.compact_scale**3 * denominator**6
        if scale_x == 0 or scale_y == 0:
            raise ZeroDivisionError("the compact Weierstrass chart is singular")
        result = (x / scale_x, y / scale_y)
        t = self.raw_to_published_base(u)
        if result[1] ** 2 != result[0] ** 3 + self.published_A(t) * result[0] + self.published_B(t):
            raise ArithmeticError("the raw-to-published point missed the compact curve")
        return t, result

    def published_to_raw_point(self, t, point):
        t = QQ(t)
        u = self.published_to_raw_base(t)
        x, y = map(QQ, point)
        denominator = self.mobius_c * u + self.mobius_d
        result = (
            self.compact_scale**2 * denominator**4 * x,
            self.compact_scale**3 * denominator**6 * y,
        )
        if result[1] ** 2 != result[0] ** 3 + self.raw_A(u) * result[0] + self.raw_B(u):
            raise ArithmeticError("the published-to-raw point missed the q12 child")
        return u, result

    def parent_to_raw(self, s, point):
        s = QQ(s)
        px, py = map(QQ, point)
        A = QQ(self.parent_A(s))
        B = QQ(self.parent_B(s))
        if py**2 != px**3 + A * px + B:
            raise ArithmeticError("the input point misses the 4A1 parent fibre")
        Z = QQ(self.HZ(s))
        if Z == 0:
            raise ZeroDivisionError("the horizontal section needs another parent chart")
        hx = QQ(self.HX(s)) / Z**2
        hy = QQ(self.HY(s)) / Z**3
        if px == hx:
            raise ZeroDivisionError("the input lies on a vertical chord denominator")
        slope = (py + hy) / (px - hx)
        restrictions = [
            QQ(AA(s)) + QQ(BB(s)) * Z * slope
            for AA, BB in ((self.AA0, self.BB0), (self.AA1, self.BB1))
        ]
        if restrictions[1] == 0:
            raise ZeroDivisionError("the raw child parameter is infinite")
        u = -restrictions[0] / restrictions[1]
        aa = QQ(self.AA0(s)) + u * QQ(self.AA1(s))
        bb = QQ(self.BB0(s)) + u * QQ(self.BB1(s))
        square = QQ(self._eval_bivariate(self.square_factor, s, u))
        if bb == 0 or square == 0:
            raise ZeroDivisionError("the q12 chord needs another quartic chart")
        W = bb**2 * (2 * px + hx - slope**2) / square
        if W**2 != QQ(self._eval_bivariate(self.quartic, s, u)):
            raise ArithmeticError("the parent point missed the q12 binary quartic")
        if s == 0:
            raise ZeroDivisionError("the pointed quartic origin needs its projective chart")
        q0 = QQ(self.q0(u))
        d = QQ(self.qd(u))
        c = QQ(self.qc(u))
        a1 = QQ(self.qa1(u))
        a3 = QQ(self.qa3(u))
        b2 = QQ(self.qb2(u))
        xg = (2 * q0 * (W + q0) + d * s) / s**2
        yg = (
            4 * q0**2 * (W + q0)
            + 2 * q0 * (d * s + c * s**2)
            - d**2 * s**2 / (2 * q0)
        ) / s**3
        raw = (
            9 * (xg + b2 / 12),
            27 * (yg + (a1 * xg + a3) / 2),
        )
        if raw[1] ** 2 != raw[0] ** 3 + self.raw_A(u) * raw[0] + self.raw_B(u):
            raise ArithmeticError("the pointed-quartic image missed the raw q12 child")
        return u, raw, W

    def raw_to_parent(self, u, point):
        u = QQ(u)
        X, Y = map(QQ, point)
        if Y**2 != X**3 + self.raw_A(u) * X + self.raw_B(u):
            raise ArithmeticError("the input point misses the raw q12 child")
        q0 = QQ(self.q0(u))
        d = QQ(self.qd(u))
        a1 = QQ(self.qa1(u))
        a2 = QQ(self.qa2(u))
        a3 = QQ(self.qa3(u))
        b2 = QQ(self.qb2(u))
        xg = X / 9 - b2 / 12
        yg = Y / 27 - (a1 * xg + a3) / 2
        if yg == 0:
            raise ZeroDivisionError("the pointed-quartic inverse needs another chart")
        s = 2 * q0 * (xg + a2) / yg
        if s == 0:
            raise ZeroDivisionError("the pointed-quartic inverse reached its origin")
        W = (xg * s**2 - d * s) / (2 * q0) - q0
        if W**2 != QQ(self._eval_bivariate(self.quartic, s, u)):
            raise ArithmeticError("the raw child inverse missed the q12 quartic")
        Z = QQ(self.HZ(s))
        if Z == 0:
            raise ZeroDivisionError("the inverse chord needs another parent chart")
        AA = QQ(self.AA0(s)) + u * QQ(self.AA1(s))
        BB = QQ(self.BB0(s)) + u * QQ(self.BB1(s))
        if BB == 0:
            raise ZeroDivisionError("the inverse q12 chord has BB=0")
        slope = -AA / (BB * Z)
        square = QQ(self._eval_bivariate(self.square_factor, s, u))
        hx = QQ(self.HX(s)) / Z**2
        hy = QQ(self.HY(s)) / Z**3
        px = (W * square / BB**2 - hx + slope**2) / 2
        py = slope * (px - hx) - hy
        if py**2 != px**3 + self.parent_A(s) * px + self.parent_B(s):
            raise ArithmeticError("the q12 chord inverse missed the parent fibre")
        restrictions = [
            QQ(AAi(s)) + QQ(BBi(s)) * Z * slope
            for AAi, BBi in ((self.AA0, self.BB0), (self.AA1, self.BB1))
        ]
        if restrictions[1] == 0 or u != -restrictions[0] / restrictions[1]:
            raise ArithmeticError("the reconstructed parent point has the wrong pencil value")
        return s, (px, py), W

    def parent_to_published(self, s, point):
        u, raw, W = self.parent_to_raw(s, point)
        t, published = self.raw_to_published_point(u, raw)
        return t, published, {"raw_u": u, "quartic_W": W, "raw_point": raw}

    def published_to_parent(self, t, point):
        u, raw = self.published_to_raw_point(t, point)
        s, parent, W = self.raw_to_parent(u, raw)
        return s, parent, {"raw_u": u, "quartic_W": W, "raw_point": raw}


def load_control_points():
    cas = ROOT / "elliptic-curves/cas"
    ecsearch = ROOT / "elliptic-curves"
    sys.path[:0] = [str(ecsearch), str(cas)]
    from elkies_rank25 import POINTS as rank25
    from elkies_rank26 import POINTS as rank26
    from elkies_rank27 import POINTS as rank27
    from elkies_rank28 import POINTS as rank28
    from verify_icarm_curve394_rank21 import PUBLIC_POINTS as curve394

    return {
        "-2/377": rank25,
        "-308/251": rank26,
        "2456/135": rank27,
        "-9529/5471": rank28,
        "3/8": curve394,
    }


def compact_minimal_to_affine_points(parameter, points):
    """Undo the exact public global-minimal change and dehomogenize."""

    cas = ROOT / "elliptic-curves/cas"
    ecsearch = ROOT / "elliptic-curves"
    sys.path[:0] = [str(ecsearch), str(cas)]
    from ecsearch.q12o5867_specialization import (
        evaluate_projective_specialization,
        global_minimal_model_with_change,
        load_q12o5867_data,
    )
    from elliptic_candidate_record import target_point_to_source

    data = load_q12o5867_data(
        PUBLISHED_MODEL,
        ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json",
    )
    numerator, denominator = primitive_pair(parameter)
    specialization = evaluate_projective_specialization(data, numerator, denominator)
    minimal_model, change, _metadata = global_minimal_model_with_change(
        specialization.model
    )
    affine_points = []
    for point in points:
        projective = target_point_to_source(point, change)
        projective = tuple(QQ(value) for value in projective)
        affine = (
            projective[0] / QQ(denominator) ** 4,
            projective[1] / QQ(denominator) ** 6,
        )
        affine_points.append(affine)
    return tuple(QQ(value) for value in minimal_model), change.to_record(), affine_points


def controls_payload(factory: Q12O5867PointFactory):
    controls = json.loads(CONTROL_CERTIFICATE.read_text())
    visibility = json.loads(VISIBILITY_CERTIFICATE.read_text())
    invisible_labels = {
        row["parameter"]: set(row["canonical_complement_labels"])
        for row in visibility["visibility"]
    }
    public_points = load_control_points()
    rows = []
    for fibre in controls["fibres"]:
        parameter = fibre["parameter"]
        indices = [
            int(index) - 1
            for index in fibre["public_complement"]["source_point_indices_one_based"]
        ]
        selected = [public_points[parameter][index] for index in indices]
        minimal_model, minimal_change, affine_points = compact_minimal_to_affine_points(
            QQ(parameter), selected
        )
        records = []
        for label_index, (source_index, minimal_point, affine_point) in enumerate(
            zip(indices, selected, affine_points), start=1
        ):
            parent_base, parent_point, intermediate = factory.published_to_parent(
                QQ(parameter), affine_point
            )
            roundtrip_t, roundtrip_point, forward = factory.parent_to_published(
                parent_base, parent_point
            )
            if roundtrip_t != QQ(parameter) or roundtrip_point != affine_point:
                raise ArithmeticError("a public control failed the exact birational round trip")
            curve = factory.parent_curve(parent_base)
            P = curve(parent_point)
            height = P.height(precision=128)
            normalized_parent_base = factory.parent_crossratio_base(parent_base)
            complement_label = f"public-complement-Q{label_index}"
            records.append(
                {
                    "public_complement_label": complement_label,
                    "source_public_point_index_one_based": source_index + 1,
                    "public_minimal_point": point_record(minimal_point),
                    "published_affine_point": point_record(affine_point),
                    "parent_base": qtext(parent_base),
                    "parent_base_projective_height": max(
                        abs(int(parent_base.numerator())),
                        int(parent_base.denominator()),
                    ),
                    "parent_crossratio_base": qtext(normalized_parent_base),
                    "parent_crossratio_base_projective_height": max(
                        abs(int(normalized_parent_base.numerator())),
                        int(normalized_parent_base.denominator()),
                    ),
                    "parent_point": point_record(parent_point),
                    "parent_point_canonical_height_128bit": str(height),
                    "parent_point_is_torsion": bool(P.has_finite_order()),
                    "rational_bisection_atlas_invisible": (
                        complement_label in invisible_labels[parameter]
                    ),
                    "raw_q12_parameter": qtext(intermediate["raw_u"]),
                    "quartic_ordinate": qtext(intermediate["quartic_W"]),
                    "exact_parent_and_child_equation_checks": True,
                    "exact_forward_inverse_roundtrip": True,
                    "specialized_MW13_representation": {
                        "status": "NOT_YET_COMPUTED",
                        "boundary": (
                            "The exact parent point and canonical height are recorded. "
                            "Membership and coordinates in a specialized thirteen-section "
                            "basis require a separately pinned equation-level MW13 basis."
                        ),
                    },
                }
            )
        rows.append(
            {
                "label": fibre["label"],
                "published_parameter": parameter,
                "minimal_model": [qtext(value) for value in minimal_model],
                "minimal_change_source_to_public": [str(value) for value in minimal_change],
                "public_complement_dimension": len(indices),
                "points": records,
            }
        )
    return {
        "schema": "elkies-k3.q12o5867-genus-one-point-factory-controls.v1",
        "status": "PASS_EXACT_Q12O5867_BIRATIONAL_POINT_MAP_AND_CONTROL_ROUNDTRIPS",
        "inputs": {
            str(path.relative_to(ROOT)): file_sha256(path)
            for path in (
                PARENT_MODEL,
                HORIZONTAL,
                Q12_MODEL,
                MATCH,
                PUBLISHED_MODEL,
                CONTROL_CERTIFICATE,
                VISIBILITY_CERTIFICATE,
            )
        },
        "map": {
            "domain": "P1229-pointed 4A1/MW13 parent",
            "codomain": "published compact rootless R17 fibration",
            "forward": "(s,x_parent,y_parent) -> (t,x_R17,y_R17)",
            "inverse": "(t,x_R17,y_R17) -> (s,x_parent,y_parent)",
            "normalized_parent_base": (
                "z=s*(r1-r2)/(r1*(s-r2)), mapping the ordered finite parent "
                "I2 supports (0,r1,r2) to (0,1,infinity)"
            ),
            "parent_i2_supports_0_r1_r2": [
                "0",
                qtext(factory.parent_i2_r1),
                qtext(factory.parent_i2_r2),
            ],
            "chart_boundary": (
                "Affine common open only: zero section, pointed origin, infinite base, "
                "and displayed denominator-zero loci require projective companion charts."
            ),
        },
        "controls": rows,
        "rank28_focus": {
            "published_parameter": "-9529/5471",
            "public_complement_dimension": 11,
            "rational_bisection_visible_dimension": 1,
            "rational_bisection_invisible_dimension": 10,
            "invisible_labels": sorted(invisible_labels["-9529/5471"]),
            "outcome": (
                "All ten invisible directions transport exactly, but none has a "
                "small parent cross-ratio base or small canonical height in this "
                "calibration; this does not prove absence of a simpler MW13 word."
            ),
        },
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/q12o5867_genus_one_point_factory.sage --mode controls "
            "--output artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json"
        ),
        "proof_boundary": (
            "The birational formulas, all 42 public-complement transports, parent/child "
            "equation identities, and forward/inverse round trips are exact. Canonical "
            "heights are 128-bit numerical evaluations. Specialized MW13 coordinates "
            "and positive-rank fibre enumeration are not yet certified by this artifact."
        ),
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("controls", "forward", "inverse"), default="controls")
    parser.add_argument("--output", type=Path, default=GENERATED)
    parser.add_argument("--parent-base")
    parser.add_argument("--parent-x")
    parser.add_argument("--parent-y")
    parser.add_argument("--published-t")
    parser.add_argument("--published-x")
    parser.add_argument("--published-y")
    return parser.parse_args()


def main():
    args = parse_args()
    factory = Q12O5867PointFactory()
    if args.mode == "controls":
        payload = controls_payload(factory)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        point_count = sum(len(row["points"]) for row in payload["controls"])
        print(
            "Q12O5867POINTFACTORY|controls={}|roundtrips=exact|status={}|output={}".format(
                point_count, payload["status"], args.output
            )
        )
        return
    if args.mode == "forward":
        if None in (args.parent_base, args.parent_x, args.parent_y):
            raise ValueError("forward mode requires --parent-base, --parent-x, --parent-y")
        t, point, intermediate = factory.parent_to_published(
            QQ(args.parent_base), (QQ(args.parent_x), QQ(args.parent_y))
        )
        payload = {
            "published_t": qtext(t),
            "published_affine_point": point_record(point),
            "raw_q12_parameter": qtext(intermediate["raw_u"]),
            "raw_q12_point": point_record(intermediate["raw_point"]),
            "quartic_ordinate": qtext(intermediate["quartic_W"]),
        }
    else:
        if None in (args.published_t, args.published_x, args.published_y):
            raise ValueError("inverse mode requires --published-t, --published-x, --published-y")
        s, point, intermediate = factory.published_to_parent(
            QQ(args.published_t), (QQ(args.published_x), QQ(args.published_y))
        )
        payload = {
            "parent_base": qtext(s),
            "parent_point": point_record(point),
            "raw_q12_parameter": qtext(intermediate["raw_u"]),
            "raw_q12_point": point_record(intermediate["raw_point"]),
            "quartic_ordinate": qtext(intermediate["quartic_W"]),
        }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
