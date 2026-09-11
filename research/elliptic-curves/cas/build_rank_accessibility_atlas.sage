#!/usr/bin/env sage-python
"""Rebuild the Curve302/11952 finite-atlas accessibility comparison.

This is a retrospective measurement program, not a point search.  It fixes a
target-blind chart bank before looking at any target and serializes the exact
inputs for every chart and every parameter evaluation.  ``D`` always means the
displayed independent point lattice; it is not asserted to be all of E(Q).

The main payload is gzip-compressed because every exact parameter numerator,
denominator, residual cancellation, and chart normalisation is retained.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import gzip
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from itertools import combinations, product
import json
from pathlib import Path
import sys

from sage.all import GF, QQ, RR, ZZ, RealField, EllipticCurve, matrix, pari, vector
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
OUT = ART / "rank_accessibility_atlas_curve302_11952_v1.json.gz"
SUMMARY = ART / "rank_accessibility_atlas_curve302_11952_v1.summary.json"
PARENT302 = ART / "curve302_recovered_mw17_parent_v1.json"
CERT302 = ART / "icarm_curve302_rank31_v1.json.gz"
LIVE11952 = ROOT / "artifacts" / "local" / "elliptic-curves" / "broad-rank-v1" / "runtime" / "research" / "broad-cases" / "b-7bb187bc9254e81c6283"
SEED11952 = LIVE11952 / "batch-000/generic/b-7bb187bc9254e81c6283/seed-M17.json"
TERM11952 = LIVE11952 / "batch-002/search-00/terminal.json"
VERIFY11952 = LIVE11952 / "batch-002/search-00/verified.json"

sys.path.insert(0, str(CAS))
import icarm_curve302 as curve302  # noqa: E402
import half_lattice_pointed_sieve as pointed  # noqa: E402

geometry = SourceFileLoader(
    "rank_accessibility_geometry", str(CAS / "prospective_half_lattice.sage")
).load_module()


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical(value) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()


def qtext(value) -> str:
    return str(QQ(value))


def fraction(value):
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def point_record(point) -> list[str]:
    return [qtext(point[0]), qtext(point[1])]


def word_key(word) -> str:
    return ",".join(map(str, word))


def add_words(*words):
    return tuple(sum(entries) for entries in zip(*words))


def neg_word(word):
    return tuple(-entry for entry in word)


def unit_word(dimension: int, index: int):
    return tuple(int(position == index) for position in range(dimension))


def binary_rank(rows, width: int) -> int:
    return int(matrix(GF(2), rows, ncols=width).rank())


def log2_exact(integer) -> str:
    integer = abs(int(integer))
    if integer < 1:
        raise ValueError("height input must be positive")
    return "%.70f" % RealField(256)(integer).log(2)


def primitive_integral_x_map(n_coefficients, f_coefficients):
    """One fixed primitive integral pair (N,F), before evaluating a target."""
    from math import gcd, lcm
    values = tuple(n_coefficients) + tuple(f_coefficients)
    denominator = lcm(*(int(value.denominator()) for value in values))
    integers = [int(value * denominator) for value in values]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    if not content:
        raise ArithmeticError("zero x-map pair")
    return tuple(value // content for value in integers[:5]), tuple(value // content for value in integers[5:]), denominator, content


def evaluate_homogeneous(coefficients, numerator: int, denominator: int) -> int:
    return sum(int(coefficient) * numerator**i * denominator ** (4 - i) for i, coefficient in enumerate(coefficients))


def short_vector_summary(form):
    """Certified-for-the-rounded-form short-vector census through an LLL bound.

    Canonical heights are numerical, so the word ``successive'' below is not a
    claim of exact real-arithmetic certification.  We retain a high-resolution
    integral rounding and the complete fplll enumeration through the maximum
    LLL-basis norm; this is enough to give reproducible numerical minima and
    concrete short quotient words without pretending that they are a theorem
    about an exact height pairing.
    """
    from fpylll import Enumeration, EvaluatorStrategy, GSO, IntegerMatrix
    scale = ZZ(10) ** 12
    rounded = matrix(ZZ, [[ZZ((entry * scale).round()) for entry in row] for row in form.rows()])
    U = matrix(ZZ, pari(rounded).qflllgram()).transpose()
    if abs(U.det()) != 1:
        raise ArithmeticError("LLL quotient change of basis is not unimodular")
    reduced = U * rounded * U.transpose()
    gso = GSO.Mat(IntegerMatrix.from_matrix([list(map(int, row)) for row in reduced.rows()]), gram=True, float_type="dd", update=True)
    cap = 100000
    enumeration = Enumeration(gso, strategy=EvaluatorStrategy.BEST_N_SOLUTIONS, nr_solutions=cap)
    solutions = enumeration.enumerate(0, form.nrows(), float(max(reduced.diagonal())) + 1.0, 0)
    if len(solutions) >= cap:
        raise ArithmeticError("short-vector census hit its declared cap")
    rows = []
    for rounded_norm, coordinates in solutions:
        reduced_word = vector(ZZ, [ZZ(round(value)) for value in coordinates])
        word = reduced_word * U
        if not any(word):
            continue
        height = vector(RealField(384), word) * form * vector(RealField(384), word)
        rows.append({"word_in_exceptional_quotient_basis": list(map(int, word)), "height": decimal_text(height), "rounded_norm": str(ZZ(round(rounded_norm)))})
    rows.sort(key=lambda row: RealField(384)(row["height"]))
    independent, minima = [], []
    for row in rows:
        trial = matrix(ZZ, independent + [row["word_in_exceptional_quotient_basis"]])
        if trial.rank() > len(independent):
            independent.append(row["word_in_exceptional_quotient_basis"])
            minima.append(row)
    if len(minima) != form.nrows():
        raise ArithmeticError("LLL enumeration did not span the quotient")
    return {
        "method": "LLL then complete fplll enumeration through the maximum rounded LLL-basis norm",
        "rounded_height_scale": str(scale),
        "rounded_enumeration_solution_count": len(rows),
        "successive_minima_numerical": minima,
        "short_vectors_numerical": rows[: min(32, len(rows))],
        "boundary": "These are reproducible numerical canonical-height minima for the stated rounded form, not exact symbolic Neron--Tate heights.",
    }


def decimal_text(value) -> str:
    return format(value, "f")


def standard_bank(generators):
    """Signed units and signed pair sums/differences: exactly 2*r^2 words."""
    words = []
    for source in generators:
        words.extend((source, neg_word(source)))
    for left, right in combinations(generators, 2):
        words.extend((add_words(left, right), add_words(left, neg_word(right)),
                      add_words(neg_word(left), right), add_words(neg_word(left), neg_word(right))))
    if len(words) != 2 * len(generators) ** 2 or len(set(words)) != len(words):
        raise ArithmeticError("native target-blind anchor bank is malformed")
    return tuple(words)


def generic_control_bank(generic, count: int):
    """Fixed lexicographic signed generic shells, with the requested cardinality.

    The first 578 charts are exactly the generic signed unit/pair bank.  Later
    charts append support-three shells.  Selection uses only generic words and
    the declared count, never a target or an exceptional point.
    """
    dimension = len(generic)
    representatives = []
    for support_size in range(1, dimension + 1):
        for support in combinations(range(dimension), support_size):
            for signs in product((-1, 1), repeat=support_size - 1):
                coefficients = [0] * dimension
                coefficients[support[0]] = 1
                for index, sign in zip(support[1:], signs):
                    coefficients[index] = sign
                word = tuple(sum(coefficients[i] * generic[i][j] for i in range(dimension))
                             for j in range(len(generic[0])))
                representatives.append(word)
                if 2 * len(representatives) == count:
                    answer = tuple(item for rep in representatives for item in (rep, neg_word(rep)))
                    if len(set(answer)) != count:
                        raise ArithmeticError("generic-only control has duplicate words")
                    return answer
    raise ArithmeticError("generic-only shell stream exhausted")


class Atlas:
    def __init__(self, label, model, displayed, generic_words, exceptional_words, bindings, rank_claim):
        self.label = label
        self.model = tuple(QQ(item) for item in model)
        self.E = EllipticCurve(QQ, self.model)
        self.displayed = tuple(self.E(point) for point in displayed)
        self.dimension = len(self.displayed)
        self.generic_words = tuple(tuple(map(int, word)) for word in generic_words)
        self.exceptional_words = tuple(tuple(map(int, word)) for word in exceptional_words)
        self.bindings = bindings
        self.rank_claim = rank_claim
        self._point_cache = {}
        self._chart_cache = {}
        self._norm_cache = {}
        self.gram, self.height_asymmetry = geometry.canonical_height_gram(
            self.model, [tuple(map(QQ, point.xy())) for point in self.displayed]
        )
        self.H = matrix(RealField(384), [[RealField(384)(str(value)) for value in row] for row in self.gram])

    def point(self, word):
        word = tuple(map(int, word))
        key = word_key(word)
        if key not in self._point_cache:
            answer = self.E(0)
            for coefficient, source in zip(word, self.displayed):
                if coefficient:
                    answer += coefficient * source
            self._point_cache[key] = answer
        return self._point_cache[key]

    def chart(self, word):
        word = tuple(map(int, word))
        key = word_key(word)
        if key not in self._chart_cache:
            anchor = self.point(word)
            if anchor.is_zero():
                raise ArithmeticError("anchor at infinity")
            a, b = map(QQ, anchor.xy())
            A, B = self.model[3:]
            raw_f = (-3 * a * a - 4 * A, -8 * b, -6 * a, QQ(0), QQ(1))
            raw_n = (a**3 + 4 * B, 4 * a * b, 6 * a * a + 4 * A, 4 * b, a)
            integral_n, integral_f, denominator, content = primitive_integral_x_map(raw_n, raw_f)
            chart = pointed.make_chart(tuple(map(fraction, self.model)), (fraction(a), fraction(b))).record()
            self._chart_cache[key] = {
                "word_in_D": list(word),
                "anchor": point_record(anchor),
                "raw_F_Q_ascending": [qtext(value) for value in raw_f],
                "raw_N_Q_ascending": [qtext(value) for value in raw_n],
                "primitive_integral_x_map": {
                    "N_ascending": list(map(str, integral_n)), "F_ascending": list(map(str, integral_f)),
                    "common_denominator_before_content": str(denominator), "content_removed": str(content),
                },
                "normalization": chart,
            }
        return self._chart_cache[key]

    def norm(self, word):
        key = word_key(word)
        if key not in self._norm_cache:
            v = vector(RealField(384), word)
            self._norm_cache[key] = v * self.H * v
        return self._norm_cache[key]

    def dot(self, left, right):
        return vector(RealField(384), left) * self.H * vector(RealField(384), right)

    def measurement(self, target_word, anchors):
        """Exact parameter/cancellation data plus numeric Neron--Tate residuals."""
        target = self.point(target_word)
        x, y = map(QQ, target.xy())
        target_norm = self.norm(target_word)
        rows = []
        for anchor_word in anchors:
            chart = self.chart(anchor_word)
            anchor = self.point(anchor_word)
            a, b = map(QQ, anchor.xy())
            denominator = x - a
            if not denominator:
                # This can only occur when a held-out target actually equals an
                # anchor.  The projective parameter is retained explicitly.
                rows.append({"anchor_word_in_D": list(anchor_word), "parameter": "infinity"})
                continue
            parameter = (y + b) / denominator
            numerator, denominator = map(int, (parameter.numerator(), parameter.denominator()))
            M = max(abs(numerator), denominator)
            # Homogenisations of N_Q and F_Q at (n,d).
            fixed = chart["primitive_integral_x_map"]
            F = evaluate_homogeneous(tuple(map(int, fixed["F_ascending"])), numerator, denominator)
            N = evaluate_homogeneous(tuple(map(int, fixed["N_ascending"])), numerator, denominator)
            from math import gcd
            cancellation = gcd(abs(N), abs(F))
            residual_word = add_words(tuple(2 * value for value in target_word), neg_word(anchor_word))
            residual_height = 4 * target_norm - 4 * self.dot(target_word, anchor_word) + self.norm(anchor_word)
            rows.append({
                "anchor_word_in_D": list(anchor_word),
                "parameter": {"numerator": str(numerator), "denominator": str(denominator), "height_log2": log2_exact(M)},
                "residual_word_in_D": list(residual_word),
                "residual_canonical_height": decimal_text(residual_height),
                "raw_x_map": {"N_homogeneous": str(N), "F_homogeneous": str(F)},
                "finite_cancellation_gcd": str(cancellation),
                "lambda_finite_bits": "%.70f" % (RealField(256)(cancellation).log(2) / 4),
                "real_distortion_bits": "%.70f" % (-RealField(256)(max(abs(N), abs(F))).log(2) / 4 + RealField(256)(M).log(2)),
                "normalization_maximum_coefficient_bits": chart["normalization"]["maximum_coefficient_bits"],
            })
        finite = [row for row in rows if row["parameter"] != "infinity"]
        if not finite:
            raise ArithmeticError("all chart parameters were projective")
        def best(field, numeric=float):
            return min(finite, key=lambda row: numeric(row[field]))
        return {
            "all_anchor_measurements": rows,
            "minima": {
                "kappa_parameter_height_bits": best("parameter", lambda value: float(value["height_log2"])),
                "r_residual_canonical_height": best("residual_canonical_height", float),
                "m_normalized_quartic_coefficient_bits": best("normalization_maximum_coefficient_bits", int),
                "lambda_finite_cancellation_bits": best("lambda_finite_bits", float),
            },
        }

    def quotient_geometry(self):
        combined = matrix(ZZ, [list(word) for word in (*self.generic_words, *self.exceptional_words)]).transpose()
        smith = [abs(int(value)) for value in combined.smith_form()[0].diagonal()]
        if combined.rank() != self.dimension:
            raise ArithmeticError("generic and exceptional words do not span D")
        generic = matrix(ZZ, [list(word) for word in self.generic_words]).transpose()
        generic_smith = [abs(int(value)) for value in generic.smith_form()[0].diagonal()]
        T = combined.change_ring(RealField(384))
        G = T.transpose() * self.H * T
        r = len(self.generic_words)
        A, B, C = G[:r, :r], G[:r, r:], G[r:, r:]
        quotient = C - B.transpose() * A.inverse() * B
        return {
            "D_rank": self.dimension,
            "M17_rank": len(self.generic_words),
            "known_exceptional_quotient_rank": len(self.exceptional_words),
            "M17_smith_in_D": generic_smith,
            "combined_basis_smith_in_D": smith,
            "M17_is_primitive_in_displayed_D": generic_smith == [1] * len(self.generic_words),
            "claim_boundary": "Primitivity is inside the displayed point lattice D, not a saturation assertion in E(Q).",
            "height_gram_in_D": [[decimal_text(value) for value in row] for row in self.H.rows()],
            "height_gram_in_M17_plus_exceptional_basis": [[decimal_text(value) for value in row] for row in G.rows()],
            "schur_complement_known_exceptional_quotient_form": [[decimal_text(value) for value in row] for row in quotient.rows()],
            "successive_minima_and_short_vectors": short_vector_summary(quotient),
            "exceptional_mod_2_images_in_quotient_basis": [list(unit_word(len(self.exceptional_words), i)) for i in range(len(self.exceptional_words))],
        }

    def build(self):
        if len(self.generic_words) != 17:
            raise ArithmeticError("comparison is defined only for an M17 source")
        geometry = self.quotient_geometry()
        atlas = {}
        panels = []
        # Native S_k uses the first k directions, and E_k,... are held out.
        # This makes every known exceptional point a target at its admission
        # stage, without leaking it into that stage's native bank.
        for k in range(len(self.exceptional_words)):
            subgroup = self.generic_words + self.exceptional_words[:k]
            native = standard_bank(subgroup)
            control = generic_control_bank(self.generic_words, len(native))
            native_ids, control_ids = [], []
            for prefix, bank, ids in (("native", native, native_ids), ("generic_control", control, control_ids)):
                for word in bank:
                    key = word_key(word)
                    if key not in atlas:
                        atlas[key] = self.chart(word)
                    ids.append(key)
            for target_index in range(k, len(self.exceptional_words)):
                target = self.exceptional_words[target_index]
                panels.append({
                    "stage": k,
                    "subgroup_rank": len(subgroup),
                    "target_exceptional_index_zero_based": target_index,
                    "target_word_in_D": list(target),
                    "native_bank": {"kind": "signed_units_and_pairs", "chart_count": len(native), "chart_ids": native_ids,
                                    "measurement": self.measurement(target, native)},
                    "generic_only_control": {"kind": "fixed_lexicographic_generic_shells", "chart_count": len(control), "chart_ids": control_ids,
                                             "measurement": self.measurement(target, control)},
                })
        return {
            "label": self.label,
            "curve": [qtext(value) for value in self.model],
            "displayed_D_points": [point_record(point) for point in self.displayed],
            "M17_specialization": [{"word_in_D": list(word), "point": point_record(self.point(word))} for word in self.generic_words],
            "exceptional_directions": [{"word_in_D": list(word), "point": point_record(self.point(word))} for word in self.exceptional_words],
            "rank_claim": self.rank_claim,
            "input_hashes": self.bindings,
            "canonical_height_precision_decimal_digits": 110,
            "canonical_height_maximum_asymmetry": str(self.height_asymmetry),
            "known_quotient": geometry,
            "chart_atlas": atlas,
            "nested_target_blind_panels": panels,
            "measurement_conventions": {
                "r": "minimum displayed Neron--Tate height of 2P-Q over the finite bank; numerical at the stated precision",
                "kappa": "minimum exact reduced parameter height log2(max(|n|,d))",
                "m": "minimum maximum coefficient bit length among the bank's fixed normalized quartics",
                "lambda": "finite cancellation (1/4) log2 gcd(N_Q(n,d),F_Q(n,d)); retained exactly through the gcd",
                "normalization_terms": "Each anchor records raw F_Q,N_Q and the factor-free integral/Gauss normalisation record.",
            },
        }


def curve302_inputs():
    parent = read(PARENT302)
    generic = matrix(ZZ, parent["basis_embedding_in_public_D"])
    generic_words = tuple(tuple(int(generic[row, column]) for row in range(31)) for column in range(17))
    basis = list(generic.change_ring(GF(2)).columns())
    exceptional = []
    for index in range(31):
        axis = vector(GF(2), unit_word(31, index))
        if matrix(GF(2), basis + [axis]).rank() > len(basis):
            basis.append(axis)
            exceptional.append(unit_word(31, index))
    expected = (0, 2, 3, 5, 8, 12, 13, 16, 18, 21, 25, 28, 29, 30)
    if tuple(word.index(1) for word in exceptional) != expected:
        raise ArithmeticError("Curve302 mod-2 complement changed")
    bindings = {str(path.relative_to(ROOT)): digest(path) for path in (PARENT302, CERT302, CAS / "icarm_curve302.py", CAS / "half_lattice_pointed_sieve.py", CAS / "prospective_half_lattice.sage", Path(__file__))}
    return Atlas(
        "Curve302 known M17 to displayed D31 (+14)", curve302.short_coefficients(), curve302.SHORT_POINTS,
        generic_words, exceptional, bindings,
        {"known_rank_lower_bound": 31, "interpretation": "D31 is a certified displayed independent lattice; no rank upper bound asserted."},
    )


def curve11952_inputs():
    seed, terminal, verified = read(SEED11952), read(TERM11952), read(VERIFY11952)
    if seed["points"] != terminal["points"][:17] or terminal["rank_lower_bound"] != 25:
        raise ArithmeticError("11952 M17 is not the frozen prefix of D25")
    if not str(verified.get("status", "")).startswith("PASS"):
        raise ArithmeticError("11952 terminal witness lacks its independent verification")
    rows = [row for signature in terminal["proof"]["signatures"] for row in signature["rows"]]
    if binary_rank(rows, 25) != 25 or binary_rank([[entry for entry in row[:17]] for row in rows], 17) != 17:
        raise ArithmeticError("11952 finite-reduction independence payload changed")
    generic = tuple(unit_word(25, index) for index in range(17))
    exceptional = tuple(unit_word(25, index) for index in range(17, 25))
    bindings = {str(path.relative_to(ROOT)): digest(path) for path in (SEED11952, TERM11952, VERIFY11952, CAS / "half_lattice_pointed_sieve.py", CAS / "prospective_half_lattice.sage", Path(__file__))}
    rank_claim = {
        "known_rank_lower_bound": 25,
        "interpretation": "D25 is a certified displayed independent lattice, so D25/M17 is the known rank-8 exceptional quotient. It is not the full exceptional quotient and does not give an upper bound on rank.",
        "finite_reduction_independence": {
            "mod_2_combined_rank_D25": 25,
            "mod_2_combined_rank_M17": 17,
            "no_rational_2_torsion_prime": terminal["proof"]["no_rational_2_torsion_prime"],
            "selected_primes": [signature["prime"] for signature in terminal["proof"]["signatures"]],
        },
    }
    return Atlas("11952 at t=921/653: known M17 to displayed D25 (+8)", seed["curve"], terminal["points"], generic, exceptional, bindings, rank_claim)


def compact_summary(payload):
    curves = []
    for curve in payload["curves"]:
        records = []
        for panel in curve["nested_target_blind_panels"]:
            def minima(bank):
                data = bank["measurement"]["minima"]
                return {key: value["anchor_word_in_D"] for key, value in data.items()}
            records.append({"stage": panel["stage"], "target_exceptional_index_zero_based": panel["target_exceptional_index_zero_based"],
                            "native_winners": minima(panel["native_bank"]), "generic_control_winners": minima(panel["generic_only_control"])})
        curves.append({"label": curve["label"], "rank_claim": curve["rank_claim"], "known_quotient": curve["known_quotient"], "winner_index": records})
    return {"schema": "elliptic-curves.rank-accessibility-atlas-summary.v1", "payload_sha256": sha256(canonical(payload)).hexdigest(), "curves": curves,
            "boundary": "The compressed payload retains all exact anchors and measurements. This summary contains only navigation data, not a rank inference."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="recompute and compare both deterministic payloads")
    args = parser.parse_args()
    inputs = [curve302_inputs(), curve11952_inputs()]
    payload = {
        "schema": "elliptic-curves.rank-accessibility-atlas.v1",
        "status": "COMPLETE_FINITE_ATLAS_MEASUREMENT",
        "claim_boundary": "No chart result is a new point search, rank admission, saturation statement in E(Q), or rank upper bound. 11952 remains rank at least 25.",
        "software": {"sage": SAGE_VERSION},
        "policy": {"anchor_banks": "native signed units/pairs; equal-cardinality generic-only lexicographic shells", "target_leakage": "Each target is held out of its native S_k bank.", "chart_coordinate": "t_Q(P)=(y(P)+y(Q))/(x(P)-x(Q)) on the short model."},
        "curves": [item.build() for item in inputs],
    }
    data = canonical(payload)
    compressed = gzip.compress(data, mtime=0)
    summary = canonical(compact_summary(payload))
    if args.check:
        if not OUT.is_file() or OUT.read_bytes() != compressed or not SUMMARY.is_file() or SUMMARY.read_bytes() != summary:
            raise SystemExit("RANK_ACCESSIBILITY_ATLAS_CHECK|status=FAIL")
        print("RANK_ACCESSIBILITY_ATLAS_CHECK|status=PASS")
        return
    ART.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(compressed)
    SUMMARY.write_bytes(summary)
    print("RANK_ACCESSIBILITY_ATLAS|status=PASS|payload_sha256=%s|output=%s" % (sha256(data).hexdigest(), OUT))


if __name__ == "__main__":
    main()
