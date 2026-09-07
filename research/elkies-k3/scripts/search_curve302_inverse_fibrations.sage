#!/usr/bin/env sage-python
"""Checkpointed, rank-unrestricted norm-eight extension of the 302 inverse search.

Search budgets select equally from every minimum-vector multiplicity stratum;
they never select a desired MW rank.  Only exact specialization witnesses can
enter the overlap ranking.  Modular survival is not an overlap measurement.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
SCRIPTS = ROOT / "elkies-k3/scripts"
SOURCES = {
    "103b2": "103b2-norm8",
    "08f72": "08f72-alternate-norm8",
}
PRIMES = (1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051,
          1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109,
          1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187,
          1193, 1201, 1213, 1217, 1223, 1229, 1231)


def load(name, path):
    return SourceFileLoader(name, str(path)).load_module()


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def write(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def projective_no_root(polynomial, degree, prime):
    """Use FIXED homogeneous degree, retaining infinity after degree drops."""
    if polynomial.degree() > degree:
        raise ArithmeticError("comparison exceeds declared projective degree")
    if not polynomial or polynomial[degree] == 0:
        return False
    t = polynomial.parent().gen()
    return polynomial.gcd(t**prime - t).degree() == 0


def modular_context(model, prime, helpers):
    field = GF(prime)
    ring = PolynomialRing(field, "u")
    rational = lambda x: helpers.reduce_rational(x, field)
    a = ring([rational(x) for x in model["weierstrass_model"]["A_coefficients_low_to_high"]])
    b = ring([rational(x) for x in model["weierstrass_model"]["B_coefficients_low_to_high"]])
    if a.degree() != 8 or b.degree() != 12:
        raise ArithmeticError("source degrees drop")
    curve = EllipticCurve(ring.fraction_field(), [a, b])
    basis = [curve(helpers.polynomial_from_record(r["X"], ring, field),
                   helpers.polynomial_from_record(r["Y"], ring, field))
             for r in model["sections"]["records"]]
    return field, ring, a, curve, basis


def comparison_for_word(context, word, c4, c6, helpers, chord):
    field, ring, a, curve, basis = context
    trace = sum((int(c)*p for c, p in zip(word, basis) if c), curve(0))
    if trace.is_zero():
        raise ArithmeticError("trace is zero after reduction")
    frame = chord.trace_chord_frame(trace[0], trace[1], ring)
    h, nx, ny, m0 = (frame[k] for k in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("finite pole degree is not two")
    lr = PolynomialRing(field, "t")
    br = PolynomialRing(lr, "u")
    h, nx, ny, m = br(h), br(nx), br(ny), br(m0) + lr.gen()*br(h)**2
    numerator = m**4 - 6*m**2*nx - 8*m*ny - 3*nx**2 - 4*br(a)*h**4
    quartic, remainder = numerator.quo_rem(h**6)
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("residual chord is not a quartic")
    i, j = helpers.binary_quartic_invariants(quartic, lr)
    if i.degree() > 8 or j.degree() > 12:
        raise ArithmeticError("Jacobian exceeds K3 degree bounds")
    ta = -27*helpers.reduce_rational(c4, field)
    tb = -54*helpers.reduce_rational(c6, field)
    return (-27*i)**3 * tb**2 - ta**3 * (-27*j)**2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-stratum", type=int, default=128)
    parser.add_argument("--seconds", type=int, default=1800)
    parser.add_argument("--output", type=Path, default=GEN / "elkies-k3-curve302-inverse-fibration-extension-v1.json")
    parser.add_argument("--checkpoint-dir", type=Path, default=ROOT / "artifacts/local/elkies-k3/curve302-inverse-extension-v1")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.per_stratum < 1 or args.seconds < 1:
        parser.error("budgets must be positive")
    helper_path = SCRIPTS / "screen_icarm_curve398_norm8_a1_fibrations.sage"
    chord_path = SCRIPTS / "construct_elkies_2026_bisections.sage"
    public_path = ROOT / "elliptic-curves/cas/icarm_curve302.py"
    helpers = load("inverse302_helpers", helper_path)
    chord = load("inverse302_chord", chord_path)
    public = load("inverse302_public", public_path)
    target = EllipticCurve(QQ, [QQ(str(a)) for a in public.GENERAL_WEIERSTRASS_COEFFICIENTS])
    points = [target(QQ(str(x)), QQ(str(y))) for x, y in public.POINTS]
    if len(points) != 31:
        raise ArithmeticError("expected the complete public 31-point configuration")
    c4, c6 = target.c_invariants()
    started = time.monotonic()
    sources = []
    inputs = {relative(p): digest(p) for p in (Path(__file__), helper_path, chord_path, public_path)}
    for label, prefix in SOURCES.items():
        mp = GEN / f"elkies-k3-r17-norm12-orbit{label}-direct-fibration-v1.json"
        tp = GEN / f"elkies-k3-r17-norm12-{prefix}-pencil-priority-v1.tsv"
        cp = tp.with_suffix(".json")
        model, census = json.loads(mp.read_text()), json.loads(cp.read_text())
        if model["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
            raise ArithmeticError("source section basis is not certified saturated")
        if census["priority_table_sha256"] != digest(tp):
            raise ArithmeticError("norm-eight census table hash mismatch")
        gram = matrix(ZZ, model["sections"]["height_gram"])
        rows = list(csv.DictReader(tp.open(), delimiter="\t"))
        if len(rows) != census["class_count"]:
            raise ArithmeticError("census class count mismatch")
        counts = Counter()
        selected = []
        for row in rows:
            m = int(row["minimal_unoriented_count"])
            if counts[m] >= args.per_stratum:
                continue
            word = tuple(map(int, row["section_basis_w"].split()))
            v = vector(ZZ, word)
            if len(word) != 17 or v*gram*v != 8:
                raise ArithmeticError("trace has wrong source height")
            selected.append({"priority_rank": int(row["priority_rank"]),
                             "minimum_unoriented_count": m,
                             "trace_word": list(word), "status": "UNKNOWN_UNTESTED",
                             "certified_overlap_rank": None})
            counts[m] += 1
        source_inputs = {relative(p): digest(p) for p in (mp, tp, cp)}
        inputs.update(source_inputs)
        identity = {"inputs": {**inputs}, "source": label, "per_stratum": args.per_stratum,
                    "primes": list(PRIMES), "selected_priority_ranks": [r["priority_rank"] for r in selected]}
        identity_hash = sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        checkpoint = args.checkpoint_dir / f"{label}.json"
        completed = []
        if args.resume and checkpoint.exists():
            saved = json.loads(checkpoint.read_text())
            if saved["identity_hash"] != identity_hash:
                raise ArithmeticError("resume inputs/protocol changed")
            selected, completed = saved["records"], saved["completed_primes"]
        elif checkpoint.exists():
            raise FileExistsError(f"use --resume or a fresh checkpoint directory: {checkpoint}")
        for prime in PRIMES:
            if prime in completed:
                continue
            pending = [r for r in selected if r["status"] != "EXCLUDED_NO_RATIONAL_PARAMETER"]
            if not pending or time.monotonic() - started >= args.seconds:
                break
            try:
                if target.discriminant().numerator() % prime == 0:
                    raise ArithmeticError("target bad reduction")
                context = modular_context(model, prime, helpers)
            except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                print(f"INVERSE302|source={label}|prime={prime}|unavailable={error}", flush=True)
                continue
            finished_pass = True
            for position, row in enumerate(pending):
                if time.monotonic() - started >= args.seconds:
                    finished_pass = False
                    break
                try:
                    comparison = comparison_for_word(context, row["trace_word"], c4, c6, helpers, chord)
                    if projective_no_root(comparison, 24, prime):
                        row.update(status="EXCLUDED_NO_RATIONAL_PARAMETER", witness_prime=prime,
                                   homogeneous_comparison_coefficients_mod_p=[int(comparison[i]) for i in range(25)])
                    else:
                        row["status"] = "UNKNOWN_MODULAR_SURVIVOR"
                except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                    row["status"] = "UNKNOWN_UNAVAILABLE_REDUCTION"
                    row["last_unavailable_reason"] = str(error)
                if (position + 1) % 2048 == 0:
                    write(checkpoint, {"identity_hash": identity_hash, "identity": identity,
                                       "completed_primes": completed, "records": selected})
                    print(f"INVERSE302|source={label}|prime={prime}|processed={position+1}/{len(pending)}", flush=True)
            if finished_pass:
                completed.append(prime)
            write(checkpoint, {"identity_hash": identity_hash, "identity": identity,
                               "completed_primes": completed, "records": selected})
            remaining = sum(r["status"] != "EXCLUDED_NO_RATIONAL_PARAMETER" for r in selected)
            print(f"INVERSE302|source={label}|prime={prime}|remaining={remaining}/{len(selected)}|seconds={time.monotonic()-started:.1f}", flush=True)
        sources.append({"source_chart": label, "available_classes": len(rows),
                        "selected_by_minimum_multiplicity": dict(sorted(counts.items())),
                        "completed_primes": completed, "records": selected})
    records = [r for s in sources for r in s["records"]]
    remaining = sum(r["status"] != "EXCLUDED_NO_RATIONAL_PARAMETER" for r in records)
    result = {
        "schema": "elkies-k3.curve302-inverse-fibration-extension.v1",
        "status": "BOUNDED_EXACT_NO_SPECIALIZATION" if not remaining else "UNKNOWN_SURVIVORS_REQUIRE_EXACT_RECONSTRUCTION",
        "target": {"curve_id": 302, "public_points_checked_on_curve": len(points),
                   "known_subgroup_rank": 31, "whole_curve_exact_rank": "UNKNOWN",
                   "ainvs": list(map(str, target.a_invariants()))},
        "protocol": {"per_multiplicity_stratum": args.per_stratum, "prime_chain": list(PRIMES),
                     "wall_seconds_limit_per_invocation": args.seconds,
                     "generic_rank_filter": None, "source_picard_rank": 19,
                     "fixed_source_generic_rank_ceiling": 17,
                     "all_fibrations_exhausted": False,
                     "selection": "First equation-cost-ordered representatives equally across every minimum multiplicity; no target-rank filter."},
        "selected_count": len(records), "excluded_count": len(records)-remaining,
        "unresolved_count": remaining, "sources": sources, "overlap_ranking": [],
        "overlap_definition": "dim_Q(span_Q(sp_t(M_candidate)) intersect span_Q(P1,...,P31)), after an exact Q-isomorphism to E302",
        "overlap_gate": "Require exact rational family sections, parameter, Q-isomorphism, specialized group-law coordinate identities in the 31-point span, and exact rational matrix rank. Without witnesses overlap stays null, never zero.",
        "proof_boundary": "A no-root homogeneous degree-24 comparison over F_p excludes rational specialization of that declared chord pencil. Minimum multiplicity is a sampling stratum, not a newly certified generic rank. This searches two fixed-NS charts, not inverse NS reconstruction from the specialized height lattice. No parent family, absence of all parents, or rank upper bound is inferred.",
        "inputs": inputs, "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": "sage -python elkies-k3/scripts/search_curve302_inverse_fibrations.sage --per-stratum " + str(args.per_stratum) + " --seconds " + str(args.seconds) + " --resume",
    }
    write(args.output, result)
    print(f"INVERSE302|selected={len(records)}|excluded={len(records)-remaining}|unresolved={remaining}|output={relative(args.output)}", flush=True)


if __name__ == "__main__":
    main()
