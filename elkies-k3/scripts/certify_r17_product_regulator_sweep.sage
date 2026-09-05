#!/usr/bin/env sage-python
"""Exact regulator obstruction for the four remaining product twists.

Reuse the pinned 131/137 polynomials. Only 0f82c:025be needs p=151,
because its p=131 analytic rank is two. No section or Selmer search runs.
Rank-one local/BSD normalization is shared with the original certificate;
all input hashes and independent Magma controls are checked locally.
"""

from __future__ import annotations

import argparse
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from sage.all import PolynomialRing, QQ
from sage.version import version


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
SHARED = Path(__file__).resolve().with_name("certify_r17_product_19bad_083ad_rank_zero.sage")
CLASSIFICATION = RESULTS / "elkies-k3-r17-all17-product-toric-frobenius-campaign-v1.json"
OUTPUT = RESULTS / "elkies-k3-r17-product-regulator-sweep-v1.json"
CONTROLS = RESULTS / "r17-product-regulator-sweep-controls"
TARGETS = ("11ee2:0c36e", "0c10b:17a1a", "0f82c:025be", "11ae6:0f82c")
EXTRA_PRIMES = {"0f82c:025be": (151,)}


def load_shared():
    spec = importlib.util.spec_from_loader(
        "r17_product_regulator_shared", SourceFileLoader("r17_product_regulator_shared", str(SHARED))
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def analytic_rank(frobenius):
    """Use exact polynomial division, including ranks other than one."""
    ring = PolynomialRing(QQ, "T")
    T = ring.gen()
    p = int(frobenius["prime"])
    F = ring(frobenius["elliptic_L"]["frobenius_characteristic_coefficients_low_to_high"])
    if F.degree() != 28 or not F.is_monic():
        raise ArithmeticError("not a complete degree-28 characteristic polynomial")
    L = F.reverse()
    moments = [str(-L[1]), str(L[1]**2 - 2*L[2])]
    if L[0] != 1 or moments != frobenius["elliptic_L"]["power_sums_n1_n2"]:
        raise ArithmeticError("L-polynomial orientation or moment mismatch")
    quotient, rank = L, 0
    while quotient(QQ(1)/p) == 0:
        quotient, remainder = quotient.quo_rem(1-p*T)
        if remainder:
            raise ArithmeticError("nonexact rank division")
        rank += 1
    if quotient(QQ(1)/p) <= 0:
        raise ArithmeticError("nonpositive normalized leading coefficient")
    stored = sum(int(hit["multiplicity"]) for hit in
                 frobenius["elliptic_L"]["cyclotomic_hits_after_T_equals_pZ"]
                 if hit["order"] == 1)
    if rank != stored:
        raise ArithmeticError("fixed-factor multiplicity mismatch")
    return rank, L


def compare_rank_one(reductions, squareclass):
    eligible = [row for row in reductions if row["analytic_rank"] == 1]
    if len(eligible) < 2:
        raise ArithmeticError("two rank-one good reductions are required for this gate")
    comparisons = []
    first = eligible[0]
    for second in eligible[1:]:
        ratio = (QQ(first["regulator_times_sha_over_torsion_squared_if_rank_one"]) /
                 QQ(second["regulator_times_sha_over_torsion_squared_if_rank_one"]))
        is_square = bool(ratio.is_square())
        witness = next(({"prime": int(p), "valuation": int(e)}
                        for p, e in ratio.factor() if e % 2), None)
        if is_square != (witness is None):
            raise ArithmeticError("square test and valuation witness disagree")
        comparisons.append({
            "primes": [first["prime"], second["prime"]],
            "regulator_ratio": str(ratio), "ratio_squareclass": squareclass(ratio),
            "is_rational_square": is_square, "odd_valuation_witness": witness,
        })
    return comparisons


def self_test():
    """Cheap regressions for the two consequential rank/squareclass branches."""
    shared = load_shared()

    def row(prime, value):
        return {"prime": prime, "analytic_rank": 1,
                "regulator_times_sha_over_torsion_squared_if_rank_one": str(value)}

    # Unknown specialization indices multiply heights by rational squares.
    compatible = compare_rank_one([row(131, 18), row(137, 50)], shared.squareclass)
    assert compatible[0]["regulator_ratio"] == "9/25"
    assert compatible[0]["is_rational_square"]
    assert compatible[0]["odd_valuation_witness"] is None
    incompatible = compare_rank_one([row(131, 18), row(137, 75)], shared.squareclass)
    assert not incompatible[0]["is_rational_square"]
    assert incompatible[0]["ratio_squareclass"] == 6
    witness = incompatible[0]["odd_valuation_witness"]
    assert QQ("6/25").valuation(witness["prime"]) == witness["valuation"]
    # Rank two must neither supply a height squareclass nor count as a second
    # usable reduction, even if a bogus squareclass was attached to its row.
    try:
        compare_rank_one([row(137, 2), {"prime": 131, "analytic_rank": 2,
                                     "regulator_squareclass_if_rank_one": 2}], shared.squareclass)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("rank-two reduction accepted by the rank-one gate")
    path = RESULTS / "elkies-k3-r17-product-alternate-orbit-0f82c--alternate-orbit-025be-p131-toric-frobenius-v1.json"
    frobenius = json.loads(path.read_text())
    assert analytic_rank(frobenius)[0] == 2
    frobenius["elliptic_L"]["power_sums_n1_n2"][0] = "0"
    try:
        analytic_rank(frobenius)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("incorrect L-polynomial moment accepted")
    print("PRODUCTREGULATOR|square_rescaling=PASS|nonsquare_witness=PASS|rank_two_rejection=PASS|moment_tamper_rejection=PASS")


def build_payload(export_magma=False):
    shared = load_shared()
    inputs = {}

    def pin(path):
        value = shared.digest(path)
        inputs[shared.relative(path)] = value
        return value

    def read_pinned(path):
        pin(path)
        return json.loads(path.read_text())

    pin(Path(__file__).resolve())
    pin(SHARED)
    direct = read_pinned(shared.DIRECT)
    pairs = read_pinned(shared.PAIRS)
    height = read_pinned(shared.HEIGHT)
    campaign = read_pinned(CLASSIFICATION)
    if direct["status"] != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ArithmeticError("uncertified direct model")
    if height["status"] != "PASS_EXACT_PRODUCT_SURVIVOR_GALOIS_HEIGHT_GATE":
        raise ArithmeticError("uncertified characteristic-zero component/height gate")
    height_by_key = {row["pair_key"]: row for row in height["records"]}
    campaign_by_key = {row["pair_key"]: row for row in campaign["targets"]}
    records, controls = [], []
    for tag in TARGETS:
        key = ":".join("alternate-orbit-"+part for part in tag.split(":"))
        pair = next(row for row in pairs["pairs"] if row["pair_key"] == key)
        gate = height_by_key[key]
        if gate["arithmetic_rank_gate"]["QQ_u_mw_rank_upper_bound"] != 1:
            raise ArithmeticError("arithmetic rank-at-most-one gate missing")
        if any(row["two_division_factor_degrees_over_residue_field"] != [3]
               for row in gate["branch_places"]):
            raise ArithmeticError("rational component gate missing")
        reductions = []
        for prime in (131, 137) + EXTRA_PRIMES.get(tag, ()):
            path = RESULTS / f"elkies-k3-r17-product-{key.replace(':', '--')}-p{prime}-toric-frobenius-v1.json"
            frobenius = read_pinned(path)
            if (frobenius["pair_key"] != key or frobenius["prime"] != prime or
                frobenius["status"] != "PASS_COMPLETE_FROBENIUS_PICARD_BOUND" or
                frobenius["good_reduction"]["status"] != "PASS"):
                raise ArithmeticError("uncertified or misidentified Frobenius input")
            # Bind the complete polynomial to the direct model and point-count audit.
            for name, expected in frobenius["inputs"].items():
                if pin(ROOT / name) != expected:
                    raise ArithmeticError(f"stale Frobenius input: {name}")
            audit_paths = [ROOT / name for name in frobenius["inputs"] if "audit" in name]
            if len(audit_paths) != 1:
                raise ArithmeticError("unique independent point-count audit required")
            audit = read_pinned(audit_paths[0])
            audit_target = next(row for row in audit["targets"] if row["pair_key"] == key)
            if ([QQ(c) for c in audit_target["product_quartic_coefficients_low_to_high"]] !=
                [QQ(c) for c in pair["product_quartic_coefficients_low_to_high"]]):
                raise ArithmeticError("twist differs from the audited character, including its constant")
            audit_reduction = next(row for row in audit_target["reductions"] if row["prime"] == prime)
            if (list(map(QQ, audit_reduction["elliptic_L_frobenius_power_sums_n1_n2"])) !=
                list(map(QQ, frobenius["elliptic_L"]["power_sums_n1_n2"]))):
                raise ArithmeticError("independent moment audit disagrees")
            if prime in (131, 137):
                anchor = next(row for row in campaign_by_key[key]["reductions"] if row["prime"] == prime)
                if anchor["certificate_sha256"] != pin(path):
                    raise ArithmeticError("Frobenius certificate differs from campaign anchor")
            else:
                raw_dir = CONTROLS / f"toric-{tag.replace(':', '-')}-p{prime}"
                for extension in ("input", "output"):
                    raw_path = raw_dir / f"toric-controlled-reduction.{extension}"
                    if pin(raw_path) != frobenius["software"][f"ToricControlledReduction_{extension}_sha256"]:
                        raise ArithmeticError("retained toric replay bytes differ from certificate")
                pin(raw_dir / "input-certificate.json")
                pin(raw_dir / "toric-controlled-reduction.log")
            rank, L = analytic_rank(frobenius)
            if rank != 1:
                # A hypothetical section does not force rank equality here, and
                # a rank-two regulator would not determine an individual height.
                reductions.append({
                    "prime": prime, "analytic_rank": rank,
                    "L_coefficients_low_to_high": shared.coefficients(L),
                    "good_reduction": True,
                    "rank_one_squareclass_test": "NOT_APPLICABLE",
                    "regulator_squareclass_if_rank_one": None,
                    "reason": "This reduction cannot force algebraic rank = analytic rank = 1.",
                    "frobenius_certificate": shared.relative(path),
                })
                continue
            row, magma = shared.reduction_record(prime, direct["weierstrass_model"], pair, frobenius, key)
            actual = sorted((r["base_place_degree"], sorted(r["residual_two_division_factor_degrees"]))
                            for r in row["branch_places"])
            expected = sorted((r["branch_degree"], sorted(r["cubic_factor_degrees_over_residue_field"]))
                              for r in frobenius["boundary"]["branch_factor_records"])
            if actual != expected:
                raise ArithmeticError("local factorization differs from certified boundary factor")
            row["frobenius_certificate"] = shared.relative(path)
            row["rank_one_squareclass_test"] = "APPLICABLE"
            reductions.append(row)
            job_path = CONTROLS / f"{tag.replace(':', '--')}-p{prime}.m"
            xml_path = job_path.with_suffix(".xml")
            controls.append((job_path, xml_path, magma, row, tag))

        comparisons = compare_rank_one(reductions, shared.squareclass)
        excluded = any(not row["is_rational_square"] for row in comparisons)
        records.append({
            "pair_key": key, "short_label": tag, "reductions": reductions,
            "comparisons": comparisons,
            "classification": "rank 0" if excluded else "regulator-compatible with rank 1",
            "status": "PROVED_ARITHMETIC_PRODUCT_TWIST_RANK_ZERO" if excluded else "REGULATOR_COMPATIBLE_RANK_UNKNOWN",
            "rank_over_QQ_u": {"lower": 0, "upper": 0 if excluded else 1,
                                 "status": "PROVED" if excluded else "UNKNOWN"},
            "rank_over_QQbar_u": {"lower": 0, "upper": 2, "status": "UNKNOWN"},
            "section_solving_eligible": not excluded,
            "rational_section_boxes": "ALL_POSITIVE_HEIGHTS_PROVED_EMPTY" if excluded else "UNKNOWN",
        })

    for job_path, xml_path, magma, row, tag in controls:
        if export_magma:
            CONTROLS.mkdir(parents=True, exist_ok=True)
            job_path.write_text(magma)
            continue
        if job_path.read_text() != magma:
            raise ArithmeticError("stale independent-control input")
        raw = ET.fromstring(xml_path.read_text())
        if raw.findtext("headers/warning"):
            raise ArithmeticError("independent control has a warning or error")
        lines = [line.text or "" for line in raw.findall("results/line")]
        expected = f"ANALYTIC_INFORMATION <1, 2, {row['regulator_times_sha_over_torsion_squared_if_rank_one']}>"
        if (expected not in lines or "COMPLETE" not in lines or
            f"TARGET|{tag}|p={row['prime']}" not in lines):
            raise ArithmeticError("independent local/BSD normalization does not agree")
        row["independent_magma_control"] = {
            "version": raw.findtext("headers/version"), "result": expected,
            "boundary": "Uses the supplied certified L-polynomial; independently checks local/BSD normalization only.",
        }
        pin(job_path)
        pin(xml_path)

    compatible = [row["pair_key"] for row in records if row["section_solving_eligible"]]
    return {
        "schema": "elkies-k3.r17-product-regulator-sweep.v1",
        "status": "PASS_EXACT_REGULATOR_OBSTRUCTION_SWEEP",
        "targets": records,
        "rank_zero_count": len(records)-len(compatible),
        "regulator_compatible_count": len(compatible),
        "explicit_section_solving_queue": compatible,
        "new_frobenius_computations": [{
            "short_label": "0f82c:025be", "prime": 151,
            "reason": "The existing p=131 polynomial has analytic rank two; only p=137 supplies a rank-one regulator constraint.",
            "replay": "bash elkies-k3/scripts/run_r17_product_toric_frobenius_extra_prime.sh 'alternate-orbit-0f82c:alternate-orbit-025be' 151",
        }],
        "logic": [
            "The torsion-free 48I1 double cover makes every nonzero QQ(u)-section nontorsion.",
            "Good surface reduction preserving the geometric trivial lattice preserves its positive Shioda height.",
            "A nonzero reduction at an analytic-rank-one prime forces rank equality, hence refined BSD and finite square-order Sha.",
            "The section height has squareclass p^(chi-1)*L_star/product(c_v), with chi=4 and Tamagawa factors taken once per closed place.",
            "Two unequal squareclasses exclude every positive rational section height.",
            "Matching squareclasses are only a necessary rank-one condition and keep arithmetic rank UNKNOWN in [0,1].",
        ],
        "proof_boundary": "Arithmetic QQ(u) ranks only. Geometric ranks remain UNKNOWN in [0,2]; no full Selmer computation or section search is claimed.",
        "theoretical_references": [
            {"url": "https://math.stanford.edu/~conrad/BSDseminar/refs/Ulmer.pdf",
             "location": "Theorems 6.2.6, 6.3.1; Proposition 6.3.3"},
            {"path": "elkies-k3/R17_PRODUCT_19BAD_083AD_ARITHMETIC_RANK_ZERO_2026-09-05.md",
             "location": "Why this is an unconditional obstruction (height preservation and BSD normalization)"},
        ],
        "software": version, "inputs": inputs,
        "reproducing_command": "sage -python elkies-k3/scripts/certify_r17_product_regulator_sweep.sage --check",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--export-magma", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.self_test:
        if args.check or args.export_magma:
            parser.error("run --self-test separately from certificate replay/export")
        self_test()
        return
    if args.check and args.export_magma:
        parser.error("--check is read-only; export Magma separately")
    payload = build_payload(args.export_magma)
    if args.export_magma:
        print(f"PRODUCTREGULATOR|magma_inputs=EXPORTED|directory={CONTROLS}")
        return
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale or missing regulator-sweep certificate")
    else:
        args.output.write_text(rendered)
    for row in payload["targets"]:
        print(f"PRODUCTREGULATOR|pair={row['short_label']}|classification={row['classification']}")
    print(f"PRODUCTREGULATOR|rank_zero={payload['rank_zero_count']}|compatible={payload['regulator_compatible_count']}|status={'PASS' if args.check else 'WROTE'}")


if __name__ == "__main__":
    main()
