#!/usr/bin/env sage-python
"""Prepare all 104 MW16 finalists; benchmark only nine frozen size-selected charts.

No new parameters, adaptive waves, descent, or exceptional points. All 856
existing deepest masks are measured. Candidate checkpoints precede a separate
tiny benchmark, so point-search outcomes cannot affect coordinate selection.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import EllipticCurve, QQ, ZZ, version

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))
import mw16_model_size as ms
from mod2_reduction_independence import (
    combined_mod2_rank, find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)

ART = ROOT / "artifacts/generated-results/elliptic-curves"
INPUT = ART / "a1_mw16_target_free_parameter_candidates_h300_v1.json"
MASKS = ART / "a1_mw16_target_free_parameter_search_h300_v1.json"
OUTPUT = ART / "mw16_short_models_h300_v1.json.gz"
SUMMARY = ART / "mw16_short_models_h300_summary_v1.json"
BENCHMARK = ART / "mw16_short_models_chart_benchmark_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/mw16-short-models-square-content"
LEGACY = CAS / "run_curve385_iterated_half_lattice_search.sage"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, indent=2, sort_keys=True)+"\n").encode()
    if path.suffix == ".gz":
        data = gzip.compress(data, mtime=0)
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_bytes(data)
    temp.replace(path)


def read(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == ".gz" else data)


def prepare(candidate, masks, legacy):
    E = EllipticCurve(QQ, candidate["raw_short_model"])
    points = [ms.read_point(E, p) for p in candidate["raw_generic_points"]]
    if len(points) != 16 or not all(points):
        raise ArithmeticError("expected sixteen finite sections")
    # PARI's rational minimal-model algorithm is cheap on these fibres; do
    # not use the number-field global_minimal_model path (factors Delta).
    from cysignals.alarm import alarm, cancel_alarm
    alarm(2)
    try:
        minimal = E.minimal_model()
    finally:
        cancel_alarm()
    phi_min, min_points = ms.transport(E, minimal, points)
    short = EllipticCurve(QQ, [-minimal.c4()/48, -minimal.c6()/864])
    integral_short = minimal.short_weierstrass_model()
    square = EllipticCurve(QQ, [0, minimal.b2()/4, 0, minimal.b4()/2, minimal.b6()/4])
    phi_sq, sq_points = ms.transport(E, square, points)
    r = min(sq_points, key=lambda p: (ms.bits(p[0]), p[0]))[0]
    centered = EllipticCurve(QQ, [0, square.a2()+3*r, 0,
        square.a4()+2*square.a2()*r+3*r*r,
        square.a6()+square.a4()*r+square.a2()*r*r+r**3])
    variants = {"raw": E, "minimal": minimal, "rational_short": short,
                "integral_short": integral_short, "completed_square": square,
                "section_centered": centered}
    metrics = []
    for name, curve in variants.items():
        phi, transported = ms.transport(E, curve, points)
        metrics.append({"name": name, "model": list(map(ms.qstr, curve.ainvs())),
                        "map_from_raw": ms.map_record(phi),
                        "model_size": ms.size(curve.ainvs()),
                        "section_size": ms.size([v for p in transported for v in p.xy()]),
                        "section_identities_checked": 16, "round_trips_checked": 16})
    # A short model is required by the existing finite-reduction checker.
    name = min(("raw", "rational_short", "integral_short"),
               key=lambda n: (*ms.size(variants[n].ainvs()).values(), n))
    S = variants[name]
    phi, pts = ms.transport(E, S, points)
    model = tuple(Fraction(str(v)) for v in S.ainvs())
    signatures = find_mod2_reduction_certificate(
        model, tuple(tuple(Fraction(str(v)) for v in p.xy()) for p in pts), prime_bound=1000)
    if combined_mod2_rank(signatures, 16) != 16:
        raise ArithmeticError("finite-reduction independence did not certify rank sixteen")
    independence = legacy.signature_record(signatures, 16)
    independence["no_rational_two_torsion_prime"] = find_two_torsion_certificate_prime(model)
    gram = tuple(tuple(int(2*QQ(v)) for v in row) for row in candidate["generic_height_gram"])
    oracle = legacy.CosetOracle(gram)
    charts = []
    for mask in masks:
        norm, representative, error = oracle.solve(tuple((mask >> i) & 1 for i in range(16)))
        if norm != 23 or error > 1.e-6:
            raise ArithmeticError("frozen deepest-mask norm failed")
        Q = sum((v*p for v, p in zip(representative, pts)), S(0))
        raw_Q = sum((v*p for v, p in zip(representative, points)), E(0))
        if phi(raw_Q) != Q or not Q:
            raise ArithmeticError("chart base-point group-law transport failed")
        raw_f, f = ms.pointed(E, raw_Q), ms.pointed(S, Q)
        direct_g, direct_scale = ms.normalize(f)
        selected = ms.select_chart(S, Q, pts)
        # Pure x translations cancel from the pointed polynomial exactly.
        q_sq = phi_sq(raw_Q)
        center_phi = square.isomorphism_to(centered)
        # An isomorphism may choose the negative y sign; account for it.
        center_f = ms.pointed(centered, center_phi(q_sq))
        if center_f not in (ms.pointed(square, q_sq), ms.pointed(square, q_sq)(-ms.z)):
            raise ArithmeticError("x-translation invariance failed")
        charts.append({"mask": mask, "representative": list(representative),
            "twice_generic_norm": norm, "raw_base_point": ms.point_record(raw_Q),
            "base_point": ms.point_record(Q), "raw_quartic_size": ms.size([raw_f[i] for i in range(5)]),
            "short_quartic_size": ms.size([f[i] for i in range(5)]),
            "short_integral_quartic_size": ms.size([direct_g[i] for i in range(5)]),
            "selected": selected})
    j = E.j_invariant()
    return {"candidate_id": candidate["candidate_id"], "presentation_id": candidate["presentation_id"],
        "parameter": candidate["parameter"], "status": "PASS_EXACT_MODEL_SECTION_TRANSPORT",
        "global_minimal_model": list(map(ms.qstr, minimal.ainvs())),
        "raw_to_minimal": ms.map_record(phi_min),
        "minimal_sections": list(map(ms.point_record, min_points)),
        "selected_short_model_name": name, "selected_short_model": list(map(ms.qstr, S.ainvs())),
        "raw_to_selected_short": ms.map_record(phi),
        "selected_short_sections": list(map(ms.point_record, pts)),
        "model_trials": metrics, "independence": independence,
        "j_numerator_bits": int(abs(j.numerator()).nbits()),
        "j_denominator_bits": int(j.denominator().nbits()),
        "any_integral_binary_quartic_coefficient_bit_lower_bound": ms.integral_quartic_bit_lower_bound(j),
        "charts": charts}


def search_chart(candidate, row, chart, mode, height, timeout):
    E = EllipticCurve(QQ, candidate["raw_short_model"])
    S = EllipticCurve(QQ, row["selected_short_model"])
    phi = E.isomorphism_to(S)
    Q = ms.read_point(S, chart["base_point"])
    if mode == "raw_direct":
        curve, q, map_to = E, ms.read_point(E, chart["raw_base_point"]), E.isomorphism_to(E)
    else:
        curve, q, map_to = S, Q, phi
    if mode in ("selected_direct", "selected_reduced"):
        record = chart["selected"]
    else:
        g, scale = ms.normalize(ms.pointed(curve, q))
        record = {"matrix_a_b_c_d": ["1", "0", "0", "1"],
            "ordinate_scale": ms.qstr(scale),
            "integral_coefficients_ascending": [ms.qstr(g[i]) for i in range(5)]}
    g = ms.R(record["integral_coefficients_ascending"])
    poly = str(g).replace("z", "x")
    reduced = mode == "selected_reduced"
    setup = 'C=hyperellred([g,0],&m);' if reduced else 'C=[g,0];'
    inverse = ('d=m[2][2,1]*p[1]+m[2][2,2];if(d==0,print("POLE"),'
               'p=[(m[2][1,1]*p[1]+m[2][1,2])/d,(m[1]*p[2]+subst(m[3],x,p[1]))/d^2];'
               'if(p[2]^2!=subst(g,x,p[1]),error("reduction inverse failed"));print("POINT|",p[1],"|",p[2]))') if reduced else 'print("POINT|",p[1],"|",p[2])'
    program = (f'g={poly};print("PHASE|reduction");gettime();{setup}'
        'print("REDMS|",gettime());print("PHASE|search");gettime();'
        f'V=hyperellratpoints(C,{height});print("SEARCHMS|",gettime());print("COUNT|",#V);'
        f'for(i=1,#V,p=V[i];{inverse});print("DONE");quit\n')
    started = time.monotonic()
    try:
        run = subprocess.run(["gp", "-q", "-s", "1000000000"], input=program,
            capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as error:
        stdout = error.stdout or b""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        return {"mode": mode, "status": "bounded_search_timeout" if "PHASE|search" in stdout else "bounded_reduction_timeout",
                "wall_seconds": time.monotonic()-started, "stdout": stdout}
    answer = {"mode": mode, "wall_seconds": time.monotonic()-started,
              "integral_coefficient_size": ms.size([g[i] for i in range(5)])}
    if run.returncode or "***" in run.stderr or "DONE" not in run.stdout.splitlines():
        return {**answer, "status": "backend_failure", "stdout": run.stdout, "stderr": run.stderr}
    points, poles, markers = [], 0, {}
    for line in run.stdout.splitlines():
        if line.startswith("POINT|"):
            _, x, y = line.split("|")
            p = ms.quartic_point_to_source(record, x, y, curve, q, map_to)
            if p is None:
                poles += 1
            else:
                points.append(ms.point_record(p))
        elif line == "POLE":
            poles += 1
        elif "|" in line:
            key, value = line.split("|", 1)
            markers[key] = value
    return {**answer, "status": "bounded_search_complete",
            "reduction_milliseconds": int(markers["REDMS"]),
            "search_milliseconds": int(markers["SEARCHMS"]),
            "signed_affine_points_reported": int(markers["COUNT"]),
            "points_transported_to_raw": points, "known_infinity_poles": poles}


def benchmark(payload, candidates, output):
    # Freeze all choices before the first search: median candidate by its
    # worst selected chart, then median chart, within each presentation.
    selected = []
    for presentation in sorted({r["presentation_id"] for r in payload["candidates"]}):
        rows = sorted((r for r in payload["candidates"] if r["presentation_id"] == presentation),
            key=lambda r: (max(c["selected"]["maximum_bits"] for c in r["charts"]), r["candidate_id"]))
        row = rows[len(rows)//2]
        charts = sorted(row["charts"], key=lambda c: (c["selected"]["maximum_bits"], c["mask"]))
        selected.append((row, charts[len(charts)//2]))
    result = {"schema": "mw16.short-model-tiny-benchmark.v1", "status": "RUNNING",
        "input_sha256": digest(OUTPUT), "height_bound": 10000, "timeout_seconds": 2,
        "selection": "per presentation: median candidate by worst selected size, then median chart by selected size; ties by identifier/mask",
        "frozen_charts": [{"candidate_id": r["candidate_id"], "mask": c["mask"]} for r, c in selected],
        "results": [], "claim_boundary": "Different coordinates define different bounded boxes. No quotient rank is inferred from point counts or bounded misses."}
    write(output, result)
    lookup = {c["candidate_id"]: c for c in candidates}
    for row, chart in selected:
        trials = [search_chart(lookup[row["candidate_id"]], row, chart, mode, 10000, 2)
                  for mode in ("raw_direct", "short_direct", "selected_direct", "selected_reduced")]
        result["results"].append({"candidate_id": row["candidate_id"], "mask": chart["mask"], "trials": trials})
        write(output, result)
        print(f"BENCHMARK|{row['candidate_id']}|{[t['status'] for t in trials]}", flush=True)
    result["status"] = "COMPLETE_TINY_ARITHMETIC_BENCHMARK"
    result["status_counts"] = dict(Counter(t["status"] for r in result["results"] for t in r["trials"]))
    write(output, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="recompute deterministic transport/size certificate; do not rerun benchmark")
    parser.add_argument("--benchmark-only", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--maximum-candidates", type=int, default=104, help="bounded local pilot; full certificate requires all 104")
    args = parser.parse_args()
    if not 1 <= args.maximum_candidates <= 104:
        parser.error("candidate limit must be 1..104")
    inputs, previous = read(INPUT), read(MASKS)
    if inputs["status"] != "PASS_TARGET_FREE_A1_MW16_PARAMETER_CANDIDATES":
        raise ArithmeticError("candidate source is not passing")
    candidates = inputs["candidates"][:args.maximum_candidates]
    mask_lookup = {r["candidate_id"]: r["deepest_masks"] for r in previous["results"]}
    hashes = {relative(p): digest(p) for p in (INPUT, MASKS, LEGACY, Path(__file__),
        CAS / "mw16_model_size.py", CAS / "mod2_reduction_independence.py")}
    if args.benchmark_only:
        payload = read(OUTPUT)
        if payload["inputs"] != hashes:
            raise ArithmeticError("prepared certificate is stale")
        benchmark(payload, candidates, BENCHMARK)
        return
    legacy = SourceFileLoader("mw16_size_legacy", str(LEGACY)).load_module()
    payload = {"schema": "mw16.short-model-transport.v1", "status": "PASS_EXACT_MODEL_SIZE_AUDIT",
        "inputs": hashes, "software": {"sage": version(), "pari_gp": subprocess.run(["gp", "--version-short"], capture_output=True, text=True, check=True).stdout.strip()},
        "selection_rule": "short Weierstrass model: minimum maximum/total coefficient bits; chart: minimum maximum/total integral coefficient bits over frozen 31-transform menu; names break ties",
        "limits": {"candidate_count": len(candidates), "minimal_model_timeout_seconds": 2, "prime_bound_independence": 1000,
                   "new_parameter_sweep": False, "adaptive_wave": False, "selmer_calls": 0},
        "claim_boundary": "Exact model and section transport plus rank-16 independence only. No new direction, upper rank bound, saturation or negative rank evidence.",
        "candidates": []}
    for i, c in enumerate(candidates):
        checkpoint = LOCAL / (c["candidate_id"]+".json")
        if checkpoint.exists() and not args.check:
            cached = read(checkpoint)
            if cached["inputs"] != hashes:
                raise ArithmeticError(f"stale checkpoint: {checkpoint}")
            row = cached["result"]
        else:
            row = prepare(c, mask_lookup[c["candidate_id"]], legacy)
            if not args.check:
                write(checkpoint, {"inputs": hashes, "result": row})
        payload["candidates"].append(row)
        print(f"MODEL_SIZE|{i+1}/{len(candidates)}|{c['candidate_id']}|quartic_bits={max(t['selected']['maximum_bits'] for t in row['charts'])}", flush=True)
    if len(candidates) != 104:
        write(LOCAL / "pilot.json.gz", payload)
        return
    if args.check:
        if payload != read(OUTPUT):
            raise ArithmeticError("deterministic certificate differs")
        print("PASS_EXACT_REPLAY_ALL_104_MODELS_1664_SECTIONS_856_CHARTS")
        return
    write(OUTPUT, payload)
    rows = payload["candidates"]
    def extent(values):
        values = sorted(values)
        return {"minimum": values[0], "median": values[len(values)//2], "maximum": values[-1]}
    summary = {"schema": payload["schema"], "status": payload["status"], "certificate": relative(OUTPUT), "certificate_sha256": digest(OUTPUT),
        "candidate_count": len(rows), "section_count": 16*len(rows), "chart_count": sum(len(r["charts"]) for r in rows),
        "raw_weierstrass_bits": extent(next(t for t in r["model_trials"] if t["name"] == "raw")["model_size"]["maximum_bits"] for r in rows),
        "minimal_weierstrass_bits": extent(ms.size(r["global_minimal_model"])["maximum_bits"] for r in rows),
        "raw_quartic_bits": extent(c["raw_quartic_size"]["maximum_bits"] for r in rows for c in r["charts"]),
        "selected_quartic_bits": extent(c["selected"]["maximum_bits"] for r in rows for c in r["charts"]),
        "universal_integral_quartic_lower_bound_bits": extent(r["any_integral_binary_quartic_coefficient_bit_lower_bound"] for r in rows),
        "selected_transform_counts": dict(Counter(c["selected"]["name"] for r in rows for c in r["charts"])),
        "claim_boundary": payload["claim_boundary"]}
    write(SUMMARY, summary)
    if not args.prepare_only:
        benchmark(payload, candidates, BENCHMARK)


if __name__ == "__main__":
    main()
