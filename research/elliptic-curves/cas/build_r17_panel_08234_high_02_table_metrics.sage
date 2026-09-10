#!/usr/bin/env sage-python
"""Certify the minimal-model table metrics for r17-panel-08234-high-02."""

import argparse
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, QQ, ZZ, RealField, pari
from sage.version import version


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
import certify_compact_r17_candidates as cert

CURVE_ID = "r17-panel-08234-high-02"
PANEL_ID = "08234-high-02"
PANEL = ROOT / "artifacts/generated-results/elliptic-curves/r17_60_panel_results_v1.json"
CONDUCTORS = ROOT / "artifacts/generated-results/elliptic-curves/inventory321_conductor_snapshot_v4.json"
OUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_panel_08234_high_02_table_metrics_v1.json"


def run(check=False):
    panel = cert.read(PANEL)
    selected = next(row for row in panel["results"] if row["id"] == PANEL_ID)
    packet = selected["packet"]
    conductor_record = cert.read(CONDUCTORS)["records"][CURVE_ID]["result"]
    conductor = conductor_record["certificate"]
    if (selected["status"] != "PASS_INDEPENDENT_R17_60_CASE"
            or selected["rank_lower_bound"] != 24
            or packet["rank_lower_bound"] != 24
            or conductor_record["status"] != "PASS_INDEPENDENT_CONDUCTOR_REPLAY"
            or conductor["status"] != "EXACT"
            or conductor["curve"] != packet["curve"]):
        raise ArithmeticError("rank/conductor source gate failed")

    source = EllipticCurve(QQ, [QQ(value) for value in packet["curve"]])
    pari_raw, pari_change = pari.ellinit(list(source.a_invariants())).ellminimalmodel()
    model = [QQ(pari_raw[i]) for i in range(5)]
    minimal = EllipticCurve(QQ, model)
    if list(minimal.global_minimal_model().a_invariants()) != model:
        raise ArithmeticError("Sage and PARI minimal models differ")
    iso = source.isomorphism_to(minimal)
    if list(iso.tuple()) != [QQ(value) for value in pari_change]:
        raise ArithmeticError("Sage and PARI changes of variables differ")

    c4, c6 = minimal.c_invariants()
    discriminant = minimal.discriminant()
    if (c4, c6, discriminant) != tuple(QQ(pari_raw[i]) for i in (9, 10, 11)):
        raise ArithmeticError("exact invariant replay differs")

    invariant_gcd = ZZ(c4).gcd(ZZ(c6)).abs()
    gcd_primes = invariant_gcd.prime_divisors()
    local_minimality = []
    for prime in sorted(set(gcd_primes + [ZZ(2), ZZ(3)])):
        data = minimal.local_data(prime, algorithm="generic", proof=True)
        if data.discriminant_valuation() != discriminant.valuation(prime):
            raise ArithmeticError("local minimality replay failed")
        local_minimality.append({
            "prime": str(prime),
            "minimal_discriminant_valuation": int(data.discriminant_valuation()),
        })

    conductor_local = conductor["local_data"]
    conductor_primes = [ZZ(row["prime"]) for row in conductor_local]
    saved_minimal_discriminant = ZZ(1)
    saved_conductor = ZZ(1)
    local_replay = []
    for row, prime in zip(conductor_local, conductor_primes):
        data = minimal.local_data(prime, algorithm="generic", proof=True)
        exponent = int(data.conductor_valuation())
        valuation = int(data.discriminant_valuation())
        if (exponent != int(row["conductor_exponent"])
                or valuation != int(row["minimal_discriminant_valuation"])):
            raise ArithmeticError("minimal local data differs from conductor certificate")
        saved_minimal_discriminant *= prime ** valuation
        saved_conductor *= prime ** exponent
        local_replay.append({
            "prime": str(prime),
            "conductor_exponent": exponent,
            "minimal_discriminant_valuation": valuation,
        })
    if (saved_minimal_discriminant != abs(discriminant)
            or saved_conductor != ZZ(conductor["exact_conductor"])):
        raise ArithmeticError("exact local products differ")

    transported = []
    for pair in packet["points"]:
        point = source([QQ(value) for value in pair])
        image = iso(point)
        if image.is_zero() or image not in minimal:
            raise ArithmeticError("point transport failed")
        transported.append([str(image[0]), str(image[1])])

    real = RealField(160)
    naive = real(max(abs(c4) ** 3, c6 ** 2)).log()
    faltings = real(minimal.faltings_height(stable=False, prec=160))
    faltings_low = real(minimal.faltings_height(stable=False, prec=96))
    if abs(faltings - faltings_low) > real("1e-22"):
        raise ArithmeticError("Faltings-height precision replay failed")

    result = {
        "schema": "elliptic-curves.selected-table-metrics.v1",
        "status": "PASS_EXACT_MINIMAL_MODEL_AND_POINT_TRANSPORT",
        "id": CURVE_ID,
        "rank_lower_bound": 24,
        "source_ainvs": [str(value) for value in source.a_invariants()],
        "minimal_ainvs": [str(value) for value in model],
        "original_to_minimal_isomorphism": [str(value) for value in iso.tuple()],
        "c4": str(c4),
        "c6": str(c6),
        "discriminant": str(discriminant),
        "log_abs_discriminant": float(real(abs(discriminant)).log()),
        "naive_height": float(naive),
        "faltings_height": float(faltings),
        "invariant_gcd": str(invariant_gcd),
        "minimality_local_checks": local_minimality,
        "conductor": str(saved_conductor),
        "conductor_local_replay": local_replay,
        "points": transported,
        "point_transports": len(transported),
        "sage_version": version,
        "precision_bits": [96, 160],
        "definitions": {
            "log": "natural logarithm",
            "naive_height": "log max(abs(c4)^3,c6^2) on the global minimal model",
            "faltings_height": "-1/2 log(period lattice area) on the global minimal model; Sage stable=False",
        },
        "sources": {
            str(path.relative_to(ROOT)): cert.hashed(path)
            for path in (Path(__file__), PANEL, CONDUCTORS)
        },
        "claim_boundary": "Exact global minimal equation, discriminant, source-model isomorphism, local conductor replay, and 24 point transports. Heights and logarithms are numerical approximations checked at 96 and 160 bits. The saved independent rank lower bound and exact conductor are reused; no exact-rank, novelty, or record claim.",
    }
    if check:
        if cert.read(OUT) != result:
            raise ArithmeticError("selected table-metric replay differs")
    elif OUT.exists():
        if cert.read(OUT) != result:
            raise FileExistsError("preserve selected table-metric certificate")
    else:
        cert.write(OUT, result)
    print(
        "R17_PANEL_TABLE_METRICS_PASS"
        f"|id={CURVE_ID}|rank=24|points={len(transported)}"
        f"|logN={real(saved_conductor).log():.12f}"
        f"|naive={naive:.12f}|faltings={faltings:.12f}"
        f"|logDelta={real(abs(discriminant)).log():.12f}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    run(parser.parse_args().check)
