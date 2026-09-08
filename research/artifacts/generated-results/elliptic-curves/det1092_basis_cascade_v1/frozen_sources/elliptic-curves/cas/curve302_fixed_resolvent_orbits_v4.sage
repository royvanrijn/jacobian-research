#!/usr/bin/env sage -python
"""v4 replay: v3 determinant slice with a Sage-10.9 squareclass extractor."""

from __future__ import annotations

import argparse
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path

from sage.all import NumberField, QQ


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
OUTPUT = ART / "curve302_fixed_resolvent_orbits_v4.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-fixed-resolvent-orbits-v4"
V3 = HERE.with_name("curve302_fixed_resolvent_orbits_v3.sage")

spec = importlib.util.spec_from_loader("curve302_fixed_resolvent_orbits_v3", SourceFileLoader("curve302_fixed_resolvent_orbits_v3", str(V3)))
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)


def read(path): return json.loads(path.read_text())
def put_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True); handle.write("\n")


def extract_squareclass_sage109(poly, raw):
    """Exact resolvent-field identification in the Sage version pinned by v3."""
    resolvent = v3.v2.base.classical_resolvent(poly)
    if not resolvent.is_irreducible():
        return {"status": "REJECTED_REDUCIBLE_CLASSICAL_RESOLVENT"}
    Kc, Kr = NumberField(resolvent, "u"), NumberField(raw, "z")
    outcome = Kc.is_isomorphic(Kr, isomorphism_maps=True)
    if not outcome or outcome[0] is not True:
        return {"status": "UNKNOWN_CUBIC_FIELD_IDENTIFICATION"}
    # For the norm-square quartic chart, the resolvent roots are
    # alpha_i - Tr(alpha)/2.  Its u^2 coefficient is therefore Tr(alpha)/2.
    trace_alpha = QQ(resolvent[2])
    # Sage 10.9 returns a singleton list containing the image of the
    # source generator, rather than a callable homomorphism.
    alpha = Kr(outcome[1][0]) + trace_alpha
    return {"status": "PASS", "alpha_raw_power_basis": [QQ(alpha.polynomial()[i]) for i in range(3)],
            "classical_resolvent_ascending": [str(value) for value in resolvent.list()],
            "identification": "Exact NumberField(resolvent)-to-NumberField(raw) Sage isomorphism."}


def compute():
    # The determinant search and every later gate are unchanged.  This only
    # replaces the Sage API call on the first field-bearing branch.
    v3.v2.base.extract_squareclass = extract_squareclass_sage109
    result = v3.compute()
    result["schema"] = "elliptic-curves.curve302-fixed-resolvent-orbits.v4"
    result["extractor_runtime_fix"] = "Use NumberField(polynomial, name), the Sage-10.9 API, in the dormant classical-resolvent branch."
    result["bindings"] = {v3.v2.base.relative(path): v3.v2.base.digest(path) for path in (v3.PROTOCOL, V3, v3.V2, v3.v2.V1, v3.v2.base.MAXIMUM_FORM, v3.v2.base.STRICT, v3.v2.base.FILTRATION, v3.v2.base.CURVE, v3.v2.base.HIDDEN_SOLVER, HERE)}
    return result


def capture():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    WORK.mkdir(parents=True, exist_ok=True)
    put_new(OUTPUT, compute())
    print(read(OUTPUT)["status"], flush=True)


def check():
    stored=read(OUTPUT)
    assert stored["bindings"] == {v3.v2.base.relative(path): v3.v2.base.digest(path) for path in (v3.PROTOCOL, V3, v3.V2, v3.v2.V1, v3.v2.base.MAXIMUM_FORM, v3.v2.base.STRICT, v3.v2.base.FILTRATION, v3.v2.base.CURVE, v3.v2.base.HIDDEN_SOLVER, HERE)}
    if stored["status"] != "PASS_NEW_STRICT_CHARACTER_HANDED_OFF": assert compute() == stored
    print("PASS v4 fixed-resolvent replay", flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=("capture","check"));args=parser.parse_args()
    (capture if args.mode=="capture" else check)()
