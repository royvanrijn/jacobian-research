#!/usr/bin/env python3
"""Bounded shared arithmetic, local-feature and known-subspace descent commands.

Discovery writes portable arithmetic facts and exact witnesses. --verify replays
them, including into an empty cache, without conic/local/point searches. The
frozen blinded BNF experiment is deliberately a separate programme.
"""
from dataclasses import asdict
import argparse
import json
from pathlib import Path
import shutil
import sys

CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))
from research_runtime.store import FactStore, checkpoint, digest
from research_runtime.supervisor import Limits, run


def worker(request, store, limits, *, retained=None):
    from sage.all import pari
    from research_runtime.arithmetic import ArithmeticContext, CurveModel, TwoTorsionContext, rationals
    from research_runtime.finite_reduction import ReductionCache
    from research_runtime.mw_state import MWState
    from research_runtime.sage_arithmetic import SageArithmetic
    from research_runtime.sage_subspace import PointClassWitness, SageSubspaceBackend
    from research_runtime.subspace import SubspaceDescent
    pari.allocatemem(min(64_000_000, limits.pari_stack_bytes), limits.pari_stack_bytes, silent=True)
    if retained is not None:
        if retained["request_hash"] != digest(request):
            raise ValueError("witness request mismatch")
        store.import_snapshot(retained["arithmetic_facts"])
    arithmetic, reductions = SageArithmetic(store), ReductionCache(store)
    mode = request.get("mode", "features")
    if mode not in ("features", "subspace", "complete-selmer"):
        raise ValueError("unknown arithmetic pipeline mode")

    def prepare(spec):
        if "two_torsion" in spec:
            return arithmetic.prepare_congruent(spec["model"], TwoTorsionContext(**spec["two_torsion"]),
                spec["generator"], factor_primes=spec.get("factor_primes", ()), discover=retained is None)
        return arithmetic.prepare(spec["model"], factor_primes=spec.get("factor_primes", ()),
                                  discover=retained is None)

    def point_state(context, spec):
        state = MWState.empty(context, cache=reductions, primes=spec["good_primes"],
                              no_two_torsion_prime=spec["no_two_torsion_prime"])
        for point in spec.get("points", []):
            previous = state.rank
            state = state.adjoin(point, cache=reductions, extra_primes=spec.get("extra_primes", ()))
            if state.rank != previous+1:
                raise ArithmeticError("supplied global basis is ambiguous at the declared primes")
        return state

    result = {"schema": "elliptic-curves.arithmetic-pipeline.v1", "request": request,
              "request_hash": digest(request), "limits": asdict(limits), "mode": mode}
    if mode == "features":
        # Scheduling needs neither global minimalization nor the complete
        # discriminant factorization. Derive the raw division polynomial by
        # exact formulas, then ask only for the selected small local features.
        from sage.all import PolynomialRing, QQ, matrix
        spec = request["curve"]
        model = CurveModel(spec["model"])
        algebra = (TwoTorsionContext(**spec["two_torsion"]) if "two_torsion" in spec else
                   TwoTorsionContext(model.two_division_polynomial))
        if "two_torsion" in spec:
            ring = PolynomialRing(QQ, "theta")
            K = ring.quotient(ring(list(map(QQ, algebra.polynomial))), "a")
            alpha = 4*K(list(map(QQ, spec["generator"])))
            if ring(list(map(QQ, model.two_division_polynomial)))(alpha) != 0 or not matrix(QQ, 3, 3,
                lambda i, j: (alpha**j).lift()[i]).det():
                raise ArithmeticError("feature field identification does not match the curve")
        features = arithmetic.fast_features(algebra, primes=request.get("feature_primes", [2, 3, 5, 7]))
        result.update({"curve_model": list(model.coefficients),
                       "two_torsion_context": {"polynomial": list(algebra.polynomial), "labels": list(algebra.labels)},
                       "features": features, "status": "LOCAL_FEATURES_ONLY"})
        if retained is not None and retained["features"] != features:
            raise ArithmeticError("retained local features do not replay")
    elif mode == "complete-selmer":
        context = prepare(request["curve"])
        requirement = request["requirement"]
        field_support=request.get('field_factor_primes',
            (2,*context.bad_primes) if context.two_torsion.polynomial==context.minimal_model.two_division_polynomial else ())
        arithmetic.field(context.two_torsion, factor_primes=field_support, discover=retained is None)
        selmer = arithmetic.full_selmer(context, requirement=requirement, discover=retained is None)
        result.update({"context": context.record(), "selmer": selmer,
                       "status": "RETAINED_COMPLETE_SELMER"})
        if retained is not None:
            # This mode validates retained setup/result integrity; it is not a
            # substitute for an independent upper-bound certificate replayer.
            if retained["selmer"] != selmer:
                raise ArithmeticError("retained complete-Selmer result changed")
            result["status"] = "RETAINED_RESULT_INTEGRITY_ONLY"
    else:
        if retained is None:
            source = prepare(request["source"])
            context = prepare(request["target"]) if "target" in request else source
            state = point_state(source, request["source"])
        else:
            state = MWState.from_record(retained["source_state"], cache=reductions)
            context = ArithmeticContext.from_record(retained["context"])
            # Recheck input binding, not just the certificate's internal hashes.
            if state.arithmetic != prepare(request["source"]) or context != (
                prepare(request["target"]) if "target" in request else state.arithmetic):
                raise ArithmeticError("retained state belongs to another requested model")
            if state.basis != tuple(rationals(point) for point in request["source"].get("points", [])):
                raise ArithmeticError("retained source point basis changed")
        arithmetic.field(context.two_torsion, factor_primes=request.get("field_factor_primes", ()),
                         discover=retained is None)
        global_witness = PointClassWitness(state, reductions)
        backend = SageSubspaceBackend(arithmetic, context, global_witness,
                                      local_candidate_cap=request.get("local_candidate_cap", 40000))
        descent = SubspaceDescent(context, global_witness.classes, backend)
        witness = descent.run(retained=None if retained is None else retained["descent"],
                              include_ct=not request.get("local_only", False))
        result.update({"context": context.record(), "source_state": state.record(),
                       "descent": witness, "status": witness["status"]})
        if request.get("search_masks"):
            if retained is None:
                spec = request.get("target", request["source"])
                initial = point_state(context, spec)
                pruning=[]
                final = descent.search(initial, witness, iter(request["search_masks"]), limits=limits,
                    pruning_audit=pruning,
                    search=lambda s, c, g, m, lim: backend.search_cover(s, c, g, m, lim,
                        cache=reductions, height=request["point_height"]))
                result.update({"initial_search_state": initial.record(), "final_search_state": final.record(),
                               "search_pruning":pruning,
                               "bounded_search_completeness_replayed": False})
            else:
                initial = MWState.from_record(retained["initial_search_state"], cache=reductions)
                final = MWState.from_record(retained["final_search_state"], cache=reductions)
                if initial.arithmetic != context or final.arithmetic != context or final.basis[:initial.rank] != initial.basis:
                    raise ArithmeticError("invalid retained point-search state transition")
                result.update({"initial_search_state": initial.record(), "final_search_state": final.record(),
                               "search_pruning":retained.get('search_pruning',[]),
                               "bounded_search_completeness_replayed": False})
    result["arithmetic_facts"] = store.snapshot()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--verify", type=Path, help="replay retained arithmetic and mathematical witnesses")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--wall-seconds", type=float, default=60)
    parser.add_argument("--rss-bytes", type=int, default=1_073_741_824)
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument("--sage", default=shutil.which("sage") or "sage")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if bool(args.request) == bool(args.verify):
        parser.error("choose exactly one of --request and --verify")
    if args.output.exists():
        parser.error("output already exists; preserve it and choose a new path")
    limits = Limits(args.wall_seconds, args.rss_bytes, pari_stack_bytes=args.pari_stack_bytes)
    if args.worker:
        retained = json.loads(args.verify.read_text()) if args.verify else None
        request = retained["request"] if retained else json.loads(args.request.read_text())
        result = worker(request, FactStore(args.cache_dir), limits, retained=retained)
        checkpoint(args.output, result)
        print(result["status"], flush=True)
        return
    command = [args.sage, "-python", str(Path(__file__).resolve()), "--worker",
               "--verify" if args.verify else "--request", str((args.verify or args.request).resolve()),
               "--output", str(args.output.resolve()), "--cache-dir", str(args.cache_dir.resolve()),
               "--wall-seconds", str(args.wall_seconds), "--rss-bytes", str(args.rss_bytes),
               "--pari-stack-bytes", str(args.pari_stack_bytes)]
    outcome = run(command, limits=limits, log_path=args.output.with_suffix(".log"),
                  result_path=args.output, checkpoint_path=args.output.with_suffix(".supervisor.json"))
    print(json.dumps({"outcome": outcome["outcome"], "output": str(args.output)}, sort_keys=True))
    if outcome["outcome"] != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
