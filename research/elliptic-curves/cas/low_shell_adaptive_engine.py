"""Target-blind, intrinsic low-shell worker shared by cascade controllers.

The parent module supplies ``protocol()``, ``D``, ``SEED`` and ``ROW``.  This
module intentionally depends only on that frozen wave state; it never reads a
catalogue, a prior point campaign, or a residual-direction diagnostic.
"""

from __future__ import annotations

from pathlib import Path

import retained_native19_trial_v3 as base


control = None


def bind(module) -> None:
    global control
    control = module
    base.ROW, base.D, base.SEED = module.ROW, module.D, module.SEED


def sources() -> dict[str, str]:
    root = Path(__file__).resolve().parents[2] if control is None else control.ROOT
    return {
        **base.sources(),
        str(Path(__file__).relative_to(root)): base.cert.hashed(Path(__file__)),
    }


def _initial(cache):
    if control is None:
        raise ArithmeticError("low-shell engine is not bound")
    base.ROW, base.D, base.SEED = control.ROW, control.D, control.SEED
    return base.initial(cache)


def _expected(maps: dict) -> list[dict]:
    # This is the only ordering choice: the smallest nonzero sampled intrinsic
    # parity norms, then parity mask.  The map builder records all 2048 rows.
    return sorted(maps["sample"], key=lambda row: (row["metric_norm"], row["parity"]))[:49]


def worker() -> None:
    p = control.protocol()
    maps, out = control.cert.read(control.D / "maps.json"), control.D / "result.json"
    cache = base.ReductionCache(base.MemoryFactStore())
    seed, state = _initial(cache)
    model = tuple(map(control.cert.F, seed["curve"]))
    if out.exists():
        raise FileExistsError("preserve low-shell point attempt")
    if maps["status"] != "COMPLETE_DECLARED_MAPS" or maps["protocol_hash"] != base.digest(p) or len(maps["rows"]) != 49:
        raise ArithmeticError("frozen low-shell maps are incomplete")
    if maps["centres"] != _expected(maps):
        raise ArithmeticError("maps are not the intrinsic low-shell roster")
    data = {key: seed[key] for key in ("family", "parameter", "curve", "generic_points")}
    data.update(protocol_hash=base.digest(p), maps_sha256=control.cert.hashed(control.D / "maps.json"), initial_state=state.record(), initial_dimension=control.ROW["initial_rank"], centres=maps["centres"], metric_gram=maps["metric_gram"], charts=[], status="RUNNING", rank_lower_bound=state.rank, final_state=state.record(), arithmetic_facts=cache.store.snapshot())
    control.checkpoint(out, data)
    for index, mapping in enumerate(maps["rows"]):
        state, archive = base.rotate(state)
        archive_path = control.D / "states" / f"{index:03}.json"
        control.checkpoint(archive_path, archive)
        representative = mapping["centre"]["representative"] + [0] * (state.rank - control.ROW["initial_rank"])
        search = base.PointedQuarticSearch(state=state, centre={"coefficients": representative}, coordinate_policy=mapping["coordinate_policy"])
        transcript, points = base.backend.execute(search, mapping, p["height"], p["seconds_per_chart"], p["gp_sha256"])
        compression = base.compress(model, state.basis, representative, points)
        for point_index in compression["kept_indices"]:
            state = state.adjoin(points[point_index], cache=cache)
        final = state.record()
        data["charts"].append({"index": index, "centre": mapping["centre"], "archive_path": str(archive_path.relative_to(control.ROOT)), "archive_sha256": control.cert.hashed(archive_path), "search": transcript, "admission_compression": compression, "admission_observations": final["state"]["observations"], "state_key": state.key, "rank_lower_bound": state.rank})
        data.update(final_state=final, rank_lower_bound=state.rank, arithmetic_facts=cache.store.snapshot())
        control.checkpoint(out, data)
        print("LOW-SHELL CASCADE", index + 1, transcript["status"], "rank", state.rank, flush=True)
    data["status"] = "COMPLETE_DECLARED_ADAPTIVE_ATTEMPT"
    control.checkpoint(out, data)


def replay() -> None:
    p = control.protocol()
    data, maps = control.cert.read(control.D / "result.json"), control.cert.read(control.D / "maps.json")
    cache = base.ReductionCache(base.MemoryFactStore())
    seed, state = _initial(cache)
    model = tuple(map(control.cert.F, seed["curve"]))
    if data["protocol_hash"] != base.digest(p) or data["maps_sha256"] != control.cert.hashed(control.D / "maps.json") or data["initial_state"] != state.record():
        raise ArithmeticError("low-shell worker is detached from its frozen state")
    if maps["centres"] != _expected(maps) or [row["centre"] for row in maps["rows"]] != _expected(maps):
        raise ArithmeticError("low-shell ordering does not replay")
    if [row["parity"] for row in maps["sample"]] != control.masks(p):
        raise ArithmeticError("low-shell parity sample differs")
    for index, row in enumerate(data["charts"]):
        mapping = maps["rows"][index]
        state, archive = base.rotate(state)
        archive_path = control.ROOT / row["archive_path"]
        if row["index"] != index or row["centre"] != mapping["centre"] or control.cert.hashed(archive_path) != row["archive_sha256"] or control.cert.read(archive_path) != archive:
            raise ArithmeticError("low-shell archived state differs")
        representative = mapping["centre"]["representative"] + [0] * (state.rank - control.ROW["initial_rank"])
        search = base.PointedQuarticSearch(state=state, centre={"coefficients": representative}, coordinate_policy=mapping["coordinate_policy"])
        transcript = row["search"]
        points = base.backend.replay(search, mapping, transcript)
        compression = base.compress(model, state.basis, representative, points)
        if compression != row["admission_compression"]:
            raise ArithmeticError("low-shell admission compression differs")
        for point_index in compression["kept_indices"]:
            state = state.adjoin(points[point_index], cache=cache)
        if state.key != row["state_key"] or state.rank != row["rank_lower_bound"]:
            raise ArithmeticError("low-shell state transition differs")
    if data["status"] != "COMPLETE_DECLARED_ADAPTIVE_ATTEMPT" or len(data["charts"]) != 49 or state.record() != data["final_state"]:
        raise ArithmeticError("low-shell final state differs")
    base.checked_rank(model, state.basis, state.reductions.primes, state.no_two_torsion_prime)
    print("REPLAYED LOW-SHELL CASCADE", len(data["charts"]), "rank >=", state.rank, flush=True)
