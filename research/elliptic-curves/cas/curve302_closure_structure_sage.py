"""Sage-only helpers for Curve302 closure-structure experiments."""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import importlib.machinery
import json
import math
from pathlib import Path

from curve302_closure_structure_core import (
    DIM, GENERIC_RANK, SCALE, EXPECTED_DIRECTIONS, atomic, binary_rank,
    is_strict_mod2, primitive_vector, read, require, rref_binary,
    spearman, stable_float, vector_mask_mod2,
)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves"
LANDSCAPE = ART / "curve302_exceptional_subgroup_landscape_v1.json"
VISIBILITY = ART / "curve302_residual_visibility_geometry_v1.json"
M24 = ART / "curve302_recovered_followup_wave_03_mod2_v1.json"
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
FILTRATION = ART / "curve302_recovered_quotient_local_filtration_v1.json"
GEOMETRY_SOURCE = CAS / "prospective_half_lattice_v3.sage"
TWO_SEED_LOCAL = LOCAL / "curve302-seeded-v3-amplifier-v1"
PANEL_LOCAL = LOCAL / "curve302-seed-universality-panel-v1"


def load_module(path, name):
    return importlib.machinery.SourceFileLoader(name, str(path)).load_module()


def seed_folder(seed):
    return (TWO_SEED_LOCAL if seed in {"recovered-strict-02", "recovered-strict-03"} else PANEL_LOCAL) / seed


def exact_linear_combination(model, points, coefficients):
    """Exact Sage group law; deliberately independent of the pointed-sieve C++ worker."""
    from sage.all import EllipticCurve, QQ, ZZ
    require(len(points) == len(coefficients), "point/coefficient length mismatch")
    E = EllipticCurve(QQ, list(map(QQ, model)))
    total = E(0)
    for point, coefficient in zip(points, coefficients):
        require(int(coefficient) == coefficient, "nonintegral exact group coefficient")
        total += ZZ(coefficient) * E([QQ(str(point[0])), QQ(str(point[1]))])
    if total.is_zero():
        return None
    x, y = total.xy()
    return F(str(x)), F(str(y))


def fixed_basis_context():
    from sage.all import GF, QQ, RealField, ZZ, matrix, vector
    import icarm_curve302 as curve
    geometry = load_module(GEOMETRY_SOURCE, "closure_lab_geometry")
    visibility, m24, parent = read(VISIBILITY), read(M24), read(PARENT)
    require(visibility["status"] == "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC", "visibility input not passed")
    model = tuple(F(v) for v in m24["curve"])
    base17 = tuple(tuple(F(v) for v in p) for p in m24["independent_points"][:GENERIC_RANK])
    directions = visibility["directions"]
    names = [row["id"] for row in directions]
    require(tuple(names) == EXPECTED_DIRECTIONS, "fixed basis roster changed")

    public = tuple(tuple(F(v) for v in point) for point in curve.SHORT_POINTS)
    require(len(public) == 31, "public rank-31 point roster changed")
    generic_public = matrix(ZZ, parent["basis_embedding_in_public_D"])
    require(generic_public.dimensions() == (31, 17), "generic public embedding changed")
    for column in range(GENERIC_RANK):
        recovered = exact_linear_combination(model, public, generic_public.column(column))
        require(recovered == base17[column], "generic public embedding no longer matches M17 basis")

    target_public = []
    targets = []
    for entry in directions:
        word = vector(ZZ, entry["public_word"])
        require(len(word) == 31, "target public word has wrong width")
        target = exact_linear_combination(model, public, word)
        require(target is not None, "target public word became infinity")
        for row in entry["exact_pointed_quartic_rows"]:
            translated = tuple(F(v) for v in row["exact_target_point_short_model"])
            inverse = exact_linear_combination(model, base17, [-int(value) for value in row["target_translation_m17_word"]])
            reconstructed = translated if inverse is None else exact_linear_combination(model, (translated, inverse), (1, 1))
            require(reconstructed == target, "target chart and public word encode different directions")
        target_public.append(word)
        targets.append(target)

    fixed_to_public = generic_public.augment(matrix(ZZ, [list(v) for v in target_public]).transpose())
    require(fixed_to_public.dimensions() == (31, 31) and abs(fixed_to_public.det()) == 1,
            "fixed diagnostic basis is not unimodular in public D")
    ambient_points = (*base17, *targets)

    height, asymmetry = geometry.canonical_height_gram(model, ambient_points)
    ambient = matrix(ZZ, [[int((value * 1_000_000).to_integral_value()) for value in row] for row in height])
    require(ambient.is_positive_definite(), "rounded 31-dimensional metric is not positive definite")
    expected_hash = read(LANDSCAPE)["fixed_basis"]["rounded_metric_sha256"]
    actual_hash = sha256(json.dumps([list(map(int, row)) for row in ambient.rows()], separators=(",", ":")).encode()).hexdigest()
    require(actual_hash == expected_hash, "reconstructed rounded metric differs from landscape")

    strict_public = read(FILTRATION)["local_filtration_mod_2"]["strict_kernel_public_words"]
    inv2 = fixed_to_public.change_ring(GF(2)).inverse()
    strict_q_rows = []
    for word in strict_public:
        fixed = inv2 * vector(GF(2), word)
        strict_q_rows.append(vector_mask_mod2([int(fixed[GENERIC_RANK+i]) for i in range(DIM)]))
    strict_signature = rref_binary(strict_q_rows)
    require(len(strict_signature) == 10, "strict quotient dimension is not ten")

    real = RealField(192)
    high = matrix(real, [[real(str(value).strip("()")) for value in row] for row in height])
    return {
        "geometry": geometry, "model": model,
        "base17": base17, "targets": tuple(targets), "ambient_points": ambient_points,
        "ambient": ambient, "asymmetry": str(asymmetry), "names": names,
        "strict_signature": strict_signature, "height_real": high,
    }


def qnorm(form, values):
    from sage.all import vector
    v = vector(form.base_ring(), values)
    return v * form * v


def nearest_integer(value):
    lo = int(math.floor(float(value)))
    return min(range(lo-2, lo+4), key=lambda k: abs(value-k))


def relation_experiment(folder):
    from sage.all import QQ
    ctx = fixed_basis_context()
    ambient = ctx["ambient"].change_ring(QQ)
    names = ctx["names"]
    A = ambient[:GENERIC_RANK, :GENERIC_RANK]
    B = ambient[:GENERIC_RANK, GENERIC_RANK:]
    C = ambient[GENERIC_RANK:, GENERIC_RANK:]
    quotient = C - B.transpose() * A.inverse() * B
    require(quotient.is_positive_definite(), "quotient Schur form is not positive definite")

    landscape = read(LANDSCAPE)
    states = landscape["subset_states"]
    predicted, observed, pair_rows = [], [], []
    for i in range(DIM):
        for j in range(DIM):
            if i == j:
                continue
            optimum = 2 * quotient[i,j] / quotient[i,i]
            k0 = nearest_integer(optimum)
            value, k = min(
                ((4*quotient[j,j]-4*k*quotient[i,j]+k*k*quotient[i,i], k) for k in range(k0-3, k0+4)),
                key=lambda pair: (pair[0], abs(pair[1]), pair[1]),
            )
            before = int(states[0]["retained_numerators"][j])
            after = int(states[1 << i]["retained_numerators"][j])
            pdrop = float(4*quotient[j,j] - value)
            odrop = before - after
            predicted.append(pdrop); observed.append(odrop)
            pair_rows.append({
                "added": names[i], "target": names[j], "best_k": k,
                "projected_before": str(4*quotient[j,j]), "projected_after": str(value),
                "projected_drop": stable_float(pdrop), "retained_drop": odrop,
            })

    candidates = {}
    def add(indices, choices):
        for coeffs in product(choices, repeat=len(indices)):
            values = [0] * DIM
            for index, coefficient in zip(indices, coeffs):
                values[index] = coefficient
            values = primitive_vector(values)
            if any(values):
                candidates[values] = qnorm(quotient, values)
    for size in (2, 3):
        for indices in combinations(range(DIM), size):
            add(indices, (-2, -1, 1, 2))
    for indices in combinations(range(DIM), 4):
        add(indices, (-1, 1))
    shortest = sorted(candidates.items(), key=lambda row: (row[1], row[0]))[:512]
    cohorts = landscape["fixed_basis"]["cohorts"]
    schema_counts = Counter()
    short_rows = []
    for vector, norm in shortest:
        support = [i for i, value in enumerate(vector) if value]
        schema = (
            tuple(sorted(abs(vector[i]) for i in support)),
            tuple(sorted("local" if "local" in cohorts[names[i]] else "strict" for i in support)),
        )
        schema_counts[str(schema)] += 1
        short_rows.append({
            "vector": list(vector),
            "symbolic": " + ".join(f"{vector[i]}*{names[i]}" for i in support),
            "support": [names[i] for i in support],
            "norm": str(norm), "scaled_norm": stable_float(float(norm)/SCALE),
            "mod2_strict": is_strict_mod2(vector_mask_mod2(vector), ctx["strict_signature"]),
            "schema": str(schema),
        })

    closure = read(Path(folder)/"closure-laws.json")
    trigger_relations = []
    for threshold_name, block in closure["thresholds"].items():
        for target_name, trigger_block in block["minimal_triggers"].items():
            target = names.index(target_name)
            for trigger_names in trigger_block["named"]:
                if len(trigger_names) > 4:
                    continue
                indices = [names.index(name) for name in trigger_names]
                best = None
                for coeffs in product(range(-2,3), repeat=len(indices)):
                    vector = [0] * DIM
                    vector[target] = 2
                    for index, coefficient in zip(indices, coeffs):
                        vector[index] -= coefficient
                    norm = qnorm(quotient, vector)
                    key = (norm, sum(abs(x) for x in coeffs), coeffs)
                    if best is None or key < best[0]:
                        best = (key, coeffs, vector)
                trigger_relations.append({
                    "threshold": threshold_name, "target": target_name, "trigger": trigger_names,
                    "best_coefficients": list(best[1]), "near_relation_vector": best[2],
                    "quotient_norm": str(best[0][0]), "scaled_norm": stable_float(float(best[0][0])/SCALE),
                })

    result = {
        "schema": "curve302-quotient-near-relations.v1",
        "status": "PASS_QUOTIENT_RELATION_ANALYSIS",
        "direction_ids": names,
        "rounded_metric_scale": SCALE,
        "schur_quotient": [[str(value) for value in row] for row in quotient.rows()],
        "strict_quotient_rref_bitmasks": list(ctx["strict_signature"]),
        "pairwise_projected_vs_retained": {"spearman": stable_float(spearman(predicted, observed)), "rows": pair_rows},
        "short_primitive_combinations": short_rows,
        "short_schema_counts": dict(schema_counts.most_common()),
        "trigger_near_relations": trigger_relations,
        "boundary": (
            "The Schur form is exact for the entrywise-1e6-rounded 31-point height metric, not exact canonical heights. "
            "Every nonzero short vector is a near-relation, not a Mordell-Weil dependency."
        ),
    }
    atomic(Path(folder)/"quotient-relations.json", result)
    print("CLOSURE_LAB_EXP2|spearman={}|shortest={}".format(
        result["pairwise_projected_vs_retained"]["spearman"], short_rows[0]["symbolic"]), flush=True)


def recognize_point(ctx, point):
    from sage.all import matrix
    geometry, model = ctx["geometry"], ctx["model"]
    ambient_points, high = ctx["ambient_points"], ctx["height_real"]
    real = high.base_ring()
    pairing = geometry.canonical_height_gram(model, (*ambient_points, point))[0]
    column = matrix(real, 31, 1, [real(str(pairing[i][31]).strip("()")) for i in range(31)])
    coefficients = high.solve_right(column)
    for denominator in range(1, 65):
        word = [int(round(float(denominator*coefficients[i,0]))) for i in range(31)]
        error = max(abs(denominator*coefficients[i,0]-word[i]) for i in range(31))
        if error > real("1e-24"):
            continue
        if exact_linear_combination(model, (point,), (denominator,)) == exact_linear_combination(model, ambient_points, word):
            common = denominator
            for value in word:
                common = math.gcd(common, abs(value))
            if common > 1:
                denominator //= common
                word = [value // common for value in word]
            return {"denominator": denominator, "word": word, "max_numeric_error": str(error)}
    return None


def stage_audit_path(folder, stage):
    return folder / "replay-M17" / f"epoch-{int(stage['epoch']):02d}" / stage["audit"]


def trajectory_experiment(folder):
    ctx = fixed_basis_context()
    names, base17 = ctx["names"], ctx["base17"]
    runs = []
    convergence = defaultdict(Counter)
    transition_counter = Counter()
    cache = {}

    def recognized(point):
        if point not in cache:
            cache[point] = recognize_point(ctx, point)
        return cache[point]

    for seed in names:
        source = seed_folder(seed)
        seed_input = read(source/"seed-input.json")
        verified = read(source/"seeded-verified.json")
        terminal = read(source/"replay-M17/terminal.json")
        require(verified.get("status") == "PASS_INDEPENDENT_SEEDED_V3_REPLAY", "seed replay not passed: " + seed)
        require(int(verified["initial_rank"]) == 18 and int(terminal["initial_rank"]) == 18, "seed did not start at 18")
        initial = tuple(tuple(F(v) for v in point) for point in seed_input["points"])
        require(initial[:17] == base17 and len(initial) == 18, "seed prefix differs: " + seed)
        seed_word = recognized(initial[17])
        require(seed_word is not None and seed_word["denominator"] == 1, "seed not integral in fixed D: " + seed)
        quotient_seed = seed_word["word"][17:]
        index = names.index(seed)
        require(sum(abs(x) for x in quotient_seed) == 1 and abs(quotient_seed[index]) == 1,
                "seed does not represent its named quotient axis: " + seed)

        quotient_rows = [vector_mask_mod2(quotient_seed)]
        run_signatures = {1: str(rref_binary(quotient_rows))}
        previous = initial
        stages = []
        for stage in terminal["stages"]:
            audit = read(stage_audit_path(source, stage))
            points = tuple(tuple(F(v) for v in point) for point in audit["independent_points"])
            require(points[:len(previous)] == previous, "V3 basis prefix changed: " + seed)
            new_rows = []
            for point in points[len(previous):]:
                record = recognized(point)
                require(record is not None, "gained point outside displayed D rational span: " + seed)
                quotient_word = record["word"][17:]
                row = {
                    "denominator": record["denominator"], "quotient_word": quotient_word,
                    "integral": record["denominator"] == 1,
                    "support": [names[i] for i, value in enumerate(quotient_word) if value],
                    "max_numeric_error": record["max_numeric_error"],
                }
                require(row["integral"], "gained point has nonintegral fixed-D coordinates: " + seed)
                primitive = primitive_vector(quotient_word)
                bitmask = vector_mask_mod2(quotient_word)
                quotient_rows.append(bitmask)
                row.update({
                    "primitive_quotient_word": list(primitive),
                    "symbolic": " + ".join(f"{primitive[i]}*{names[i]}" for i in range(DIM) if primitive[i]),
                    "mod2_bitmask": bitmask,
                    "mod2_strict": is_strict_mod2(bitmask, ctx["strict_signature"]),
                })
                transition_counter[str(primitive)] += 1
                new_rows.append(row)
            qrank = binary_rank(quotient_rows)
            require(qrank == len(points)-GENERIC_RANK, "mod2 quotient rank disagrees with certified rank: " + seed)
            signature = str(rref_binary(quotient_rows))
            run_signatures[qrank] = signature
            stages.append({
                "epoch": stage["epoch"], "before": stage["before"], "after": stage["after"],
                "new": new_rows, "quotient_mod2_rank": qrank, "quotient_subspace_rref": signature,
            })
            previous = points
        for rank, signature in run_signatures.items():
            convergence[rank][signature] += 1
        runs.append({
            "seed": seed, "final_rank": terminal["final_rank_lower_bound"],
            "stop_reason": terminal["stop_reason"], "charts": terminal["charts"], "stages": stages,
        })
        print("CLOSURE_TRAJECTORY|seed={}|final={}|stages={}".format(seed, terminal["final_rank_lower_bound"], len(stages)), flush=True)

    convergence_rows = []
    for rank in sorted(convergence):
        counts = convergence[rank]
        top = counts.most_common(5)
        total = sum(counts.values())
        convergence_rows.append({
            "quotient_rank": rank, "runs_present": total, "distinct_subspaces": len(counts),
            "largest_basin": top[0][1], "largest_basin_fraction": top[0][1]/total,
            "top_subspaces": [{"rref": key, "count": value} for key, value in top],
        })
    result = {
        "schema": "curve302-actual-v3-closure-trajectories.v1",
        "status": "PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED",
        "direction_ids": names,
        "runs": runs,
        "convergence_by_quotient_rank": convergence_rows,
        "repeated_transition_vectors": [{"primitive": key, "count": value} for key, value in transition_counter.most_common()],
        "boundary": (
            "High-precision heights only propose fixed-basis coordinates; exact Sage elliptic group arithmetic certifies every accepted word. "
            "This is retrospective reconciliation of completed searches, not a prospective selector."
        ),
    }
    atomic(Path(folder)/"trajectories.json", result)
    print("CLOSURE_LAB_EXP3|runs=14|repeated_vectors={}".format(sum(value > 1 for value in transition_counter.values())), flush=True)
