"""Pure helpers for the Curve302 closure-structure experiments."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
import math
import os
from pathlib import Path

GENERIC_RANK = 17
DIM = 14
FULL = (1 << DIM) - 1
SCALE = 1_000_000
EXPECTED_DIRECTIONS = (
    "recovered-local-01", "recovered-local-02", "recovered-local-03", "recovered-local-04",
    "recovered-strict-01", "recovered-strict-02", "recovered-strict-03",
    "residual-strict-01", "residual-strict-02", "residual-strict-03", "residual-strict-04",
    "residual-strict-05", "residual-strict-06", "residual-strict-07",
)


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def atomic(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("." + path.name + ".tmp-" + str(os.getpid()))
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bit_indices(mask, dimension=DIM):
    return [i for i in range(dimension) if mask >> i & 1]


def mask_names(mask, names):
    return [names[i] for i in bit_indices(mask, len(names))]


def validate_landscape(source, dimension=DIM):
    require(source.get("status") == "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE", "landscape is not passed")
    names = source["fixed_basis"]["direction_ids"]
    require(len(names) == dimension and len(set(names)) == dimension, "wrong direction roster")
    states = source["subset_states"]
    require(len(states) == 1 << dimension, "incomplete subset landscape")
    for i, row in enumerate(states):
        require(row["state_mask"] == i, "subset landscape reordered")
        require(len(row["retained_numerators"]) == dimension, "wrong retained-cost width")
        for j, value in enumerate(row["retained_numerators"]):
            require((value is None) == bool(i >> j & 1), "retained edge nullity mismatch")
    for mask in range(1 << dimension):
        for added in range(dimension):
            if mask >> added & 1:
                continue
            child = mask | (1 << added)
            for target in range(dimension):
                if child >> target & 1:
                    continue
                require(int(states[child]["retained_numerators"][target]) <= int(states[mask]["retained_numerators"][target]),
                        "retained costs are not monotone")
    return names, states


def closure_mask(states, start, threshold, dimension=DIM):
    mask = int(start)
    while True:
        costs = states[mask]["retained_numerators"]
        add = [i for i in range(dimension) if not (mask >> i & 1) and int(costs[i]) <= int(threshold)]
        if not add:
            return mask
        for i in add:
            mask |= 1 << i


def closure_waves(states, start, threshold, dimension=DIM):
    mask = int(start)
    waves = []
    while True:
        costs = states[mask]["retained_numerators"]
        add = [i for i in range(dimension) if not (mask >> i & 1) and int(costs[i]) <= int(threshold)]
        if not add:
            return mask, waves
        waves.append(add)
        for i in add:
            mask |= 1 << i


def arrival_thresholds(states, start, dimension=DIM):
    size = 1 << dimension
    full = size - 1
    inf = 1 << 250
    best = [inf] * size
    best[start] = 0
    arrivals = [0 if start >> i & 1 else inf for i in range(dimension)]
    for pop in range(dimension + 1):
        for mask in range(size):
            if mask.bit_count() != pop or (mask & start) != start or best[mask] == inf:
                continue
            for target in range(dimension):
                if mask >> target & 1:
                    continue
                child = mask | (1 << target)
                candidate = max(best[mask], int(states[mask]["retained_numerators"][target]))
                if candidate < best[child]:
                    best[child] = candidate
                if candidate < arrivals[target]:
                    arrivals[target] = candidate
    require(all(x < inf for x in arrivals) and best[full] < inf, "unreachable minimax state")
    return arrivals, best[full]


def inclusion_minimal(masks):
    winners = []
    for mask in sorted(set(masks), key=lambda m: (m.bit_count(), m)):
        if not any((prior & mask) == prior for prior in winners):
            winners.append(mask)
    return winners


def minimal_triggers(states, target, threshold, dimension=DIM):
    candidates = [mask for mask in range(1 << dimension)
                  if not (mask >> target & 1) and int(states[mask]["retained_numerators"][target]) <= int(threshold)]
    return inclusion_minimal(candidates)


def closure_law_audit(closures, dimension=DIM):
    size = 1 << dimension
    extensive = idempotent = monotone = True
    for mask in range(size):
        c = closures[mask]
        extensive &= (c & mask) == mask
        idempotent &= closures[c] == c
        for i in range(dimension):
            if not (mask >> i & 1):
                monotone &= (c & ~closures[mask | (1 << i)]) == 0
    require(extensive and idempotent and monotone, "threshold operator failed closure axioms")
    antecedents = exchange_bad = anti_bad = 0
    exchange_examples = []
    anti_examples = []
    for closed in range(size):
        if closures[closed] != closed:
            continue
        outside = [i for i in range(dimension) if not (closed >> i & 1)]
        for x in outside:
            cx = closures[closed | (1 << x)]
            for y in outside:
                if x == y or not (cx >> y & 1):
                    continue
                antecedents += 1
                reverse = bool(closures[closed | (1 << y)] >> x & 1)
                if not reverse:
                    exchange_bad += 1
                    if len(exchange_examples) < 8:
                        exchange_examples.append([closed, x, y])
                if reverse:
                    anti_bad += 1
                    if len(anti_examples) < 8:
                        anti_examples.append([closed, x, y])
    return {
        "closure_axioms": {"extensive": extensive, "idempotent": idempotent, "monotone": monotone},
        "exchange_antecedents": antecedents,
        "matroid_exchange_violations": exchange_bad,
        "anti_exchange_violations": anti_bad,
        "matroid_like": antecedents > 0 and exchange_bad == 0,
        "antimatroid_like_anti_exchange": antecedents > 0 and anti_bad == 0,
        "exchange_violation_examples": exchange_examples,
        "anti_exchange_violation_examples": anti_examples,
    }


def rref_binary(rows, width=DIM):
    work = [int(r) & ((1 << width) - 1) for r in rows if int(r)]
    pivot = 0
    for column in range(width):
        bit = 1 << column
        k = next((i for i in range(pivot, len(work)) if work[i] & bit), None)
        if k is None:
            continue
        work[pivot], work[k] = work[k], work[pivot]
        for i in range(len(work)):
            if i != pivot and work[i] & bit:
                work[i] ^= work[pivot]
        pivot += 1
        if pivot == len(work):
            break
    return tuple(sorted(work[:pivot]))


def binary_rank(rows, width=DIM):
    return len(rref_binary(rows, width))


def vector_mask_mod2(values):
    mask = 0
    for i, value in enumerate(values):
        if int(value) & 1:
            mask |= 1 << i
    return mask


def primitive_vector(values):
    values = [int(v) for v in values]
    g = 0
    for value in values:
        g = math.gcd(g, abs(value))
    if g:
        values = [v // g for v in values]
    first = next((v for v in values if v), 0)
    if first < 0:
        values = [-v for v in values]
    return tuple(values)


def is_strict_mod2(mask, strict_signature):
    return binary_rank([*strict_signature, mask]) == len(strict_signature)


def stable_float(value):
    value = float(value)
    require(math.isfinite(value), "nonfinite diagnostic")
    return format(value, ".12g")


def spearman(xs, ys):
    require(len(xs) == len(ys) and len(xs) > 1, "bad Spearman vectors")
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i + 1
            while j < len(order) and values[order[j]] == values[order[i]]:
                j += 1
            rank = (i + j - 1) / 2.0
            for k in range(i, j):
                out[order[k]] = rank
            i = j
        return out
    a, b = ranks(xs), ranks(ys)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    numerator = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = sum((x - ma) ** 2 for x in a)
    db = sum((y - mb) ** 2 for y in b)
    return 0.0 if not da or not db else numerator / math.sqrt(da * db)


def trigger_report(states, names, threshold):
    answer = {}
    for target in range(len(names)):
        masks = minimal_triggers(states, target, threshold, len(names))
        answer[names[target]] = {
            "count": len(masks),
            "size_histogram": dict(sorted(Counter(m.bit_count() for m in masks).items())),
            "masks": masks,
            "named": [mask_names(mask, names) for mask in masks],
        }
    return answer
