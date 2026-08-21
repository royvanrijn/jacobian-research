#!/usr/bin/env sage
"""Discovery beam for a rootless positive frame in a fixed lattice genus.

Every retained edge records the prime and primitive p-divisible vector, so a
hit can be replayed exactly with ``find_p_neighbor_from_vec``.  This is a
bounded discovery search, not a proof that a rootless class does or does not
exist when no hit is found.
"""

from sage.all import *
from pathlib import Path
import argparse
import csv
import hashlib
import random


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def quadratic_form(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2 if row == column else gram[row, column]
            )
    form = QuadraticForm(ZZ, gram.nrows(), coefficients)
    assert form.Hessian_matrix() == gram
    return form


def reduced_key(gram):
    try:
        reduced = gram.LLL_gram()
    except Exception:
        reduced = gram
    payload = repr(tuple(tuple(map(int, row)) for row in reduced.rows()))
    return hashlib.sha256(payload.encode()).hexdigest()


def root_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return 0, 0
    roots = matrix(ZZ, result[2]).transpose()
    return roots.rank(), count


def save_gram(gram, path, metadata):
    lines = [f"# {key} = {value}" for key, value in metadata.items()]
    lines.extend(" ".join(map(str, row)) for row in gram.rows())
    Path(path).write_text("\n".join(lines) + "\n")


parser = argparse.ArgumentParser()
parser.add_argument("--gram", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--generations", type=int, default=8)
parser.add_argument("--beam", type=int, default=24)
parser.add_argument("--samples-per-parent", type=int, default=80)
parser.add_argument("--primes", default="3,5,7,11")
parser.add_argument("--seed", type=int, default=20260821)
args = parser.parse_args()

random.seed(args.seed)
set_random_seed(args.seed)
primes = [ZZ(value) for value in args.primes.split(",")]
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=False)

start_gram = load_gram(args.gram.resolve())
assert start_gram.is_square() and start_gram.is_positive_definite()
assert all(start_gram[i, i] % 2 == 0 for i in range(start_gram.nrows()))
determinant = start_gram.det()
assert all(determinant % prime for prime in primes)

start_form = quadratic_form(start_gram)
start_roots = root_data(start_gram)
start_key = reduced_key(start_gram)
start_path = out / "frame-000000.txt"
save_gram(start_gram, start_path, {
    "node": 0, "parent": "", "prime": "", "vector": "",
    "root_rank": start_roots[0], "roots": start_roots[1],
})

# Entries are (root_rank, root_count, form, key, file, node_id).
frontier = [(start_roots[0], start_roots[1], start_form,
             start_key, start_path, 0)]
seen = {start_key}
manifest = []
next_node = 1

print(
    f"ROOTLESSP|stage=start|rank={start_gram.nrows()}|det={determinant}"
    f"|root_rank={start_roots[0]}|roots={start_roots[1]}"
    f"|primes={','.join(map(str, primes))}", flush=True,
)

for generation in range(1, args.generations + 1):
    candidates = list(frontier)
    generation_best = None
    for parent_rank, parent_count, parent, parent_key, parent_path, parent_id in frontier:
        for sample in range(args.samples_per_parent):
            prime = random.choice(primes)
            try:
                vector_p = parent.find_primitive_p_divisible_vector__random(prime)
                child = parent.find_p_neighbor_from_vec(prime, vector_p)
                gram = child.Hessian_matrix()
                if gram.det() != determinant or not all(
                    gram[i, i] % 2 == 0 for i in range(gram.nrows())
                ):
                    continue
                key = reduced_key(gram)
                if key in seen:
                    continue
                seen.add(key)
                rank, count = root_data(gram)
            except Exception:
                continue

            node_id = next_node
            next_node += 1
            path = out / f"frame-{node_id:06d}.txt"
            metadata = {
                "node": node_id,
                "parent": parent_id,
                "prime": prime,
                "vector": tuple(vector_p),
                "root_rank": rank,
                "roots": count,
            }
            save_gram(gram, path, metadata)
            manifest.append({
                "node": node_id,
                "parent": parent_id,
                "generation": generation,
                "prime": prime,
                "vector": repr(tuple(vector_p)),
                "root_rank": rank,
                "roots": count,
                "frame": path.name,
            })
            item = (rank, count, child, key, path, node_id)
            candidates.append(item)
            if generation_best is None or item[:2] < generation_best[:2]:
                generation_best = item

            if rank == 0:
                rootless_path = out / "ROOTLESS.txt"
                save_gram(gram, rootless_path, metadata)
                print(
                    f"ROOTLESSP|stage=hit|generation={generation}|node={node_id}"
                    f"|parent={parent_id}|prime={prime}|file={rootless_path}",
                    flush=True,
                )
                with (out / "edges.tsv").open("w") as handle:
                    writer = csv.DictWriter(
                        handle, fieldnames=manifest[0].keys(), delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(manifest)
                raise SystemExit(0)

    # Root rank is primary.  Preserve several representatives per exact root
    # fingerprint so the beam does not collapse onto one isometry class.
    candidates.sort(key=lambda item: (item[0], item[1], item[3]))
    next_frontier = []
    per_fingerprint = {}
    for item in candidates:
        fingerprint = item[:2]
        if per_fingerprint.get(fingerprint, 0) >= 3:
            continue
        per_fingerprint[fingerprint] = per_fingerprint.get(fingerprint, 0) + 1
        next_frontier.append(item)
        if len(next_frontier) == args.beam:
            break
    frontier = next_frontier
    best = frontier[0]
    print(
        f"ROOTLESSP|stage=generation|generation={generation}"
        f"|new={next_node - 1 - sum(1 for row in manifest if row['generation'] < generation)}"
        f"|seen={len(seen)}|best_root_rank={best[0]}|best_roots={best[1]}"
        f"|best_node={best[5]}", flush=True,
    )

with (out / "edges.tsv").open("w") as handle:
    fieldnames = (manifest[0].keys() if manifest else
                  ("node", "parent", "generation", "prime", "vector",
                   "root_rank", "roots", "frame"))
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest)

print(
    f"ROOTLESSP|stage=done|seen={len(seen)}"
    f"|best_root_rank={frontier[0][0]}|best_roots={frontier[0][1]}"
    f"|best_node={frontier[0][5]}", flush=True,
)
