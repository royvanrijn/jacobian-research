#!/usr/bin/env python3
"""Checkpointed determinant-closure shards for seven-term Keller collisions.

All raw coordinate degrees <= D are allowed.  Shard p,q,r fixes the smallest
pure-x exponent of F1 and the two smallest pure-x exponents of F2.  F3 remains
unrestricted.  The residual simultaneous 2/3 swap makes these shards a cover,
rather than a disjoint cover, of all collisions.  A finite search, including
an empty balanced-support list, is a claim only when status == complete.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEGACY = HERE.parents[1] / "foundations"
sys.path.insert(0, str(LEGACY))
from jcsearch.low_degree_census import (
    Support, SupportEnumerator, ZERO, determinant_buckets, coefficient_system,
    singular_batch_status, groebner_status,
)

SCHEMA = "sparse-seven-shard.v1"

def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".research/.scratch")
    temporary.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)

def digest():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

class Shard:
    def __init__(self, degree, bound, p, q, r, compiler=None):
        if not (2 <= p <= degree and 2 <= q < r <= degree):
            raise ValueError("invalid p,q,r")
        self.degree, self.bound, self.seed = degree, bound, (p,q,r)
        self.compiler = compiler or SupportEnumerator(degree, bound)
        self.atoms = self.compiler.atoms
        self.bits = tuple(1 << i for i in range(len(self.atoms)))
        allowed = []
        for i, (row, exponent) in enumerate(self.atoms):
            if exponent[1:] == (0,0):
                n = exponent[0]
                if row == 0 and n < p:
                    continue
                if row == 1 and n < r and n != q:
                    continue
            allowed.append(i)
        self.allowed = sum(self.bits[i] for i in allowed)
        self.seed_mask = sum(
            self.bits[self.compiler.atom_index[atom]]
            for atom in ((0,(p,0,0)), (1,(q,0,0)), (1,(r,0,0)))
        )
        self.axis_masks = tuple(
            sum(self.bits[i] for i in indices) for indices in self.compiler.axis_atoms
        )
        self.completions = {}
        for exponent, entries in self.compiler.completions.items():
            if exponent == ZERO:
                continue
            masks = {
                sum(self.bits[i] for i in required)
                for required, _ in entries
                if all(self.bits[i] & self.allowed for i in required)
            }
            self.completions[exponent] = tuple(sorted(masks))

    def support(self, mask):
        return Support.from_atoms(self.atoms[i] for i in self.indices(mask))

    @staticmethod
    def indices(mask):
        while mask:
            bit = mask & -mask
            yield bit.bit_length() - 1
            mask ^= bit

    def children(self, mask):
        size = mask.bit_count()
        if size > self.bound:
            return False, ()
        # Repair every prospective singleton on the collision axis before
        # scanning determinant completions.  This necessary lower bound is safe
        # even if a determinant repair introduces further pure-axis terms.
        offenders = [
            row for row in (1,2) if (mask & self.axis_masks[row]).bit_count() == 1
        ]
        if size + len(offenders) > self.bound:
            return False, ()
        support = self.support(mask)
        singles = [
            bucket.exponent for bucket in determinant_buckets(support)
            if bucket.exponent != ZERO and len(bucket.terms) == 1
        ]
        if singles:
            if size == self.bound:
                return False, ()
            best = None
            for exponent in singles:
                options = {
                    mask | completion for completion in self.completions[exponent]
                    if (mask | completion).bit_count() <= self.bound
                    and completion & ~mask
                }
                if not options:
                    return False, ()
                if best is None or len(options) < len(best):
                    best = options
            return False, tuple(sorted(best))
        if offenders:
            row = offenders[0]
            candidates = self.axis_masks[row] & self.allowed & ~mask
            return False, tuple(mask | self.bits[i] for i in self.indices(candidates))
        assert support.collision_axis_unimodularity_possible()
        if size == self.bound:
            return True, ()
        return True, tuple(mask | self.bits[i] for i in self.indices(self.allowed & ~mask))

    def run(self, path, seconds=180, node_limit=500000, progress_every=30):
        source_hash = digest()
        if path.exists():
            state = json.loads(path.read_text())
            assert state["schema"] == SCHEMA
            assert state["source_sha256"] == source_hash, "source changed; use a new checkpoint"
            assert state["parameters"] == dict(degree=self.degree, bound=self.bound, seed=list(self.seed))
            if state["status"] == "complete":
                return state
            stack = [int(s,16) for s in state["stack"]]
            seen = {int(s,16) for s in state["seen"]}
            balanced = {int(s,16) for s in state["balanced_masks"]}
        else:
            state = dict(
                schema=SCHEMA, source_sha256=source_hash,
                legacy_sha256=hashlib.sha256((LEGACY/"jcsearch/low_degree_census.py").read_bytes()).hexdigest(),
                parameters=dict(degree=self.degree,bound=self.bound,seed=list(self.seed)),
                elapsed_seconds=0.0,
            )
            stack, seen, balanced = [self.seed_mask], set(), set()
        start = last = time.monotonic()
        initial_nodes = len(seen)

        def save():
            state.update(
                status="complete" if not stack else "incomplete",
                visited=len(seen), pending=len(stack),
                balanced_masks=[hex(s) for s in sorted(balanced)],
                supports=[self.support(s).to_json() for s in sorted(balanced)],
                stack=[hex(s) for s in stack],
                seen=[hex(s) for s in sorted(seen)] if stack else [],
                elapsed_seconds=state["elapsed_seconds"]+time.monotonic()-start,
            )
            write_json(path, state)
            print(json.dumps({k:state[k] for k in ("parameters","status","visited","pending","elapsed_seconds")} |
                             dict(balanced=len(balanced))), flush=True)

        while stack:
            if time.monotonic()-start >= seconds or len(seen)-initial_nodes >= node_limit:
                break
            mask = stack.pop()
            if mask in seen:
                continue
            seen.add(mask)
            keep, children = self.children(mask)
            if keep:
                balanced.add(mask)
            stack.extend(c for c in reversed(children) if c not in seen)
            if time.monotonic()-last > progress_every:
                # A live checkpoint does not constitute completed enumeration.
                before = state.get("elapsed_seconds",0)
                save()
                state["elapsed_seconds"] = before
                last = time.monotonic()
        save()
        return state

def read_support(rows):
    return Support(tuple(tuple(tuple(e) for e in row) for row in rows))

def solve(path, independent=False):
    state = json.loads(path.read_text())
    assert state["status"] == "complete", "do not certify an unfinished shard"
    supports = sorted({read_support(rows).canonical_under_collision_stabilizer()
                       for rows in state["supports"]})
    results = singular_batch_status(supports, timeout_seconds=300) if supports else ()
    payload = dict(
        schema="sparse-seven-exact.v1", input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        parameters=state["parameters"], support_count=len(state["supports"]),
        orbit_count=len(supports), results=results,
        all_unit=all(row["unit_ideal"] for row in results),
        independent_sympy=False,
    )
    if independent:
        replay = [groebner_status(s) for s in supports]
        assert [row["unit_ideal"] for row in replay] == [row["unit_ideal"] for row in results]
        payload.update(independent_sympy=True, sympy_results=replay)
    write_json(path.with_suffix(".exact.json"), payload)
    print(json.dumps({k:v for k,v in payload.items() if k not in ("results","sympy_results")}),flush=True)
    return payload

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree",type=int,default=7)
    parser.add_argument("--bound",type=int,default=7)
    parser.add_argument("--seed",default="2,2,3")
    parser.add_argument("--seconds",type=int,default=180)
    parser.add_argument("--nodes",type=int,default=500000)
    parser.add_argument("--output",type=Path,default=HERE/"results/p2-q2-r3.json")
    parser.add_argument("--solve",action="store_true")
    parser.add_argument("--independent",action="store_true")
    args=parser.parse_args()
    if args.solve:
        solve(args.output,args.independent)
    else:
        p,q,r=map(int,args.seed.split(","))
        shard=Shard(args.degree,args.bound,p,q,r)
        shard.run(args.output,args.seconds,args.nodes)

if __name__ == "__main__":
    main()

