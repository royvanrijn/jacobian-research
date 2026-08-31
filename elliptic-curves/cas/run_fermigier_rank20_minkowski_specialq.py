#!/usr/bin/env sage
"""BNF-free full-ideal Minkowski relation collector for a cubic field.

Unlike the older ``a+b*theta+c*theta^2`` box, this samples the actual
degree-one special-q ideal.  For every special ideal it computes a Z-basis,
applies determinant-one archimedean shape twists, LLL-reduces the resulting
three-dimensional Minkowski lattice, and enumerates short combinations of the
reduced basis.  Accepted relations retain their exact principal generators.
The Fermigier cubic remains the default, but ``--field-polynomial-ascending``
and ``--selmer-rational-primes`` make the same audited collector available to
another fixed no-rational-2-torsion elliptic curve.  Both cubic signatures are
handled: a complex embedding contributes its real and imaginary coordinates
to the three-dimensional real Minkowski lattice.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from hashlib import sha256
import json
from itertools import product
from pathlib import Path
import time


PROTOCOL = "R20MINKQ"
REPO_ROOT = Path(__file__).resolve().parents[2]
ELKIES_RANK28_LEDGER = (
    REPO_ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_bad_place_kummer_ledger_v1.json"
)


def insert_row(pivots, row):
    while row:
        pivot = row.bit_length() - 1
        previous = pivots.get(pivot)
        if previous is None:
            pivots[pivot] = row
            return True
        row ^= previous
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor-base-bound", type=int, default=262523)
    parser.add_argument(
        "--field-polynomial-ascending",
        help=(
            "comma-separated integer coefficients c0,c1,c2,1 for the cubic; "
            "defaults to the Fermigier rank-20 2-division cubic"
        ),
    )
    parser.add_argument(
        "--selmer-rational-primes",
        help=(
            "comma-separated rational primes whose ideals are killed in the "
            "S-class quotient; defaults to 2 and the cubic discriminant primes"
        ),
    )
    parser.add_argument(
        "--curve273",
        action="store_true",
        help="use the pinned ICARM curve-273 short cubic and its declared Selmer primes",
    )
    parser.add_argument(
        "--elkies-rank28",
        action="store_true",
        help=(
            "use the published Elkies t=-9529/5471 cubic, all twelve exact "
            "bad rational primes, and their proved factor hints"
        ),
    )
    parser.add_argument("--special-q-min", type=int, default=262524)
    parser.add_argument("--special-q-max", type=int, default=320000)
    parser.add_argument("--max-special-q", type=int, default=10)
    parser.add_argument("--seed-specials", default="")
    parser.add_argument(
        "--special-residue-degree",
        type=int,
        choices=(1, 2),
        default=1,
        help="residue degree of the special prime ideals sampled by the lattice collector",
    )
    parser.add_argument(
        "--special-primes-in-factor-base",
        action="store_true",
        help=(
            "treat selected special ideals as factor-base columns rather than "
            "large-prime vertices; their principal rows are then retained "
            "directly after exact smoothness testing"
        ),
    )
    parser.add_argument(
        "--allow-rational-special-multiples",
        action="store_true",
        help=(
            "diagnostic compatibility mode: retain elements divisible by every "
            "prime above a selected rational special prime; by default these "
            "canonical (p)-multiples are rejected before relation merging"
        ),
    )
    parser.add_argument(
        "--special-ideal-mode",
        choices=("single", "cycle-pairs"),
        default="single",
        help=(
            "single samples one special ideal at a time; cycle-pairs samples "
            "Q_i Q_(i+1) around a special-ideal cycle"
        ),
    )
    parser.add_argument(
        "--pair-cycle-length",
        type=int,
        default=0,
        help="number of special ideals used in cycle-pairs mode",
    )
    parser.add_argument(
        "--large-prime-merge-mode",
        choices=("all-edges", "spanning-forest"),
        default="all-edges",
        help=(
            "all-edges retains every partial edge; spanning-forest retains "
            "only connectivity edges, so a double-special cycle closes around "
            "its intended graph rather than through parallel-edge shortcuts"
        ),
    )
    parser.add_argument("--trial-prime-bound", type=int, default=10000)
    parser.add_argument("--residual-factor-limit", type=int, default=1 << 40)
    parser.add_argument("--large-prime-bound", type=int, default=1 << 70)
    parser.add_argument("--lll-scale", type=int, default=10**30)
    parser.add_argument("--lattice-combination-bound", type=int, default=2)
    parser.add_argument(
        "--shape-twists",
        default="-20:0,-10:0,-5:0,0:0,5:0,10:0,20:0,0:-5,0:5",
        help=(
            "comma-separated base-10 exponent pairs u:v; in signature (3,0) "
            "the weights are 10^u,10^v,10^(-u-v), while in signature (1,1) "
            "u+v sets the determinant-one complex-plane scale"
        ),
    )
    parser.add_argument("--relation-ledger", type=Path, required=True)
    args = parser.parse_args()
    if args.factor_base_bound < 2 or args.trial_prime_bound < 2:
        raise ValueError("factor-base and trial-prime bounds must be at least 2")
    if args.lattice_combination_bound < 1:
        raise ValueError("--lattice-combination-bound must be positive")
    if args.special_ideal_mode == "cycle-pairs" and args.pair_cycle_length < 3:
        raise ValueError("--pair-cycle-length must be at least three in cycle-pairs mode")
    if args.special_residue_degree == 2 and args.seed_specials:
        raise ValueError("--seed-specials is only defined for degree-one special ideals")
    if args.curve273 and args.elkies_rank28:
        raise ValueError("--curve273 and --elkies-rank28 are mutually exclusive")
    if (args.curve273 or args.elkies_rank28) and (
        args.field_polynomial_ascending or args.selmer_rational_primes
    ):
        raise ValueError(
            "curve presets cannot be combined with explicit field or Selmer-prime inputs"
        )

    def parse_integer_csv(text, description):
        try:
            values = [int(item.strip()) for item in text.split(",") if item.strip()]
        except ValueError as exc:
            raise ValueError(f"{description} must be comma-separated integers") from exc
        if not values:
            raise ValueError(f"{description} must not be empty")
        return values

    from sage.all import (
        NumberField,
        ComplexField,
        PolynomialRing,
        QQ,
        RealField,
        ZZ,
        factor,
        matrix,
        pari,
        prime_range,
        vector,
    )

    twists = []
    for item in args.shape_twists.split(","):
        left, right = item.split(":", 1)
        twists.append((int(left), int(right)))
    if not twists:
        raise ValueError("at least one shape twist is required")

    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    preset_selmer_primes = None
    factor_hint_primes = []
    factor_hint_certificate = None
    if args.curve273:
        from analyze_curve273_relation_pool import S_RATIONAL
        from icarm_curve273 import short_coefficients

        short = [ZZ(QQ(value)) for value in short_coefficients()]
        coefficients = [short[4], short[3], short[1], ZZ(1)]
        preset_selmer_primes = sorted(ZZ(value) for value in S_RATIONAL)
    elif args.elkies_rank28:
        certificate = json.loads(ELKIES_RANK28_LEDGER.read_text())
        factors = certificate.get("factorization", [])
        factor_hint_primes = [ZZ(record["prime"]) for record in factors]
        factor_product = ZZ(1)
        for record in factors:
            factor_product *= ZZ(record["prime"]) ** ZZ(record["exponent"])
        if not (
            certificate.get("status")
            == "COMPLETE_ALL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND"
            and certificate.get("parameter") == "-9529/5471"
            and certificate.get("factorization_product_verified") is True
            and certificate.get("factor_primality_proof_completed") is True
            and certificate.get("all_bad_place_blocks_completed") is True
            and len(factor_hint_primes) == 12
            and len(set(factor_hint_primes)) == 12
            and factor_product == ZZ(certificate["descent_cubic_discriminant"])
        ):
            raise ValueError("the exact Elkies rank-28 bad-place ledger is stale")
        coefficients = [
            ZZ(value)
            for value in certificate["descent_cubic_coefficients_ascending"]
        ]
        if len(coefficients) != 4 or coefficients[-1] != 1:
            raise ValueError("the Elkies rank-28 ledger does not define a monic cubic")
        preset_selmer_primes = sorted(factor_hint_primes)
        # NumberField construction asks PARI for a maximal order.  Supplying
        # the already proved polynomial-discriminant support prevents it from
        # refactoring the 168-digit discriminant before relation collection.
        pari.addprimes(factor_hint_primes)
        factor_hint_certificate = {
            "path": str(ELKIES_RANK28_LEDGER.resolve()),
            "sha256": sha256(ELKIES_RANK28_LEDGER.read_bytes()).hexdigest(),
            "all_factors_proved_prime": True,
            "factor_count": len(factor_hint_primes),
            "claim_boundary": (
                "The hints accelerate maximal-order and ideal arithmetic; "
                "they do not certify class-group completeness."
            ),
        }
    elif args.field_polynomial_ascending:
        coefficients = [ZZ(value) for value in parse_integer_csv(
            args.field_polynomial_ascending, "field-polynomial-ascending"
        )]
        if len(coefficients) != 4 or coefficients[-1] != 1:
            raise ValueError("field-polynomial-ascending must give a monic cubic")
    else:
        coefficients = [
            ZZ(167347710468055045100164888198438918505621536951206),
            ZZ(-5750886029903523759416717668139307),
            ZZ(0),
            ZZ(1),
        ]
    polynomial = sum(
        coefficient * x**index for index, coefficient in enumerate(coefficients)
    )
    if not polynomial.is_irreducible():
        raise ValueError("the defining cubic must be irreducible over QQ")
    field = NumberField(polynomial, "theta")
    theta = field.gen()
    if args.selmer_rational_primes:
        selmer_primes = sorted({ZZ(value) for value in parse_integer_csv(
            args.selmer_rational_primes, "selmer-rational-primes"
        )})
        if any(value < 2 or not value.is_prime(proof=False) for value in selmer_primes):
            raise ValueError("selmer-rational-primes must be rational primes")
    elif preset_selmer_primes is not None:
        selmer_primes = preset_selmer_primes
    else:
        selmer_primes = sorted(
            {ZZ(2)} | {ZZ(p) for p, _ in factor(abs(ZZ(polynomial.discriminant())))}
        )

    factor_base = []
    factor_base_by_hnf = {}
    s_columns = []

    def hnf_key(prime):
        return str(prime.pari_hnf())

    def add_factor_base(prime, is_s=False):
        key = hnf_key(prime)
        index = factor_base_by_hnf.get(key)
        if index is None:
            index = len(factor_base)
            factor_base_by_hnf[key] = index
            factor_base.append(prime)
        if is_s and index not in s_columns:
            s_columns.append(index)
        return index

    for rational_prime in selmer_primes:
        for prime in field.primes_above(rational_prime):
            add_factor_base(prime, True)
    rational_factor_base = [ZZ(p) for p in prime_range(2, args.factor_base_bound + 1)]
    for rational_prime in rational_factor_base:
        if rational_prime in selmer_primes:
            continue
        for prime in field.primes_above(rational_prime):
            add_factor_base(prime)

    by_rational_prime = {}
    for rational_prime in sorted(set(rational_factor_base) | set(selmer_primes)):
        by_rational_prime[rational_prime] = [
            (factor_base_by_hnf[hnf_key(prime)], prime)
            for prime in field.primes_above(rational_prime)
            if hnf_key(prime) in factor_base_by_hnf
        ]
    trial_primes = [p for p in rational_factor_base if p <= args.trial_prime_bound]

    real_field = RealField(256)
    complex_field = ComplexField(256)
    complex_roots = list(polynomial.change_ring(complex_field).roots(multiplicities=False))
    real_roots = sorted(root.real() for root in complex_roots if root.imag() == 0)
    if len(real_roots) == 3:
        embedding_signature = "3,0"
        complex_root = None
    elif len(real_roots) == 1:
        embedding_signature = "1,1"
        complex_root = next((root for root in complex_roots if root.imag() > 0), None)
        if complex_root is None:
            raise ArithmeticError("could not select the positive-imaginary cubic embedding")
    else:
        raise ArithmeticError("unexpected cubic archimedean signature")

    def evaluate_at_root(element, root):
        return sum(
            complex_field(coefficient) * root**power
            for power, coefficient in enumerate(element.list())
        )

    sqrt_two = real_field(2).sqrt()

    def minkowski_coordinates(element):
        if embedding_signature == "3,0":
            return [
                real_field(evaluate_at_root(element, root).real())
                for root in real_roots
            ]
        value = evaluate_at_root(element, complex_root)
        return [
            real_field(evaluate_at_root(element, real_roots[0]).real()),
            sqrt_two * real_field(value.real()),
            sqrt_two * real_field(value.imag()),
        ]
    print(
        f"{PROTOCOL}|stage=input|fb_bound={args.factor_base_bound}"
        f"|fb_columns={len(factor_base)}|S_columns={len(s_columns)}"
        f"|trial_primes={len(trial_primes)}|twists={twists}"
        f"|field={coefficients}|selmer_primes={selmer_primes}"
        f"|archimedean_signature={embedding_signature}",
        flush=True,
    )

    specials = []
    seen_specials = set()
    if args.special_residue_degree == 1:
        for item in args.seed_specials.split(","):
            if item.strip():
                q, r = item.split(":", 1)
                pair = (ZZ(q), ZZ(r))
                specials.append(pair)
                seen_specials.add(pair)
    for rational_prime in prime_range(args.special_q_min, args.special_q_max + 1):
        rational_prime = ZZ(rational_prime)
        if args.special_residue_degree == 1:
            identifiers = [
                ZZ(int(residue))
                for residue in polynomial.change_ring(ZZ.quotient(rational_prime)).roots(
                    multiplicities=False
                )
            ]
        else:
            identifiers = [
                str(prime.pari_hnf())
                for prime in field.primes_above(rational_prime)
                if int(prime.residue_class_degree()) == 2
            ]
        for identifier in identifiers:
            pair = (rational_prime, identifier)
            if pair in seen_specials:
                continue
            specials.append(pair)
            seen_specials.add(pair)
            if len(specials) >= args.max_special_q:
                break
        if len(specials) >= args.max_special_q:
            break
    if args.max_special_q and len(specials) > args.max_special_q:
        specials = specials[: args.max_special_q]
    def special_prime(q, identifier):
        for prime in field.primes_above(q):
            if int(prime.residue_class_degree()) != args.special_residue_degree:
                continue
            if args.special_residue_degree == 1:
                if int((theta - identifier).valuation(prime)):
                    return prime
            elif str(prime.pari_hnf()) == identifier:
                return prime
        return None

    special_entries = []
    for q, residue in specials:
        if args.special_primes_in_factor_base and q > args.factor_base_bound:
            raise ValueError("--special-primes-in-factor-base requires every selected q to be in the factor base")
        prime_q = special_prime(q, residue)
        if prime_q is not None:
            special_entries.append((q, residue, prime_q))
    if args.special_ideal_mode == "single":
        special_targets = [(entry,) for entry in special_entries]
    else:
        cycle_length = min(args.pair_cycle_length, len(special_entries))
        if cycle_length < 3:
            raise ValueError("fewer than three usable special ideals were found")
        cycle = special_entries[:cycle_length]
        special_targets = [
            (cycle[index], cycle[(index + 1) % cycle_length])
            for index in range(cycle_length)
        ]
    print(
        f"{PROTOCOL}|stage=specials|mode={args.special_ideal_mode}"
        f"|count={len(special_targets)}|first={specials[:5]}",
        flush=True,
    )

    def exact_row(alpha, used_rational_primes):
        row = 0
        for rational_prime in used_rational_primes:
            for index, prime in by_rational_prime[rational_prime]:
                if int(alpha.valuation(prime)) & 1:
                    row ^= 1 << index
        return row

    def residual_prime_signature(alpha, rational_prime):
        matches = []
        for prime in field.primes_above(rational_prime):
            if int(alpha.valuation(prime)) & 1:
                matches.append(prime)
        if len(matches) != 1:
            return None
        return int(rational_prime), hnf_key(matches[0])

    generators = []
    closed_relations = []
    pivots = {}
    # Quotient the collector's progress metric by the free principal rows
    # (p), as well as S.  The rows themselves are later stored by the
    # canonical-ledger augmentation step; retaining this separate basis here
    # prevents a relation rank plateau inside the rational-prime subspace from
    # looking like S-class progress.
    canonical_pivots = {}
    for rational_prime in sorted(by_rational_prime):
        canonical_row = 0
        for prime, exponent in field.ideal(rational_prime).factor():
            if int(exponent) & 1:
                canonical_row |= 1 << factor_base_by_hnf[hnf_key(prime)]
        insert_row(canonical_pivots, canonical_row)
    canonical_s_pivots = dict(canonical_pivots)
    for column in s_columns:
        insert_row(canonical_s_pivots, 1 << column)
    quotient_pivots = dict(canonical_s_pivots)
    graph = defaultdict(list)
    ROOT = ("ROOT",)
    total_candidates = 0
    total_cycles = 0
    total_partial = 0
    started = time.monotonic()

    def after_s_dimension():
        reduced = dict(pivots)
        for column in s_columns:
            insert_row(reduced, 1 << column)
        return len(factor_base) - len(reduced)

    def after_canonical_s_dimension():
        return len(factor_base) - len(quotient_pivots)

    def path_between(start, end):
        queue = deque([start])
        parent = {start: None}
        edge = {}
        while queue:
            current = queue.popleft()
            if current == end:
                row = 0
                provenance = set()
                while parent[current] is not None:
                    edge_row, generator_index = edge[current]
                    row ^= edge_row
                    if generator_index in provenance:
                        provenance.remove(generator_index)
                    else:
                        provenance.add(generator_index)
                    current = parent[current]
                return row, provenance
            for next_vertex, edge_row, generator_index in graph[current]:
                if next_vertex not in parent:
                    parent[next_vertex] = current
                    edge[next_vertex] = (edge_row, generator_index)
                    queue.append(next_vertex)
        return None, None

    def add_partial(left, right, row, generator_index):
        path, provenance = path_between(left, right)
        if path is None:
            graph[left].append((right, row, generator_index))
            graph[right].append((left, row, generator_index))
            return None, None
        if args.large_prime_merge_mode == "all-edges":
            graph[left].append((right, row, generator_index))
            graph[right].append((left, row, generator_index))
        if generator_index in provenance:
            provenance.remove(generator_index)
        else:
            provenance.add(generator_index)
        return path ^ row, provenance

    def coordinate_key(alpha):
        """Identify the unavoidable associate pair alpha and -alpha.

        Their principal ideals and valuation rows are identical. Keeping both
        makes the large-prime graph close immediately on ``-alpha^2``, an
        explicit global square that cannot improve a mod-two class quotient.
        """

        coordinates = tuple(str(QQ(value)) for value in alpha.list())
        negative = tuple(str(-QQ(value)) for value in alpha.list())
        return min(coordinates, negative)

    for special_index, members in enumerate(special_targets, 1):
        special_ideal = field.ideal(1)
        for _, _, prime_q in members:
            special_ideal *= prime_q
        ideal_basis = tuple(special_ideal.basis())
        embedding = matrix(
            real_field,
            [
                minkowski_coordinates(element)
                for element in ideal_basis
            ],
        )
        samples = {}
        bound = args.lattice_combination_bound
        combinations = [
            coefficients
            for coefficients in product(range(-bound, bound + 1), repeat=3)
            if any(coefficients)
        ]
        for u, v in twists:
            if embedding_signature == "3,0":
                weights = [
                    real_field(10) ** u,
                    real_field(10) ** v,
                    real_field(10) ** (-u - v),
                ]
            else:
                # A complex embedding is one archimedean place.  Scale its
                # real/imaginary plane uniformly and compensate on the real
                # axis, preserving the lattice covolume.
                complex_exponent = u + v
                weights = [
                    real_field(10) ** (-2 * complex_exponent),
                    real_field(10) ** complex_exponent,
                    real_field(10) ** complex_exponent,
                ]
            lattice = matrix(
                ZZ,
                [
                    [ZZ((args.lll_scale * weights[column] * embedding[row, column]).round()) for column in range(3)]
                    for row in range(3)
                ],
            )
            _, transform = lattice.LLL(transformation=True)
            for combination in combinations:
                original = vector(ZZ, combination) * transform
                alpha = sum(
                    (original[index] * ideal_basis[index] for index in range(3)), field(0)
                )
                if alpha:
                    samples.setdefault(coordinate_key(alpha), (alpha, (u, v), tuple(int(item) for item in combination)))

        before_rank = len(pivots)
        before_after_s = after_s_dimension()
        before_after_canonical_s = after_canonical_s_dimension()
        local_candidates = local_cycles = 0
        for alpha, twist, combination in samples.values():
            member_valuations = [
                int(alpha.valuation(prime_q)) for _, _, prime_q in members
            ]
            if any(valuation <= 0 for valuation in member_valuations):
                continue
            if not args.allow_rational_special_multiples:
                selected_rational_primes = {q for q, _, _ in members}
                if any(
                    all(int(alpha.valuation(prime)) > 0 for prime in field.primes_above(q))
                    for q in selected_rational_primes
                ):
                    continue
            norm = abs(ZZ(alpha.norm()))
            if norm == 0:
                continue
            local_candidates += 1
            total_candidates += 1
            cofactor = norm
            if not args.special_primes_in_factor_base:
                for q, _, _ in members:
                    while cofactor % q == 0:
                        cofactor //= q
            used = []
            for rational_prime in trial_primes:
                if cofactor % rational_prime:
                    continue
                used.append(rational_prime)
                while cofactor % rational_prime == 0:
                    cofactor //= rational_prime
                if cofactor == 1:
                    break
            if cofactor != 1 and cofactor <= args.residual_factor_limit:
                residual = ZZ(1)
                for rational_prime, exponent in factor(cofactor):
                    rational_prime = ZZ(rational_prime)
                    if rational_prime <= args.factor_base_bound:
                        used.append(rational_prime)
                    else:
                        residual *= rational_prime**ZZ(exponent)
                cofactor = residual
            row = exact_row(alpha, used)
            generator_index = len(generators)
            generators.append(
                {
                    "power_basis": list(coordinate_key(alpha)),
                    "source_special_ideals": [
                        [int(q), str(identifier)] for q, identifier, _ in members
                    ],
                    "shape_twist": list(twist),
                    "lll_combination": list(combination),
                    "norm": str(norm),
                }
            )
            vertices = []
            if not args.special_primes_in_factor_base:
                for (q, _, prime_q), valuation in zip(members, member_valuations):
                    if valuation & 1:
                        vertices.append((int(q), hnf_key(prime_q)))
            if cofactor != 1:
                if cofactor > args.large_prime_bound:
                    continue
                residual_factorization = list(factor(cofactor))
                if len(residual_factorization) != 1:
                    continue
                residual_rational_prime = ZZ(residual_factorization[0][0])
                signature = residual_prime_signature(alpha, residual_rational_prime)
                if signature is None:
                    continue
                vertices.append(signature)
            if len(vertices) == 0:
                cycle, provenance = row, {generator_index}
            elif len(vertices) == 1:
                total_partial += 1
                cycle, provenance = add_partial(ROOT, vertices[0], row, generator_index)
            elif len(vertices) == 2:
                total_partial += 1
                cycle, provenance = add_partial(vertices[0], vertices[1], row, generator_index)
            else:
                continue
            if cycle is None:
                continue
            local_cycles += 1
            total_cycles += 1
            closed_relations.append(
                {
                    "fb_parity_mask_hex": hex(cycle),
                    "generator_indices": sorted(provenance),
                    "kind": "unit_dependency" if cycle == 0 else "minkowski_lp_cycle",
                    "source": "minkowski_ideal_lattice",
                }
            )
            insert_row(pivots, cycle)
            insert_row(quotient_pivots, cycle)
        print(
            f"{PROTOCOL}|stage=special_done|index={special_index}/{len(special_targets)}"
            f"|special_ideals={[(int(q), str(identifier)) for q, identifier, _ in members]}"
            f"|samples={len(samples)}|candidates={local_candidates}"
            f"|cycles={local_cycles}|rank_gain={len(pivots)-before_rank}"
            f"|afterS_gain={before_after_s-after_s_dimension()}|rank={len(pivots)}"
            f"|afterS={after_s_dimension()}"
            f"|afterCanonicalS_gain={before_after_canonical_s-after_canonical_s_dimension()}"
            f"|afterCanonicalS={after_canonical_s_dimension()}"
            f"|seconds={time.monotonic()-started:.3f}",
            flush=True,
        )

    ledger = {
        "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
        "status": "exact_minkowski_ideal_relations_not_class_group_completion",
        "special_ideal_mode": args.special_ideal_mode,
        "special_residue_degree": args.special_residue_degree,
        "special_primes_in_factor_base": args.special_primes_in_factor_base,
        "allow_rational_special_multiples": args.allow_rational_special_multiples,
        "large_prime_merge_mode": args.large_prime_merge_mode,
        "curve_preset": (
            "icarm-273"
            if args.curve273
            else "elkies-2026-rank28" if args.elkies_rank28 else None
        ),
        "factor_hint_certificate": factor_hint_certificate,
        "field_polynomial": str(polynomial),
        "defining_polynomial_ascending": [str(value) for value in coefficients],
        "field_discriminant": str(field.discriminant()),
        "generator_coordinate_order": ["1", "theta", "theta^2"],
        "factor_base_bound": args.factor_base_bound,
        "factor_base_completion": {
            "all_prime_ideals_above_rational_primes_through": args.factor_base_bound,
            "materialized_complete_factor_base": True,
            "extra_declared_S_rational_primes": [int(p) for p in selmer_primes],
        },
        "selmer_rational_primes": [int(p) for p in selmer_primes],
        "factor_base": [
            {
                "hnf": hnf_key(prime),
                "norm": int(prime.norm()),
                "residue_degree": int(prime.residue_class_degree()),
                "rational_prime": int(prime.smallest_integer()),
            }
            for prime in factor_base
        ],
        "S_columns": s_columns,
        "collection_early_quotient": {
            "baseline": "exact rational principal rows (p) plus S columns",
            "canonical_rational_row_rank": len(canonical_pivots),
            "dimension_after_canonical_rows_and_S": after_canonical_s_dimension(),
            "interpretation": (
                "Collector progress beyond this baseline is measured exactly "
                "over GF(2), but canonical rows must be materialized with "
                "augment_bnf_free_canonical_principal_relations.py before "
                "the ledger itself carries their generators."
            ),
        },
        "generators": generators,
        "closed_relations": closed_relations,
    }
    args.relation_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.relation_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=summary|fb_columns={len(factor_base)}|S_columns={len(s_columns)}"
        f"|candidates={total_candidates}|partial={total_partial}|cycles={total_cycles}"
        f"|relation_rank={len(pivots)}|afterS={after_s_dimension()}"
        f"|afterCanonicalS={after_canonical_s_dimension()}"
        f"|ledger={args.relation_ledger}|seconds={time.monotonic()-started:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
