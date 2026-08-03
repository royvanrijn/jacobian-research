#!/usr/bin/env python3
"""Translation-orbit observability for primitive projected Hall packets.

The prime-power tomography census reduces the projected two-colour model to
primitive Graver binomials.  Unequal scaled-factorial partitions are already
separated at all scales.  For each remaining exact scalar-factorial collision
this script asks a different, Hall-specific question: can the binomial vanish
along an entire Cartesian translation orbit while every coefficient channel
stays nonzero?

For normalized Taylor coefficients ``R_i`` and ``B_i`` put

    D_R(R_i) = (i+1) R_(i+1),   D_B(B_i) = (i+1) B_(i+1).

If ``F`` is the factorially weighted endpoint binomial, its complete
translation orbit vanishes at a point exactly when all elements of the finite
locally-nilpotent derivative closure vanish there.  Singular tests this on the
coefficient torus by adjoining

    1 - T * product(R_i) * product(B_i).

A unit ideal is an exact obstruction: that translation envelope cannot hide
the primitive packet without coefficient-support loss.  A proper ideal gives
an explicit algebraic *promotion-obstruction variety* in this projected
model.  It is not by itself a Hall packet or a GVC(2) counterexample.

The modes have deliberately different logical force:

``operator`` / ``polynomial``
    Translate one colour only.  Polynomial translation is the literal Taylor
    operation in the current Hall target; colour swapping in the projected
    quotient means both one-colour answers are retained.
``diagonal``
    Apply ``D_R + D_B`` with one common translation parameter.
``independent``
    Apply every ``D_R^a D_B^b``.  This is a stronger two-mark observability
    envelope and identifies what extra marking would suffice.

Normaliz supplies the exact Graver basis and Singular supplies the exact
characteristic-zero standard bases.  The radial bound makes this a finite
experiment, not the proof of unrestricted GVC(2); that proof is the separate
Hall-envelope theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
from collections import Counter, deque
from pathlib import Path
from typing import Any

from research_binary_gvc_prime_power_tomography import (
    Relation,
    State,
    colored_graver_basis,
    exact_marked_partition_collision,
    exact_scalar_factorial_collision,
    normaliz_version,
    primitive_normalization,
    relation_digest,
    relation_record,
    relation_support_size,
    state_counts,
)


Exponent = tuple[int, ...]
Polynomial = dict[Exponent, int]

MODES = ("operator", "polynomial", "diagonal", "independent")


def multinomial(side: tuple[int, ...]) -> int:
    answer = math.factorial(sum(side))
    for multiplicity in side:
        answer //= math.factorial(multiplicity)
    return answer


def state_weight(state: State) -> int:
    """The two side multinomials; the radial factorial is fibre-common."""

    return multinomial(state.operator) * multinomial(state.polynomial)


def side_level(side: tuple[int, ...]) -> int:
    return sum(level * multiplicity for level, multiplicity in enumerate(side))


def polynomial_level_transfer(relation: Relation) -> int:
    return side_level(relation.left.polynomial) - side_level(
        relation.right.polynomial
    )


def state_exponent(state: State) -> Exponent:
    return state.operator + state.polynomial


def normalize_polynomial(polynomial: Polynomial) -> Polynomial:
    polynomial = {
        exponent: coefficient
        for exponent, coefficient in polynomial.items()
        if coefficient
    }
    if not polynomial:
        return {}
    divisor = math.gcd(*(abs(value) for value in polynomial.values()))
    first = min(polynomial)
    sign = 1 if polynomial[first] > 0 else -1
    return {
        exponent: sign * coefficient // divisor
        for exponent, coefficient in polynomial.items()
    }


def packet_binomial(relation: Relation) -> Polynomial:
    polynomial: Polynomial = {}
    for sign, state in ((1, relation.left), (-1, relation.right)):
        exponent = state_exponent(state)
        polynomial[exponent] = (
            polynomial.get(exponent, 0) + sign * state_weight(state)
        )
    answer = normalize_polynomial(polynomial)
    if len(answer) != 2:
        raise AssertionError("a primitive packet must give a binomial")
    return answer


def derivative(
    polynomial: Polynomial,
    colours: tuple[int, ...],
    width: int,
) -> Polynomial:
    """Apply the sum of the selected normalized Taylor derivations."""

    answer: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        for colour in colours:
            offset = colour * width
            for level in range(width - 1):
                source = offset + level
                multiplicity = exponent[source]
                if not multiplicity:
                    continue
                target = source + 1
                derived = list(exponent)
                derived[source] -= 1
                derived[target] += 1
                key = tuple(derived)
                answer[key] = answer.get(key, 0) + (
                    coefficient * multiplicity * (level + 1)
                )
    return normalize_polynomial(answer)


def polynomial_key(polynomial: Polynomial) -> tuple[tuple[Exponent, int], ...]:
    return tuple(sorted(polynomial.items()))


def one_derivation_closure(
    polynomial: Polynomial,
    colours: tuple[int, ...],
    width: int,
) -> tuple[Polynomial, ...]:
    answer = []
    seen = set()
    current = polynomial
    while current:
        key = polynomial_key(current)
        if key not in seen:
            seen.add(key)
            answer.append(current)
        current = derivative(current, colours, width)
    return tuple(answer)


def independent_closure(
    polynomial: Polynomial,
    width: int,
) -> tuple[Polynomial, ...]:
    """Closure under the two commuting derivations, without duplicates."""

    answer = []
    seen = set()
    queue = deque([polynomial])
    while queue:
        current = queue.popleft()
        if not current:
            continue
        key = polynomial_key(current)
        if key in seen:
            continue
        seen.add(key)
        answer.append(current)
        queue.append(derivative(current, (0,), width))
        queue.append(derivative(current, (1,), width))
    return tuple(answer)


def derivative_closure(
    polynomial: Polynomial,
    mode: str,
    width: int,
) -> tuple[Polynomial, ...]:
    if mode == "operator":
        return one_derivation_closure(polynomial, (0,), width)
    if mode == "polynomial":
        return one_derivation_closure(polynomial, (1,), width)
    if mode == "diagonal":
        return one_derivation_closure(polynomial, (0, 1), width)
    if mode == "independent":
        return independent_closure(polynomial, width)
    raise ValueError(f"unknown mode: {mode}")


def polynomial_text(polynomial: Polynomial, variable_names: tuple[str, ...]) -> str:
    terms = []
    for exponent, coefficient in sorted(polynomial.items()):
        factors = []
        absolute = abs(coefficient)
        if absolute != 1 or not any(exponent):
            factors.append(str(absolute))
        for variable, power in zip(variable_names, exponent, strict=True):
            if power == 1:
                factors.append(variable)
            elif power > 1:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors) or "1"
        if not terms:
            terms.append(monomial if coefficient > 0 else f"-{monomial}")
        else:
            terms.append(("+" if coefficient > 0 else "-") + monomial)
    return "".join(terms) if terms else "0"


def singular_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.splitlines()[0].strip()


def singular_torus_test(
    closures: list[tuple[str, str, int, tuple[Polynomial, ...]]],
    executable: str,
) -> dict[tuple[str, str], dict[str, int | bool]]:
    """Test all closures in one exact Singular session."""

    lines = ["option(redSB);", "option(noredefine);"]
    for index, (digest, mode, width, polynomials) in enumerate(closures):
        variables = tuple(
            [f"R{level}" for level in range(width)]
            + [f"B{level}" for level in range(width)]
        )
        ring_name = f"q{index}"
        lines.append(
            f"ring {ring_name}=0,({','.join(variables)},T),dp;"
        )
        generators = [polynomial_text(value, variables) for value in polynomials]
        torus_product = "*".join(variables)
        generators.append(f"1-T*{torus_product}")
        lines.append(f"ideal I={','.join(generators)};")
        lines.append("ideal G=std(I);")
        marker = f"{digest}|{mode}"
        lines.append(f'print("BEGIN|{marker}");')
        lines.append("print(reduce(1,G));")
        lines.append("print(dim(G));")
        lines.append("print(size(G));")
        lines.append(f'print("END|{marker}");')

    completed = subprocess.run(
        [executable, "-q"],
        input="\n".join(lines) + "\n",
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode:
        raise RuntimeError(
            "Singular failed:\n" + completed.stdout + completed.stderr
        )

    results: dict[tuple[str, str], dict[str, int | bool]] = {}
    output = completed.stdout.splitlines()
    cursor = 0
    while cursor < len(output):
        line = output[cursor].strip()
        cursor += 1
        if not line.startswith("BEGIN|"):
            continue
        marker = line.removeprefix("BEGIN|")
        digest, mode = marker.split("|", 1)
        values = []
        while cursor < len(output) and not output[cursor].startswith("END|"):
            stripped = output[cursor].strip()
            cursor += 1
            if stripped:
                values.append(stripped)
        if cursor >= len(output) or output[cursor].strip() != f"END|{marker}":
            raise RuntimeError(f"truncated Singular block for {marker}")
        cursor += 1
        if len(values) != 3:
            raise RuntimeError(f"unexpected Singular block {marker}: {values}")
        normal_form, dimension, basis_size = map(int, values)
        if normal_form not in (0, 1):
            raise RuntimeError(
                f"unexpected normal form of one for {marker}: {normal_form}"
            )
        results[(digest, mode)] = {
            "torus_empty": normal_form == 0,
            "dimension": dimension,
            "standard_basis_size": basis_size,
        }

    if len(results) != len(closures):
        raise RuntimeError(
            f"parsed {len(results)} Singular blocks, expected {len(closures)}"
        )
    return results


def dense_product(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            answer[left_degree + right_degree] += left_value * right_value
    return answer


def dense_power(polynomial: list[int], exponent: int) -> list[int]:
    answer = [1]
    for _ in range(exponent):
        answer = dense_product(answer, polynomial)
    return answer


def three_state_linear_syzygy_control() -> dict[str, Any]:
    """Exact control showing why pairwise Graver separation is insufficient.

    At the torus point R_0=R_1=R_2=1 and

        (B_0(t),B_1(t),B_2(t))=(t^2+3t+2,2t+3,1),

    the Hall-weighted three-state row

        R_1^2 B_1^2 - 4 R_0 R_2 B_0 B_2 - R_0^2 B_2^2

    vanishes through its complete translation orbit.  The signs are supplied
    by the C4 label which marks R_0 and R_2 by one.  At scale two there are two
    logically distinct phase laws: fixed external signs stay negative, while
    genuine character phases are squared and become positive.  Both rows are
    nonzero, so this is a one-order linear-inheritance obstruction, not an
    all-scale packet.
    """

    b0 = [2, 3, 1]
    b1 = [3, 2]
    b2 = [1]

    def row(scale: int, phase_law: str) -> list[int]:
        middle_weight = math.comb(2 * scale, scale) ** 2
        negative_phase = -1 if phase_law == "fixed" else (-1) ** scale
        terms = (
            (1, dense_power(b1, 2 * scale)),
            (negative_phase * middle_weight, dense_power(b0, scale)),
            (negative_phase, dense_power(b2, 2 * scale)),
        )
        width = max(len(value) for _, value in terms)
        answer = [0] * width
        for coefficient, value in terms:
            for degree, entry in enumerate(value):
                answer[degree] += coefficient * entry
        while len(answer) > 1 and answer[-1] == 0:
            answer.pop()
        return answer

    first = row(1, "fixed")
    character_first = row(1, "character_power")
    fixed_second = row(2, "fixed")
    character_second = row(2, "character_power")
    if (
        first != [0]
        or character_first != first
        or fixed_second == [0]
        or character_second == [0]
    ):
        raise AssertionError(
            (first, character_first, fixed_second, character_second)
        )
    return {
        "span": 2,
        "states": [
            "R1^2*B1^2",
            "R0*R2*B0*B2",
            "R0^2*B2^2",
        ],
        "character_phases": [1, -1, -1],
        "scale_1_translation_coefficients": first,
        "scale_2_fixed_phase_translation_coefficients": fixed_second,
        "scale_2_fixed_phase_first_obstruction": next(
            [degree, value]
            for degree, value in enumerate(fixed_second)
            if value
        ),
        "scale_2_character_power_translation_coefficients": character_second,
        "scale_2_character_power_first_obstruction": next(
            [degree, value]
            for degree, value in enumerate(character_second)
            if value
        ),
    }


def selected_modes(text: str) -> tuple[str, ...]:
    modes = tuple(part for part in text.split(",") if part)
    if not modes or len(set(modes)) != len(modes):
        raise ValueError("--modes must be a nonempty list without duplicates")
    unknown = set(modes) - set(MODES)
    if unknown:
        raise ValueError(f"unknown modes: {sorted(unknown)}")
    return modes


def build_experiment(
    *,
    radial_degree: int,
    minimum_support: int,
    maximum_support: int | None,
    modes: tuple[str, ...],
    normaliz: str,
    singular: str,
) -> dict[str, Any]:
    raw_basis = colored_graver_basis(radial_degree, normaliz)
    normalized = {
        primitive_normalization(relation)
        for relation in raw_basis
        if relation_support_size(relation) >= minimum_support
    }
    if maximum_support is not None:
        normalized = {
            relation
            for relation in normalized
            if relation_support_size(relation) <= maximum_support
        }
    candidates = tuple(sorted(normalized))
    candidates = tuple(
        relation
        for relation in candidates
        if all(state_counts(relation.left))
    )
    transfers = Counter(
        abs(polynomial_level_transfer(relation)) for relation in candidates
    )
    if 0 in transfers:
        raise AssertionError(
            "a genuinely mixed Graver move had zero colour-level transfer"
        )
    factorial_survivors = tuple(
        relation
        for relation in candidates
        if exact_scalar_factorial_collision(relation)
    )

    relation_data: dict[str, dict[str, Any]] = {}
    closures = []
    for relation in factorial_survivors:
        digest = relation_digest(relation)
        width = len(relation.left.operator)
        binomial = packet_binomial(relation)
        mode_data = {}
        for mode in modes:
            closure = derivative_closure(binomial, mode, width)
            mode_data[mode] = {
                "derivative_generators": len(closure),
                "largest_generator_terms": max(map(len, closure)),
            }
            closures.append((digest, mode, width, closure))
        relation_data[digest] = {
            "digest": digest,
            "relation": relation_record(relation),
            "exact_marked_partition_collision": (
                exact_marked_partition_collision(relation)
            ),
            "polynomial_level_transfer": polynomial_level_transfer(relation),
            "endpoint_weight": state_weight(relation.left),
            "modes": mode_data,
        }

    singular_results = singular_torus_test(closures, singular)
    mode_counts = {mode: Counter() for mode in modes}
    survivor_digests = {mode: [] for mode in modes}
    for digest, record in relation_data.items():
        for mode in modes:
            result = singular_results[(digest, mode)]
            record["modes"][mode].update(result)
            outcome = "torus_empty" if result["torus_empty"] else "torus_survivor"
            mode_counts[mode][outcome] += 1
            if not result["torus_empty"]:
                survivor_digests[mode].append(digest)

    return {
        "model": (
            "Projected two-colour partition configuration; a proper orbit "
            "ideal is a promotion obstruction in this surrogate, not a "
            "rank-one Cartesian Hall packet or GVC(2) counterexample."
        ),
        "parameters": {
            "radial_degree": radial_degree,
            "minimum_support": minimum_support,
            "maximum_support": maximum_support,
            "modes": list(modes),
            "normaliz": normaliz_version(normaliz),
            "singular": singular_version(singular),
            "saturation": "all R_i and B_i coefficient channels",
        },
        "linear_syzygy_control": three_state_linear_syzygy_control(),
        "summary": {
            "raw_graver_basis": len(raw_basis),
            "normalized_unresolved_candidates": len(candidates),
            "factorial_partition_obstructions": (
                len(candidates) - len(factorial_survivors)
            ),
            "exact_scalar_factorial_survivors": len(factorial_survivors),
            "exact_marked_partition_survivors": sum(
                exact_marked_partition_collision(relation)
                for relation in factorial_survivors
            ),
            "absolute_colour_level_transfer_distribution": {
                str(transfer): count
                for transfer, count in sorted(transfers.items())
            },
            "translation_modes": {
                mode: {
                    **dict(sorted(mode_counts[mode].items())),
                    "survivor_digests": sorted(survivor_digests[mode]),
                }
                for mode in modes
            },
        },
        "relations": [relation_data[key] for key in sorted(relation_data)],
    }


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def logical_hash(value: dict[str, Any]) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radial-degree", type=int, default=4)
    parser.add_argument("--minimum-support", type=int, default=5)
    parser.add_argument("--maximum-support", type=int)
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--normaliz", default="normaliz")
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    if arguments.radial_degree < 1:
        parser.error("--radial-degree must be positive")
    if arguments.minimum_support < 2:
        parser.error("--minimum-support must be at least two")
    if (
        arguments.maximum_support is not None
        and arguments.maximum_support < arguments.minimum_support
    ):
        parser.error("--maximum-support must be at least --minimum-support")
    try:
        modes = selected_modes(arguments.modes)
    except ValueError as error:
        parser.error(str(error))

    normaliz = shutil.which(arguments.normaliz)
    singular = shutil.which(arguments.singular)
    if normaliz is None:
        parser.error(f"Normaliz executable not found: {arguments.normaliz}")
    if singular is None:
        parser.error(f"Singular executable not found: {arguments.singular}")

    result = build_experiment(
        radial_degree=arguments.radial_degree,
        minimum_support=arguments.minimum_support,
        maximum_support=arguments.maximum_support,
        modes=modes,
        normaliz=normaliz,
        singular=singular,
    )
    if arguments.output is not None:
        write_json(arguments.output, result)

    summary = result["summary"]
    print(f"Normaliz: {result['parameters']['normaliz']}")
    print(f"Singular: {result['parameters']['singular']}")
    print(
        "Graver census: "
        f"raw={summary['raw_graver_basis']}, "
        f"normalized unresolved={summary['normalized_unresolved_candidates']}, "
        "factorial obstructions="
        f"{summary['factorial_partition_obstructions']}, "
        "factorial survivors="
        f"{summary['exact_scalar_factorial_survivors']}"
    )
    for mode in modes:
        counts = summary["translation_modes"][mode]
        print(
            f"{mode}: torus empty={counts.get('torus_empty', 0)}, "
            f"survivors={counts.get('torus_survivor', 0)}"
        )
    print(f"logical_sha256={logical_hash(result)}")


if __name__ == "__main__":
    main()
