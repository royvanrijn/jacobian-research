"""Exact support-saturation compiler for finitely presented modules.

The public object in this module accepts a submodule ``N`` of a finite free
module, a boundary ideal ``I``, an optional distinguished class, and an
optional normal ideal.  It delegates the exact commutative algebra to
Singular and returns a JSON-serializable certificate for

    (N : I^infinity) / N = H^0_I(F/N).

Finite jets use ``N_n = N + m^n F`` for the declared normal ideal ``m``.
For consecutive requested orders the compiler also computes the image of
the transition from the saturated order-``n+1`` quotient to order ``n``.

Associated primes have two deliberately distinct modes:

``decompose``
    Compute a module primary decomposition with ``modDec`` and report every
    associated prime.

``regularity``
    Do not claim a list of associated primes.  Instead search the boundary
    ideal for a non-zero-divisor and certify, when found, that no associated
    prime contains the boundary ideal.

The second mode is useful for large presentations where full primary
decomposition is much more expensive than the support-saturation question.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
from dataclasses import asdict, dataclass
from typing import Any, Literal, Sequence


_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


@dataclass(frozen=True)
class PolynomialRing:
    """A polynomial ring understood by Singular."""

    variables: tuple[str, ...]
    characteristic: int = 0
    ordering: str = "dp"

    def __post_init__(self) -> None:
        if not self.variables:
            raise ValueError("the polynomial ring needs at least one variable")
        if len(set(self.variables)) != len(self.variables):
            raise ValueError("ring variables must be distinct")
        if any(not _NAME.fullmatch(name) for name in self.variables):
            raise ValueError("ring variables must be valid Singular names")
        if self.characteristic < 0:
            raise ValueError("the characteristic must be nonnegative")


@dataclass(frozen=True)
class ModulePresentation:
    """A submodule ``N`` of ``F=R^rank`` given by generating vectors."""

    ring: PolynomialRing
    rank: int
    generators: tuple[tuple[str, ...], ...]
    label: str = "module"
    singular_setup: str | None = None

    def __post_init__(self) -> None:
        if self.rank < 1:
            raise ValueError("the ambient free rank must be positive")
        if any(len(row) != self.rank for row in self.generators):
            raise ValueError("every module generator must have ambient rank")
        if not self.generators and self.singular_setup is None:
            raise ValueError(
                "provide module generators or a Singular setup program"
            )


@dataclass(frozen=True)
class NormalFiltration:
    """The normal ideal and finite jet orders to compile."""

    ideal: tuple[str, ...]
    orders: tuple[int, ...]
    strategy: Literal[
        "full_saturation",
        "distinguished_class_restriction",
        "distinguished_class_colon",
    ] = "full_saturation"
    transition_annihilator: str | None = None

    def __post_init__(self) -> None:
        if not self.ideal:
            raise ValueError("a normal filtration needs a nonzero ideal list")
        if not self.orders:
            raise ValueError("a normal filtration needs at least one order")
        if any(order < 1 for order in self.orders):
            raise ValueError("jet orders must be positive")
        if tuple(sorted(set(self.orders))) != self.orders:
            raise ValueError("jet orders must be strictly increasing")
        if (
            self.strategy == "distinguished_class_colon"
            and self.transition_annihilator is None
        ):
            raise ValueError(
                "the distinguished-class colon strategy needs an "
                "annihilating boundary element"
            )


@dataclass(frozen=True)
class CompilerOptions:
    """Backend and certificate policy."""

    associated_primes: Literal["decompose", "regularity"] = "decompose"
    saturation_strategy: Literal[
        "compute", "regularity", "perfect_height"
    ] = "compute"
    timeout_seconds: int = 900
    basis_algorithm: Literal["std", "slimgb"] = "std"
    regular_search_bound: int = 8
    torsion_exponent_bound: int = 32
    singular_binary: str | None = None

    def __post_init__(self) -> None:
        if self.timeout_seconds < 1:
            raise ValueError("timeout_seconds must be positive")
        if self.regular_search_bound < 0:
            raise ValueError("regular_search_bound must be nonnegative")
        if (
            self.saturation_strategy in {"regularity", "perfect_height"}
            and self.regular_search_bound == 0
        ):
            raise ValueError(
                "the regularity strategy needs a positive search bound"
            )
        if self.torsion_exponent_bound < 1:
            raise ValueError("torsion_exponent_bound must be positive")


class SupportSaturationError(RuntimeError):
    """Raised when the exact backend rejects or cannot finish a certificate."""


def _module(rows: Sequence[Sequence[str]]) -> str:
    if not rows:
        return "module(0)"
    return "module(" + ",".join(
        "[" + ",".join(row) + "]" for row in rows
    ) + ")"


def _ideal(items: Sequence[str]) -> str:
    return "ideal(" + (",".join(items) if items else "0") + ")"


def _section_printer(name: str, singular_object: str) -> str:
    return f"""
print("@@BEGIN:{name}");
for (int section_index=1; section_index<=size({singular_object}); section_index++)
{{
  print(string({singular_object}[section_index]));
}}
print("@@END:{name}");
"""


def _same_module(first: str, second: str) -> str:
    return (
        f"(size(simplify(reduce({first},std({second})),2))==0"
        f" && size(simplify(reduce({second},std({first})),2))==0)"
    )


def _candidate_expressions(
    boundary: Sequence[str], bound: int
) -> tuple[str, ...]:
    """Return deterministic candidate elements of the boundary ideal."""

    if bound == 0:
        return ()
    candidates: list[str] = []
    if len(boundary) > 1:
        candidates.append("+".join(f"({item})" for item in boundary))
        for seed in range(2, bound + 1):
            candidates.append(
                "+".join(
                    f"{seed ** index}*({item})"
                    for index, item in enumerate(boundary)
                )
            )
    candidates.extend(boundary)
    return tuple(dict.fromkeys(candidates))


def _parse_output(stdout: str) -> dict[str, Any]:
    scalars: dict[str, str] = {}
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if line.startswith("@@BEGIN:"):
            current = line.removeprefix("@@BEGIN:")
            sections[current] = []
        elif line.startswith("@@END:"):
            ended = line.removeprefix("@@END:")
            if current != ended:
                raise SupportSaturationError(
                    f"malformed Singular section: ended {ended!r}, "
                    f"expected {current!r}"
                )
            current = None
        elif current is not None:
            if line and line != "0" and not line.startswith("//"):
                sections[current].append(line)
        elif line.startswith("@@") and "=" in line:
            key, value = line[2:].split("=", 1)
            scalars[key] = value
    if current is not None:
        raise SupportSaturationError(f"unterminated Singular section {current}")
    return {"scalars": scalars, "sections": sections}


class SupportSaturationCompiler:
    """Compile exact support and finite-jet certificates with Singular."""

    schema = "support-saturation-certificate.v2"

    def __init__(self, options: CompilerOptions | None = None) -> None:
        self.options = options or CompilerOptions()

    def compile(
        self,
        presentation: ModulePresentation,
        boundary_ideal: Sequence[str],
        distinguished_class: Sequence[str] | None = None,
        filtration: NormalFiltration | None = None,
    ) -> dict[str, Any]:
        """Compile and return one JSON-serializable exact certificate."""

        boundary = tuple(boundary_ideal)
        if not boundary:
            raise ValueError("the boundary ideal needs at least one generator")
        if distinguished_class is not None:
            distinguished = tuple(distinguished_class)
            if len(distinguished) != presentation.rank:
                raise ValueError(
                    "the distinguished class must have ambient free rank"
                )
        else:
            distinguished = None
        if (
            filtration is not None
            and filtration.strategy
            in {
                "distinguished_class_restriction",
                "distinguished_class_colon",
            }
            and distinguished is None
        ):
            raise ValueError(
                "the distinguished-class jet strategy needs a class"
            )

        input_data = {
            "presentation": asdict(presentation),
            "boundary_ideal": list(boundary),
            "distinguished_class": (
                list(distinguished) if distinguished is not None else None
            ),
            "filtration": asdict(filtration) if filtration else None,
            "options": asdict(self.options),
        }
        canonical_input = json.dumps(
            input_data, sort_keys=True, separators=(",", ":")
        )
        input_hash = hashlib.sha256(canonical_input.encode()).hexdigest()
        program = self._program(
            presentation,
            boundary,
            distinguished,
            filtration,
        )
        singular = self.options.singular_binary or shutil.which("Singular")
        if singular is None:
            raise SupportSaturationError("Singular is required")
        version_result = subprocess.run(
            [singular, "--version"],
            text=True,
            capture_output=True,
            check=True,
        )
        singular_version = version_result.stdout.splitlines()[0].strip()
        completed = self._run(singular, program)
        parsed = _parse_output(completed.stdout)
        result = self._assemble(
            input_data,
            input_hash,
            parsed,
            completed.stderr,
            singular_version,
            filtration,
        )
        return result

    def compile_distinguished_jets(
        self,
        presentation: ModulePresentation,
        boundary_ideal: Sequence[str],
        distinguished_class: Sequence[str],
        filtration: NormalFiltration,
    ) -> dict[str, Any]:
        """Compile distinguished-class jets without claiming base saturation.

        This entry point is for presentations whose untruncated standard
        basis is unavailable but whose normal truncations are computable.  It
        deliberately omits ``(N:I^infinity)/N``, associated primes, and
        regularity claims for the base module.
        """

        boundary = tuple(boundary_ideal)
        distinguished = tuple(distinguished_class)
        if not boundary:
            raise ValueError("the boundary ideal needs at least one generator")
        if len(distinguished) != presentation.rank:
            raise ValueError(
                "the distinguished class must have ambient free rank"
            )
        if filtration.strategy not in {
            "full_saturation",
            "distinguished_class_restriction",
        }:
            raise ValueError(
                "finite-jet-only compilation supports full saturation or "
                "distinguished-class restriction"
            )
        input_data = {
            "presentation": asdict(presentation),
            "boundary_ideal": list(boundary),
            "distinguished_class": list(distinguished),
            "filtration": asdict(filtration),
            "options": asdict(self.options),
            "base_saturation": "not_computed",
        }
        canonical_input = json.dumps(
            input_data, sort_keys=True, separators=(",", ":")
        )
        input_hash = hashlib.sha256(canonical_input.encode()).hexdigest()
        program = self._distinguished_jet_program(
            presentation,
            boundary,
            distinguished,
            filtration,
        )
        singular = self.options.singular_binary or shutil.which("Singular")
        if singular is None:
            raise SupportSaturationError("Singular is required")
        version_result = subprocess.run(
            [singular, "--version"],
            text=True,
            capture_output=True,
            check=True,
        )
        singular_version = version_result.stdout.splitlines()[0].strip()
        completed = self._run(singular, program)
        parsed = _parse_output(completed.stdout)
        scalars: dict[str, str] = parsed["scalars"]

        def integer(name: str) -> int:
            try:
                return int(scalars[name])
            except (KeyError, ValueError) as error:
                raise SupportSaturationError(
                    f"missing integer marker {name}: {scalars}"
                ) from error

        def boolean(name: str) -> bool:
            value = integer(name)
            if value not in (0, 1):
                raise SupportSaturationError(
                    f"marker {name} is not Boolean: {value}"
                )
            return bool(value)

        if filtration.strategy == "full_saturation":
            jets = [
                {
                    "order": order,
                    "presentation_standard_basis_size": integer(
                        f"JET_{order}_PRESENTATION_SIZE"
                    ),
                    "saturated_standard_basis_size": integer(
                        f"JET_{order}_SATURATED_SIZE"
                    ),
                    "saturation_equal": boolean(
                        f"JET_{order}_SATURATION_EQUAL"
                    ),
                    "local_cohomology_generator_count": integer(
                        f"JET_{order}_H0_GENERATORS"
                    ),
                    "boundary_annihilation_exponent": integer(
                        f"JET_{order}_BOUNDARY_EXPONENT"
                    ),
                    "distinguished_class_zero": boolean(
                        f"JET_{order}_CLASS_ZERO"
                    ),
                    "distinguished_class_in_local_cohomology": boolean(
                        f"JET_{order}_CLASS_IN_H0"
                    ),
                }
                for order in filtration.orders
            ]
            transitions = [
                {
                    "from_order": upper,
                    "to_order": lower,
                    "image_generator_count": integer(
                        f"TRANSITION_{upper}_TO_{lower}_IMAGE_GENERATORS"
                    ),
                    "surjective": boolean(
                        f"TRANSITION_{upper}_TO_{lower}_SURJECTIVE"
                    ),
                    "distinguished_class_lifts": boolean(
                        f"TRANSITION_{upper}_TO_{lower}_CLASS_LIFTS"
                    ),
                }
                for lower, upper in zip(
                    filtration.orders, filtration.orders[1:]
                )
            ]
            scope = "full finite-jet local cohomology"
        else:
            jets = [
                {
                    "order": order,
                    "presentation_standard_basis_size": integer(
                        f"JET_{order}_PRESENTATION_SIZE"
                    ),
                    "distinguished_class_zero": boolean(
                        f"JET_{order}_CLASS_ZERO"
                    ),
                    "distinguished_class_in_local_cohomology": boolean(
                        f"JET_{order}_CLASS_IN_H0"
                    ),
                    "distinguished_class_boundary_exponent": integer(
                        f"JET_{order}_CLASS_BOUNDARY_EXPONENT"
                    ),
                    "boundary_power_tests": [
                        {
                            "exponent": exponent,
                            "annihilates_class": boolean(
                                f"JET_{order}_CLASS_KILLED_BY_POWER_{exponent}"
                            ),
                        }
                        for exponent in range(
                            1, self.options.torsion_exponent_bound + 1
                        )
                    ],
                }
                for order in filtration.orders
            ]
            transitions = [
                {
                    "from_order": upper,
                    "to_order": lower,
                    "distinguished_class_lifts": boolean(
                        f"TRANSITION_{upper}_TO_{lower}_CLASS_LIFTS"
                    ),
                }
                for lower, upper in zip(
                    filtration.orders, filtration.orders[1:]
                )
            ]
            scope = "distinguished class only"
        return {
            "schema": "support-saturation-finite-jet-certificate.v1",
            "input_sha256": input_hash,
            "input": input_data,
            "backend": {
                "name": "Singular",
                "version": singular_version,
                "stderr": completed.stderr.strip(),
            },
            "base_module": {
                "status": "not_computed",
                "claims": [],
            },
            "finite_jets": {
                "strategy": filtration.strategy,
                "normal_ideal": list(filtration.ideal),
                "convention": "N_n=N+m^n F",
                "scope": scope,
                "torsion_exponent_search_bound": (
                    self.options.torsion_exponent_bound
                ),
                "jets": jets,
                "transitions": transitions,
            },
        }

    def compile_distinguished_support_witness(
        self,
        presentation: ModulePresentation,
        boundary_ideal: Sequence[str],
        distinguished_class: Sequence[str],
    ) -> dict[str, Any]:
        """Certify a nonzero base support class without full saturation.

        ``InputPresentation`` supplied by a custom Singular setup must be a
        standard basis.  A nonzero class killed by a power of the boundary
        proves that base local cohomology is nonzero, that a
        boundary-containing associated prime exists, and that no element of
        the boundary ideal is regular.  It does not enumerate the full local
        cohomology module or its associated primes.
        """

        boundary = tuple(boundary_ideal)
        distinguished = tuple(distinguished_class)
        if not boundary:
            raise ValueError("the boundary ideal needs at least one generator")
        if len(distinguished) != presentation.rank:
            raise ValueError(
                "the distinguished class must have ambient free rank"
            )
        ring = presentation.ring
        if presentation.singular_setup is None:
            setup = "\n".join(
                [
                    (
                        f"ring support_ring={ring.characteristic},"
                        f"({','.join(ring.variables)}),{ring.ordering};"
                    ),
                    (
                        f"module InputPresentation="
                        f"{self.options.basis_algorithm}("
                        f"{_module(presentation.generators)});"
                    ),
                ]
            )
        else:
            setup = presentation.singular_setup
        lines = [
            setup,
            "module N=InputPresentation;",
            f"ideal Boundary=std({_ideal(boundary)});",
            f"vector Distinguished={_module((distinguished,))}[1];",
            'print("@@CLASS_ZERO="+string(reduce(Distinguished,N)==0));',
            "int ClassBoundaryExponent=-1;",
        ]
        for exponent in range(1, self.options.torsion_exponent_bound + 1):
            lines.extend(
                [
                    f"ideal BoundaryPower{exponent}=Boundary^{exponent};",
                    (
                        f"int ClassKilled{exponent}=size(simplify(reduce("
                        f"BoundaryPower{exponent}*module(Distinguished),"
                        "N),2))==0;"
                    ),
                    (
                        "if (ClassBoundaryExponent==-1"
                        f" && ClassKilled{exponent})"
                        "{"
                        f"ClassBoundaryExponent={exponent};"
                        "}"
                    ),
                    (
                        f'print("@@CLASS_KILLED_BY_POWER_{exponent}="+'
                        f"string(ClassKilled{exponent}));"
                    ),
                ]
            )
        lines.extend(
            [
                'print("@@CLASS_BOUNDARY_EXPONENT="+'
                "string(ClassBoundaryExponent));",
                'print("@@COMPLETE=1");',
                "quit;",
            ]
        )
        input_data = {
            "presentation": asdict(presentation),
            "boundary_ideal": list(boundary),
            "distinguished_class": list(distinguished),
            "options": asdict(self.options),
            "scope": "distinguished support witness",
        }
        canonical_input = json.dumps(
            input_data, sort_keys=True, separators=(",", ":")
        )
        input_hash = hashlib.sha256(canonical_input.encode()).hexdigest()
        singular = self.options.singular_binary or shutil.which("Singular")
        if singular is None:
            raise SupportSaturationError("Singular is required")
        version_result = subprocess.run(
            [singular, "--version"],
            text=True,
            capture_output=True,
            check=True,
        )
        singular_version = version_result.stdout.splitlines()[0].strip()
        completed = self._run(singular, "\n".join(lines))
        parsed = _parse_output(completed.stdout)
        scalars: dict[str, str] = parsed["scalars"]

        def integer(name: str) -> int:
            try:
                return int(scalars[name])
            except (KeyError, ValueError) as error:
                raise SupportSaturationError(
                    f"missing integer marker {name}: {scalars}"
                ) from error

        class_zero = bool(integer("CLASS_ZERO"))
        exponent = integer("CLASS_BOUNDARY_EXPONENT")
        tests = [
            {
                "exponent": power,
                "annihilates_class": bool(
                    integer(f"CLASS_KILLED_BY_POWER_{power}")
                ),
            }
            for power in range(1, self.options.torsion_exponent_bound + 1)
        ]
        nonzero_support_class = not class_zero and exponent > 0
        if not nonzero_support_class:
            raise SupportSaturationError(
                "the proposed distinguished support witness was not "
                "certified as nonzero boundary torsion"
            )
        return {
            "schema": "support-saturation-distinguished-witness.v1",
            "input_sha256": input_hash,
            "input": input_data,
            "backend": {
                "name": "Singular",
                "version": singular_version,
                "stderr": completed.stderr.strip(),
            },
            "local_cohomology": {
                "identity": "H^0_I(F/N) = (N:I^infinity)/N",
                "status": "certified_nonzero_by_distinguished_class",
                "full_module_computed": False,
            },
            "distinguished_class": {
                "zero_in_F_mod_N": class_zero,
                "belongs_to_local_cohomology": True,
                "boundary_annihilation_exponent": exponent,
                "boundary_power_tests": tests,
            },
            "associated_primes": {
                "status": (
                    "existence_of_boundary_containing_prime_certified; "
                    "primes_not_enumerated"
                ),
                "primes": None,
            },
            "regular_elements": {
                "candidate": None,
                "certifies_no_regular_element_in_boundary": True,
                "reason": (
                    "nonzero H^0_I implies an associated prime containing I"
                ),
            },
        }

    def _distinguished_jet_program(
        self,
        presentation: ModulePresentation,
        boundary: tuple[str, ...],
        distinguished: tuple[str, ...],
        filtration: NormalFiltration,
    ) -> str:
        """Return the Singular program for finite-jet-only compilation."""

        ring = presentation.ring
        basis = self.options.basis_algorithm
        lines: list[str] = (
            ['LIB "elim.lib";']
            if filtration.strategy == "full_saturation"
            else []
        )
        if presentation.singular_setup is None:
            lines.extend(
                [
                    (
                        f"ring support_ring={ring.characteristic},"
                        f"({','.join(ring.variables)}),{ring.ordering};"
                    ),
                    f"module N={_module(presentation.generators)};",
                ]
            )
        else:
            lines.extend(
                [
                    presentation.singular_setup,
                    "module N=InputPresentation;",
                ]
            )
        lines.extend(
            [
                f"ideal Boundary=std({_ideal(boundary)});",
                f"vector Distinguished={_module((distinguished,))}[1];",
                f"ideal Normal=std({_ideal(filtration.ideal)});",
                f"module AmbientFree=freemodule({presentation.rank});",
            ]
        )
        if filtration.strategy == "full_saturation":
            for order in filtration.orders:
                lines.extend(
                    [
                        (
                            f"module JetN{order}={basis}("
                            f"N+(Normal^{order})*AmbientFree);"
                        ),
                        (
                            f"module JetS{order}=std("
                            f"sat(JetN{order},Boundary));"
                        ),
                        (
                            f"module JetH{order}=simplify("
                            f"reduce(JetS{order},JetN{order}),2);"
                        ),
                        (
                            f'print("@@JET_{order}_PRESENTATION_SIZE="+'
                            f"string(size(JetN{order})));"
                        ),
                        (
                            f'print("@@JET_{order}_SATURATED_SIZE="+'
                            f"string(size(JetS{order})));"
                        ),
                        (
                            f'print("@@JET_{order}_H0_GENERATORS="+'
                            f"string(size(JetH{order})));"
                        ),
                        (
                            f'print("@@JET_{order}_SATURATION_EQUAL="+'
                            "string("
                            + _same_module(
                                f"JetN{order}", f"JetS{order}"
                            )
                            + "));"
                        ),
                        (
                            f'print("@@JET_{order}_CLASS_ZERO="+string('
                            f"reduce(Distinguished,JetN{order})==0));"
                        ),
                        (
                            f'print("@@JET_{order}_CLASS_IN_H0="+string('
                            f"reduce(Distinguished,JetS{order})==0));"
                        ),
                        f"int JetBoundaryExponent{order}=-1;",
                    ]
                )
                for exponent in range(
                    1, self.options.torsion_exponent_bound + 1
                ):
                    lines.extend(
                        [
                            (
                                f"ideal JetH0BoundaryPower{order}_{exponent}="
                                f"Boundary^{exponent};"
                            ),
                            (
                                f"int JetH0Killed{order}_{exponent}="
                                "size(simplify(reduce("
                                f"JetH0BoundaryPower{order}_{exponent}"
                                f"*JetS{order},JetN{order}),2))==0;"
                            ),
                            (
                                f"if (JetBoundaryExponent{order}==-1"
                                f" && JetH0Killed{order}_{exponent})"
                                "{"
                                f"JetBoundaryExponent{order}={exponent};"
                                "}"
                            ),
                        ]
                    )
                lines.append(
                    f'print("@@JET_{order}_BOUNDARY_EXPONENT="+string('
                    f"JetBoundaryExponent{order}));"
                )
            for lower, upper in zip(
                filtration.orders, filtration.orders[1:]
            ):
                lines.extend(
                    [
                        (
                            f"module JetImage{lower}_{upper}=std("
                            f"JetS{upper}+JetN{lower});"
                        ),
                        (
                            f"module JetImageDefect{lower}_{upper}="
                            "simplify(reduce("
                            f"JetImage{lower}_{upper},JetN{lower}),2);"
                        ),
                        (
                            f'print("@@TRANSITION_{upper}_TO_{lower}_'
                            "IMAGE_GENERATORS="
                            f'"+string(size(JetImageDefect{lower}_{upper})));'
                        ),
                        (
                            f'print("@@TRANSITION_{upper}_TO_{lower}_'
                            'SURJECTIVE="+string('
                            + _same_module(
                                f"JetImage{lower}_{upper}",
                                f"JetS{lower}",
                            )
                            + "));"
                        ),
                        (
                            f'print("@@TRANSITION_{upper}_TO_{lower}_'
                            'CLASS_LIFTS="+string(reduce('
                            "Distinguished,"
                            f"JetImage{lower}_{upper})==0));"
                        ),
                    ]
                )
            lines.extend(['print("@@COMPLETE=1");', "quit;"])
            return "\n".join(lines)
        for order in filtration.orders:
            lines.extend(
                [
                    (
                        f"module JetN{order}={basis}("
                        f"N+(Normal^{order})*AmbientFree);"
                    ),
                    (
                        f'print("@@JET_{order}_PRESENTATION_SIZE="+'
                        f"string(size(JetN{order})));"
                    ),
                    (
                        f'print("@@JET_{order}_CLASS_ZERO="+string('
                        f"reduce(Distinguished,JetN{order})==0));"
                    ),
                    f"int JetClassExponent{order}=-1;",
                ]
            )
            for exponent in range(
                1, self.options.torsion_exponent_bound + 1
            ):
                lines.extend(
                    [
                        (
                            f"ideal JetBoundaryPower{order}_{exponent}="
                            f"Boundary^{exponent};"
                        ),
                        (
                            f"int JetClassKilled{order}_{exponent}="
                            "size(simplify(reduce("
                            f"JetBoundaryPower{order}_{exponent}"
                            "*module(Distinguished),"
                            f"JetN{order}),2))==0;"
                        ),
                        (
                            f"if (JetClassExponent{order}==-1"
                            f" && JetClassKilled{order}_{exponent})"
                            "{"
                            f"JetClassExponent{order}={exponent};"
                            "}"
                        ),
                        (
                            f'print("@@JET_{order}_CLASS_KILLED_BY_POWER_'
                            f'{exponent}="+string('
                            f"JetClassKilled{order}_{exponent}));"
                        ),
                    ]
                )
            lines.extend(
                [
                    (
                        f'print("@@JET_{order}_CLASS_BOUNDARY_EXPONENT="+'
                        f"string(JetClassExponent{order}));"
                    ),
                    (
                        f'print("@@JET_{order}_CLASS_IN_H0="+string('
                        f"JetClassExponent{order}>0));"
                    ),
                ]
            )
        for lower, upper in zip(
            filtration.orders, filtration.orders[1:]
        ):
            lines.append(
                f'print("@@TRANSITION_{upper}_TO_{lower}_'
                "CLASS_LIFTS=\"+string("
                f"JetClassExponent{upper}>0"
                f" && JetClassExponent{lower}>0));"
            )
        lines.extend(['print("@@COMPLETE=1");', "quit;"])
        return "\n".join(lines)

    def _run(
        self, singular: str, program: str
    ) -> subprocess.CompletedProcess[str]:
        process = subprocess.Popen(
            [singular, "-q"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(
                program, timeout=self.options.timeout_seconds
            )
        except subprocess.TimeoutExpired as error:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise SupportSaturationError(
                "Singular timed out after "
                f"{self.options.timeout_seconds} seconds"
            ) from error
        completed = subprocess.CompletedProcess(
            [singular, "-q"], process.returncode, stdout, stderr
        )
        if process.returncode:
            raise SupportSaturationError(stdout + stderr)
        if "? error occurred" in stdout or "? error occurred" in stderr:
            raise SupportSaturationError(stdout + stderr)
        if "@@COMPLETE=1" not in stdout:
            raise SupportSaturationError(
                "Singular did not emit the completion marker\n"
                + stdout
                + stderr
            )
        return completed

    def _program(
        self,
        presentation: ModulePresentation,
        boundary: tuple[str, ...],
        distinguished: tuple[str, ...] | None,
        filtration: NormalFiltration | None,
    ) -> str:
        ring = presentation.ring
        basis = self.options.basis_algorithm
        candidates = _candidate_expressions(
            boundary, self.options.regular_search_bound
        )
        if (
            self.options.saturation_strategy == "perfect_height"
            and presentation.rank != 1
        ):
            raise ValueError(
                "the perfect-height strategy is restricted to rank one"
            )
        lines = [
            'LIB "primdec.lib";',
            'LIB "mprimdec.lib";',
        ]
        if presentation.singular_setup is None:
            lines.extend(
                [
                    (
                        f"ring support_ring={ring.characteristic},"
                        f"({','.join(ring.variables)}),{ring.ordering};"
                    ),
                    (
                        f"module N={basis}("
                        f"{_module(presentation.generators)});"
                    ),
                ]
            )
        else:
            lines.extend(
                [
                    presentation.singular_setup,
                    f"module N={basis}(InputPresentation);",
                ]
            )
        lines.extend(
            [
            f"ideal Boundary=std({_ideal(boundary)});",
            "int FoundRegular=0;",
            ]
        )
        if self.options.saturation_strategy == "perfect_height":
            lines.extend(
                [
                    "ideal PresentationIdeal=std(ideal(N));",
                    "resolution PresentationResolution=mres("
                    "PresentationIdeal,0);",
                    "int PresentationDimension=dim(PresentationIdeal);",
                    (
                        "int PerfectCodimensionTwo=("
                        "PresentationDimension==nvars(basering)-2"
                        " && size(PresentationResolution)==3);"
                    ),
                    'print("@@PERFECT_CODIMENSION_TWO="+'
                    "string(PerfectCodimensionTwo));",
                ]
            )
        for index, candidate in enumerate(candidates, start=1):
            candidate_test = (
                f"""
if (!FoundRegular)
{{
  WasTested{index}=1;
  module RegularColon{index}={basis}(
    quotient(N,ideal(RegularCandidate{index}))
  );
  IsRegular{index}={_same_module(f"RegularColon{index}", "N")};
  if (IsRegular{index})
  {{
    FoundRegular=1;
  }}
}}
"""
                if self.options.saturation_strategy != "perfect_height"
                else f"""
if (!FoundRegular)
{{
  WasTested{index}=1;
  ideal CandidateCut{index}={basis}(
    PresentationIdeal+ideal(RegularCandidate{index})
  );
  IsRegular{index}=(
    PerfectCodimensionTwo
    && dim(CandidateCut{index})==PresentationDimension-1
  );
  if (IsRegular{index})
  {{
    FoundRegular=1;
  }}
}}
"""
            )
            lines.extend(
                [
                    f"poly RegularCandidate{index}={candidate};",
                    f"int IsRegular{index}=0;",
                    f"int WasTested{index}=0;",
                    candidate_test,
                    (
                        f'print("@@REGULAR_{index}="+'
                        f"string(IsRegular{index}));"
                    ),
                    (
                        f'print("@@REGULAR_TESTED_{index}="+'
                        f"string(WasTested{index}));"
                    ),
                    (
                        f'print("@@REGULAR_EXPR_{index}="+'
                        f"string(RegularCandidate{index}));"
                    ),
                ]
            )
        lines.append(f'print("@@REGULAR_COUNT={len(candidates)}");')
        if self.options.saturation_strategy == "compute":
            lines.extend(
                [
                    "module Saturated=std(sat(N,Boundary));",
                    (
                        "module H0Generators=simplify("
                        "reduce(Saturated,N),2);"
                    ),
                    'print("@@SATURATION_STRATEGY=compute");',
                ]
            )
        else:
            lines.extend(
                [
                    """
if (!FoundRegular)
{
  ERROR("regularity strategy found no regular boundary element");
}
""",
                    "module Saturated=N;",
                    "module H0Generators=0;",
                    (
                        'print("@@SATURATION_STRATEGY='
                        f'{self.options.saturation_strategy}");'
                    ),
                ]
            )
        lines.extend(
            [
            'print("@@AMBIENT_RANK=' + str(presentation.rank) + '");',
            'print("@@PRESENTATION_BASIS_SIZE="+string(size(N)));',
            'print("@@SATURATED_BASIS_SIZE="+string(size(Saturated)));',
            'print("@@H0_GENERATOR_COUNT="+string(size(H0Generators)));',
            (
                "int SaturationEqual="
                + _same_module("N", "Saturated")
                + ";"
            ),
            'print("@@SATURATION_EQUAL="+string(SaturationEqual));',
            _section_printer("SATURATED_PRESENTATION", "Saturated"),
            _section_printer("H0_GENERATORS", "H0Generators"),
            """
int H0BoundaryExponent=0;
if (size(H0Generators)>0)
{
  H0BoundaryExponent=-1;
  ideal H0BoundaryPower=Boundary;
""",
            (
                "for (int h0_exponent_index=1; h0_exponent_index<="
                f"{self.options.torsion_exponent_bound}; "
                "h0_exponent_index++)"
            ),
            """
  {
    if (H0BoundaryExponent==-1)
    {
      if (size(simplify(
          reduce(H0BoundaryPower*Saturated,N),2
        ))==0)
      {
        H0BoundaryExponent=h0_exponent_index;
      }
      else
      {
        H0BoundaryPower=H0BoundaryPower*Boundary;
      }
    }
  }
}
print("@@H0_BOUNDARY_EXPONENT="+string(H0BoundaryExponent));
""",
            (
                "if (size(H0Generators)>0) "
                "{ module H0Relations=std(modulo(H0Generators,N)); "
                'ideal H0Annihilator=std(annil(H0Relations)); '
                + _section_printer("H0_RELATIONS", "H0Relations")
                + _section_printer("H0_ANNIHILATOR", "H0Annihilator")
                + "}"
            ),
            ]
        )

        if self.options.associated_primes == "decompose":
            decomposition_command = (
                "primdecSY(std(ideal(N)))"
                if presentation.rank == 1
                else "modDec(N)"
            )
            lines.extend(
                [
                    f"list AmbientDecomposition={decomposition_command};",
                    (
                        'print("@@AMBIENT_ASSOCIATED_COUNT="+'
                        "string(size(AmbientDecomposition)));"
                    ),
                    """
for (int associated_index=1;
     associated_index<=size(AmbientDecomposition);
     associated_index++)
{
  ideal AssociatedPrime=std(AmbientDecomposition[associated_index][2]);
  print(
    "@@ASSOCIATED_CONTAINS_BOUNDARY_"+string(associated_index)
    +"="+string(size(simplify(reduce(Boundary,AssociatedPrime),2))==0)
  );
  print("@@BEGIN:ASSOCIATED_PRIME_"+string(associated_index));
  for (int prime_generator=1;
       prime_generator<=size(AssociatedPrime);
       prime_generator++)
  {
    print(string(AssociatedPrime[prime_generator]));
  }
  print("@@END:ASSOCIATED_PRIME_"+string(associated_index));
}
""",
                    """
int H0AssociatedCount=0;
if (size(H0Generators)>0)
{
  H0AssociatedCount=-2;
}
print("@@H0_ASSOCIATED_COUNT="+string(H0AssociatedCount));
""",
                ]
            )
        else:
            lines.extend(
                [
                    'print("@@AMBIENT_ASSOCIATED_COUNT=-1");',
                    'print("@@H0_ASSOCIATED_COUNT=-1");',
                ]
            )

        if distinguished is not None:
            lines.extend(
                [
                    f"vector Distinguished={_module((distinguished,))}[1];",
                    "vector DistinguishedRemainder=reduce(Distinguished,N);",
                    (
                        'print("@@CLASS_ZERO="+string('
                        "DistinguishedRemainder==0));"
                    ),
                    (
                        'print("@@CLASS_IN_H0="+string('
                        "reduce(Distinguished,Saturated)==0));"
                    ),
                    (
                        "module DistinguishedKernel=std("
                        "modulo(module(Distinguished),N));"
                    ),
                    (
                        "ideal DistinguishedAnnihilator=std("
                        "ideal(DistinguishedKernel));"
                    ),
                    _section_printer(
                        "DISTINGUISHED_REMAINDER",
                        "module(DistinguishedRemainder)",
                    ),
                    _section_printer(
                        "DISTINGUISHED_ANNIHILATOR",
                        "DistinguishedAnnihilator",
                    ),
                    "int DistinguishedExponent=-1;",
                    "ideal BoundaryPower=Boundary;",
                    "if (!SaturationEqual)",
                    "{",
                    (
                        "for (int exponent_index=1; exponent_index<="
                        f"{self.options.torsion_exponent_bound}; "
                        "exponent_index++)"
                    ),
                    """
{
  if (DistinguishedExponent==-1)
  {
    if (size(simplify(
        reduce(BoundaryPower*module(Distinguished),N),2
      ))==0)
    {
      DistinguishedExponent=exponent_index;
    }
    else
    {
      BoundaryPower=BoundaryPower*Boundary;
    }
  }
}
""",
                    "}",
                    (
                        'print("@@CLASS_BOUNDARY_EXPONENT="+'
                        "string(DistinguishedExponent));"
                    ),
                ]
            )

        if filtration is not None:
            lines.extend(
                [
                    f"ideal Normal=std({_ideal(filtration.ideal)});",
                    f"module AmbientFree=freemodule({presentation.rank});",
                    f'print("@@JET_COUNT={len(filtration.orders)}");',
                    f'print("@@JET_STRATEGY={filtration.strategy}");',
                ]
            )
            if filtration.strategy == "full_saturation":
                for order in filtration.orders:
                    lines.extend(
                        [
                            (
                                f"module JetN{order}={basis}("
                                f"N+(Normal^{order})*AmbientFree);"
                            ),
                            (
                                f"module JetS{order}=std("
                                f"sat(JetN{order},Boundary));"
                            ),
                            (
                                f"module JetH{order}=simplify("
                                f"reduce(JetS{order},JetN{order}),2);"
                            ),
                            (
                                f'print("@@JET_{order}_PRESENTATION_SIZE="+'
                                f"string(size(JetN{order})));"
                            ),
                            (
                                f'print("@@JET_{order}_SATURATED_SIZE="+'
                                f"string(size(JetS{order})));"
                            ),
                            (
                                f'print("@@JET_{order}_H0_GENERATORS="+'
                                f"string(size(JetH{order})));"
                            ),
                            (
                                f'print("@@JET_{order}_SATURATION_EQUAL="+'
                                "string("
                                + _same_module(
                                    f"JetN{order}", f"JetS{order}"
                                )
                                + "));"
                            ),
                            _section_printer(
                                f"JET_{order}_H0", f"JetH{order}"
                            ),
                            f"int JetBoundaryExponent{order}=0;",
                            f"""
if (size(JetH{order})>0)
{{
  JetBoundaryExponent{order}=-1;
  ideal JetBoundaryPower{order}=Boundary;
""",
                            (
                                "for (int jet_exponent_index=1; "
                                "jet_exponent_index<="
                                f"{self.options.torsion_exponent_bound}; "
                                "jet_exponent_index++)"
                            ),
                            f"""
  {{
    if (JetBoundaryExponent{order}==-1)
    {{
      if (size(simplify(
          reduce(JetBoundaryPower{order}*JetS{order},JetN{order}),2
        ))==0)
      {{
        JetBoundaryExponent{order}=jet_exponent_index;
      }}
      else
      {{
        JetBoundaryPower{order}=JetBoundaryPower{order}*Boundary;
      }}
    }}
  }}
}}
print("@@JET_{order}_BOUNDARY_EXPONENT="
  +string(JetBoundaryExponent{order}));
""",
                        ]
                    )
                    if distinguished is not None:
                        lines.extend(
                            [
                                (
                                    f'print("@@JET_{order}_CLASS_ZERO="+'
                                    "string(reduce("
                                    f"Distinguished,JetN{order})==0));"
                                ),
                                (
                                    f'print("@@JET_{order}_CLASS_IN_H0="+'
                                    "string(reduce("
                                    f"Distinguished,JetS{order})==0));"
                                ),
                            ]
                        )
                for lower, upper in zip(
                    filtration.orders, filtration.orders[1:]
                ):
                    lines.extend(
                        [
                            (
                                f"module JetImage{lower}_{upper}=std("
                                f"JetS{upper}+JetN{lower});"
                            ),
                            (
                                f"module JetImageDefect{lower}_{upper}="
                                "simplify(reduce("
                                f"JetImage{lower}_{upper},JetN{lower}),2);"
                            ),
                            (
                                f'print("@@TRANSITION_{upper}_TO_{lower}_'
                                "IMAGE_GENERATORS="
                                f'"+string(size(JetImageDefect{lower}_{upper})));'
                            ),
                            (
                                f'print("@@TRANSITION_{upper}_TO_{lower}_'
                                'SURJECTIVE="+string('
                                + _same_module(
                                    f"JetImage{lower}_{upper}",
                                    f"JetS{lower}",
                                )
                                + "));"
                            ),
                            _section_printer(
                                f"TRANSITION_{upper}_TO_{lower}_IMAGE",
                                f"JetImageDefect{lower}_{upper}",
                            ),
                        ]
                    )
            else:
                annihilator = filtration.transition_annihilator
                if annihilator is not None:
                    lines.append(
                        f"poly JetTransitionAnnihilator={annihilator};"
                    )
                for order in filtration.orders:
                    lines.extend(
                        [
                            (
                                f"module JetN{order}={basis}("
                                f"N+(Normal^{order})*AmbientFree);"
                            ),
                            (
                                f'print("@@JET_{order}_PRESENTATION_SIZE="+'
                                f"string(size(JetN{order})));"
                            ),
                            (
                                f'print("@@JET_{order}_CLASS_ZERO="+'
                                "string(reduce("
                                f"Distinguished,JetN{order})==0));"
                            ),
                            f"int JetClassExponent{order}=-1;",
                            f"ideal JetClassBoundaryPower{order}=Boundary;",
                            (
                                "for (int jet_class_exponent=1; "
                                "jet_class_exponent<="
                                f"{self.options.torsion_exponent_bound}; "
                                "jet_class_exponent++)"
                            ),
                            f"""
{{
  if (JetClassExponent{order}==-1)
  {{
    if (size(simplify(reduce(
        JetClassBoundaryPower{order}*module(Distinguished),
        JetN{order}),2))==0)
    {{
      JetClassExponent{order}=jet_class_exponent;
    }}
    else
    {{
      JetClassBoundaryPower{order}=
        JetClassBoundaryPower{order}*Boundary;
    }}
  }}
}}
print("@@JET_{order}_CLASS_BOUNDARY_EXPONENT="
  +string(JetClassExponent{order}));
print("@@JET_{order}_CLASS_IN_H0="
  +string(JetClassExponent{order}>0));
""",
                        ]
                    )
                for lower, upper in zip(
                    filtration.orders, filtration.orders[1:]
                ):
                    if filtration.strategy == "distinguished_class_colon":
                        lines.extend(
                            [
                                (
                                    f"module JetColon{upper}={basis}(quotient("
                                    f"JetN{upper},"
                                    "ideal(JetTransitionAnnihilator)));"
                                ),
                                (
                                    "module "
                                    f"JetClassTransition{lower}_{upper}="
                                    f"{basis}(JetColon{upper}+JetN{lower});"
                                ),
                                (
                                    f'print("@@TRANSITION_{upper}_TO_{lower}_'
                                    'CLASS_LIFTS="+string(reduce('
                                    "Distinguished,"
                                    f"JetClassTransition{lower}_{upper})"
                                    "==0));"
                                ),
                            ]
                        )
                    else:
                        lines.append(
                            f'print("@@TRANSITION_{upper}_TO_{lower}_'
                            "CLASS_LIFTS=\"+string("
                            f"JetClassExponent{upper}>0"
                            f" && JetClassExponent{lower}>0));"
                        )

        lines.extend(['print("@@COMPLETE=1");', "quit;"])
        return "\n".join(lines)

    def _assemble(
        self,
        input_data: dict[str, Any],
        input_hash: str,
        parsed: dict[str, Any],
        stderr: str,
        singular_version: str,
        filtration: NormalFiltration | None,
    ) -> dict[str, Any]:
        scalars: dict[str, str] = parsed["scalars"]
        sections: dict[str, list[str]] = parsed["sections"]

        def integer(name: str) -> int:
            try:
                return int(scalars[name])
            except (KeyError, ValueError) as error:
                raise SupportSaturationError(
                    f"missing integer marker {name}: {scalars}"
                ) from error

        def boolean(name: str) -> bool:
            value = integer(name)
            if value not in (0, 1):
                raise SupportSaturationError(
                    f"marker {name} is not Boolean: {value}"
                )
            return bool(value)

        regular_tests = [
            {
                "element": scalars[f"REGULAR_EXPR_{index}"],
                "tested": boolean(f"REGULAR_TESTED_{index}"),
                "regular": boolean(f"REGULAR_{index}"),
            }
            for index in range(1, integer("REGULAR_COUNT") + 1)
        ]
        regular = next(
            (
                item["element"]
                for item in regular_tests
                if item["regular"]
            ),
            None,
        )
        associated_count = integer("AMBIENT_ASSOCIATED_COUNT")
        if associated_count >= 0:
            associated = [
                sections[f"ASSOCIATED_PRIME_{index}"]
                for index in range(1, associated_count + 1)
            ]
            boundary_containing = [
                associated[index - 1]
                for index in range(1, associated_count + 1)
                if boolean(f"ASSOCIATED_CONTAINS_BOUNDARY_{index}")
            ]
            associated_status = "computed_by_module_primary_decomposition"
        else:
            associated = None
            boundary_containing = None
            associated_status = (
                "boundary_noncontainment_certified_by_regular_element"
                if regular is not None
                else "not_computed_and_no_regular_candidate_found"
            )
        h0_associated_count = integer("H0_ASSOCIATED_COUNT")
        if h0_associated_count == -2:
            h0_associated = boundary_containing
        elif h0_associated_count >= 0:
            h0_associated = [
                sections[f"H0_ASSOCIATED_PRIME_{index}"]
                for index in range(1, h0_associated_count + 1)
            ]
        else:
            h0_associated = None
        saturation_equal = boolean("SATURATION_EQUAL")
        if boundary_containing is not None:
            if saturation_equal and boundary_containing:
                raise SupportSaturationError(
                    "primary decomposition contradicts zero boundary "
                    "local cohomology"
                )
            if not saturation_equal and not boundary_containing:
                raise SupportSaturationError(
                    "primary decomposition missed the nonzero boundary "
                    "local cohomology support"
                )

        result: dict[str, Any] = {
            "schema": self.schema,
            "input_sha256": input_hash,
            "input": input_data,
            "backend": {
                "name": "Singular",
                "version": singular_version,
                "stderr": stderr.strip(),
            },
            "saturation": {
                "strategy": scalars["SATURATION_STRATEGY"],
                "equal_to_presentation": saturation_equal,
                "presentation_standard_basis_size": integer(
                    "PRESENTATION_BASIS_SIZE"
                ),
                "saturated_standard_basis_size": integer(
                    "SATURATED_BASIS_SIZE"
                ),
                "saturated_presentation": sections[
                    "SATURATED_PRESENTATION"
                ],
            },
            "local_cohomology": {
                "identity": "H^0_I(F/N) = (N:I^infinity)/N",
                "zero": saturation_equal,
                "generator_count": integer("H0_GENERATOR_COUNT"),
                "boundary_annihilation_exponent": integer(
                    "H0_BOUNDARY_EXPONENT"
                ),
                "torsion_exponent_search_bound": (
                    self.options.torsion_exponent_bound
                ),
                "generators_in_F": sections["H0_GENERATORS"],
                "relations": sections.get("H0_RELATIONS", []),
                "annihilator": sections.get("H0_ANNIHILATOR", []),
                "associated_primes": h0_associated,
            },
            "associated_primes": {
                "status": associated_status,
                "primes": associated,
                "boundary_containing_primes": boundary_containing,
            },
            "regular_elements": {
                "candidate": regular,
                "tests": regular_tests,
                "certifies_boundary_grade_at_least_one": regular is not None,
                "certifies_no_regular_element_in_boundary": bool(
                    boundary_containing
                ),
            },
        }
        if input_data["distinguished_class"] is not None:
            result["distinguished_class"] = {
                "zero_in_F_mod_N": boolean("CLASS_ZERO"),
                "belongs_to_local_cohomology": boolean("CLASS_IN_H0"),
                "boundary_annihilation_exponent": integer(
                    "CLASS_BOUNDARY_EXPONENT"
                ),
                "remainder": sections["DISTINGUISHED_REMAINDER"],
                "annihilator": sections["DISTINGUISHED_ANNIHILATOR"],
            }
        if filtration is not None:
            jets = []
            transitions = []
            if filtration.strategy == "full_saturation":
                for order in filtration.orders:
                    jet: dict[str, Any] = {
                        "order": order,
                        "presentation_standard_basis_size": integer(
                            f"JET_{order}_PRESENTATION_SIZE"
                        ),
                        "saturated_standard_basis_size": integer(
                            f"JET_{order}_SATURATED_SIZE"
                        ),
                        "saturation_equal": boolean(
                            f"JET_{order}_SATURATION_EQUAL"
                        ),
                        "local_cohomology_generator_count": integer(
                            f"JET_{order}_H0_GENERATORS"
                        ),
                        "boundary_annihilation_exponent": integer(
                            f"JET_{order}_BOUNDARY_EXPONENT"
                        ),
                        "local_cohomology_generators": sections[
                            f"JET_{order}_H0"
                        ],
                    }
                    if input_data["distinguished_class"] is not None:
                        jet["distinguished_class_zero"] = boolean(
                            f"JET_{order}_CLASS_ZERO"
                        )
                        jet[
                            "distinguished_class_in_local_cohomology"
                        ] = boolean(f"JET_{order}_CLASS_IN_H0")
                    jets.append(jet)
                for lower, upper in zip(
                    filtration.orders, filtration.orders[1:]
                ):
                    transitions.append(
                        {
                            "from_order": upper,
                            "to_order": lower,
                            "image_generator_count": integer(
                                f"TRANSITION_{upper}_TO_{lower}_"
                                "IMAGE_GENERATORS"
                            ),
                            "surjective": boolean(
                                f"TRANSITION_{upper}_TO_{lower}_SURJECTIVE"
                            ),
                            "image_generators": sections[
                                f"TRANSITION_{upper}_TO_{lower}_IMAGE"
                            ],
                        }
                    )
                exponents = [
                    jet["boundary_annihilation_exponent"] for jet in jets
                ]
                jet_summary: dict[str, Any] = {
                    "uniform_boundary_exponent_on_requested_jets": (
                        all(exponent >= 0 for exponent in exponents)
                    ),
                    "maximum_requested_jet_boundary_exponent": (
                        max(exponents)
                        if all(exponent >= 0 for exponent in exponents)
                        else None
                    ),
                }
            else:
                for order in filtration.orders:
                    jets.append(
                        {
                            "order": order,
                            "presentation_standard_basis_size": integer(
                                f"JET_{order}_PRESENTATION_SIZE"
                            ),
                            "distinguished_class_zero": boolean(
                                f"JET_{order}_CLASS_ZERO"
                            ),
                            "distinguished_class_in_local_cohomology": (
                                boolean(f"JET_{order}_CLASS_IN_H0")
                            ),
                            "distinguished_class_boundary_exponent": integer(
                                f"JET_{order}_CLASS_BOUNDARY_EXPONENT"
                            ),
                        }
                    )
                for lower, upper in zip(
                    filtration.orders, filtration.orders[1:]
                ):
                    transitions.append(
                        {
                            "from_order": upper,
                            "to_order": lower,
                            "distinguished_class_lifts": boolean(
                                f"TRANSITION_{upper}_TO_{lower}_CLASS_LIFTS"
                            ),
                            "annihilating_element": (
                                filtration.transition_annihilator
                            ),
                        }
                    )
                jet_summary = {
                    "scope": (
                        "distinguished class only; full finite-jet local "
                        "cohomology was not computed"
                    )
                }
            result["finite_jets"] = {
                "strategy": filtration.strategy,
                "normal_ideal": list(filtration.ideal),
                "convention": "N_n=N+m^n F",
                "torsion_exponent_search_bound": (
                    self.options.torsion_exponent_bound
                ),
                **jet_summary,
                "jets": jets,
                "transitions": transitions,
            }
        return result


def certificate_json(certificate: dict[str, Any]) -> str:
    """Serialize a certificate deterministically."""

    return json.dumps(certificate, indent=2, sort_keys=True) + "\n"
