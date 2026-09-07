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
from typing import Any, Literal, Mapping, Sequence


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


@dataclass(frozen=True)
class CertificateAssurance:
    """Interpretation of an exact backend computation.

    Singular performs exact arithmetic in both characteristic zero and a
    declared prime characteristic.  ``claim`` records whether that result is
    intended as an exact certificate for the target coefficient field or as
    modular evidence for a characteristic-zero question.
    """

    claim: Literal["auto", "exact", "modular"] = "auto"
    target_characteristic: int | None = None
    note: str | None = None

    def __post_init__(self) -> None:
        if self.claim not in {"auto", "exact", "modular"}:
            raise ValueError(f"unsupported assurance claim {self.claim!r}")
        if (
            self.target_characteristic is not None
            and self.target_characteristic < 0
        ):
            raise ValueError("target_characteristic must be nonnegative")


@dataclass(frozen=True)
class SupportSaturationProblem:
    """Shared input schema for support-local-cohomology calculations.

    ``support_ideal`` is the ideal ``I`` in ``H^0_I(M)``.
    ``completion_ideal`` is the ideal used for the finite tower
    ``M/completion_ideal^n M``.  Variable roles are metadata, but are
    validated against the declared ring so adapters cannot silently exchange
    parameter and normal directions.
    """

    presentation: ModulePresentation
    support_ideal: tuple[str, ...]
    completion_ideal: tuple[str, ...] = ()
    parameter_base_variables: tuple[str, ...] = ()
    normal_variables: tuple[str, ...] = ()
    jet_orders: tuple[int, ...] = ()
    distinguished_class: tuple[str, ...] | None = None
    jet_strategy: Literal[
        "full_saturation",
        "distinguished_class_restriction",
        "distinguished_class_colon",
    ] = "full_saturation"
    transition_annihilator: str | None = None
    assurance: CertificateAssurance = CertificateAssurance()

    schema = "support-saturation-input.v1"

    def __post_init__(self) -> None:
        if not self.support_ideal:
            raise ValueError("support_ideal needs at least one generator")
        variables = set(self.presentation.ring.variables)
        parameters = set(self.parameter_base_variables)
        normals = set(self.normal_variables)
        unknown = (parameters | normals) - variables
        if unknown:
            raise ValueError(
                "variable roles contain names outside the ring: "
                + ", ".join(sorted(unknown))
            )
        overlap = parameters & normals
        if overlap:
            raise ValueError(
                "parameter/base and normal variables must be disjoint: "
                + ", ".join(sorted(overlap))
            )
        if self.jet_orders and not self.completion_ideal:
            raise ValueError(
                "jet_orders require a nonempty completion_ideal"
            )
        if self.completion_ideal and not self.jet_orders:
            raise ValueError(
                "completion_ideal requires at least one requested jet order"
            )
        if tuple(sorted(set(self.jet_orders))) != self.jet_orders:
            raise ValueError("jet_orders must be strictly increasing")
        if any(order < 1 for order in self.jet_orders):
            raise ValueError("jet_orders must be positive")
        if (
            self.distinguished_class is not None
            and len(self.distinguished_class) != self.presentation.rank
        ):
            raise ValueError(
                "the distinguished class must have ambient free rank"
            )

    @classmethod
    def from_mapping(
        cls, data: Mapping[str, Any]
    ) -> "SupportSaturationProblem":
        """Parse the public JSON-compatible input schema."""

        schema = data.get("schema", cls.schema)
        if schema != cls.schema:
            raise ValueError(
                f"unsupported support-saturation input schema {schema!r}"
            )
        try:
            ring_data = data["ring"]
            module_data = data["module_presentation"]
            support_ideal = data["support_ideal"]
            completion_ideal = data["completion_ideal"]
            parameter_variables = data["parameter_base_variables"]
            normal_variables = data["normal_variables"]
        except KeyError as error:
            raise ValueError(
                f"missing support-saturation input field {error.args[0]!r}"
            ) from error
        ring = PolynomialRing(
            variables=tuple(ring_data["variables"]),
            characteristic=int(ring_data.get("characteristic", 0)),
            ordering=str(ring_data.get("ordering", "dp")),
        )
        presentation = ModulePresentation(
            ring=ring,
            rank=int(module_data["rank"]),
            generators=tuple(
                tuple(str(entry) for entry in generator)
                for generator in module_data.get("generators", ())
            ),
            label=str(module_data.get("label", "module")),
            singular_setup=module_data.get("singular_setup"),
        )
        assurance_data = data.get("assurance", {})
        assurance = CertificateAssurance(
            claim=assurance_data.get("claim", "auto"),
            target_characteristic=assurance_data.get(
                "target_characteristic"
            ),
            note=assurance_data.get("note"),
        )
        distinguished = data.get("distinguished_class")
        return cls(
            presentation=presentation,
            support_ideal=tuple(str(item) for item in support_ideal),
            completion_ideal=tuple(
                str(item) for item in completion_ideal
            ),
            parameter_base_variables=tuple(
                str(item) for item in parameter_variables
            ),
            normal_variables=tuple(str(item) for item in normal_variables),
            jet_orders=tuple(int(item) for item in data.get("jet_orders", ())),
            distinguished_class=(
                tuple(str(item) for item in distinguished)
                if distinguished is not None
                else None
            ),
            jet_strategy=data.get("jet_strategy", "full_saturation"),
            transition_annihilator=data.get("transition_annihilator"),
            assurance=assurance,
        )

    def to_mapping(self) -> dict[str, Any]:
        """Return the stable JSON-compatible problem description."""

        presentation = {
            "rank": self.presentation.rank,
            "generators": [
                list(generator)
                for generator in self.presentation.generators
            ],
            "label": self.presentation.label,
            "singular_setup": self.presentation.singular_setup,
        }
        return {
            "schema": self.schema,
            "ring": {
                "variables": list(self.presentation.ring.variables),
                "characteristic": self.presentation.ring.characteristic,
                "ordering": self.presentation.ring.ordering,
            },
            "module_presentation": presentation,
            "support_ideal": list(self.support_ideal),
            "completion_ideal": list(self.completion_ideal),
            "parameter_base_variables": list(
                self.parameter_base_variables
            ),
            "normal_variables": list(self.normal_variables),
            "jet_orders": list(self.jet_orders),
            "distinguished_class": (
                list(self.distinguished_class)
                if self.distinguished_class is not None
                else None
            ),
            "jet_strategy": self.jet_strategy,
            "transition_annihilator": self.transition_annihilator,
            "assurance": asdict(self.assurance),
        }


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


def _certificate_state(
    ring: PolynomialRing,
    assurance: CertificateAssurance | None = None,
) -> dict[str, Any]:
    """Separate exact arithmetic from the scope assigned to its result."""

    policy = assurance or CertificateAssurance()
    claim = policy.claim
    if claim == "auto":
        claim = "exact" if ring.characteristic == 0 else "modular"
    target_characteristic = policy.target_characteristic
    if target_characteristic is None:
        target_characteristic = (
            ring.characteristic if claim == "exact" else 0
        )
    if claim == "modular" and ring.characteristic == 0:
        raise ValueError(
            "a modular certificate needs a positive backend characteristic"
        )
    if claim == "exact" and target_characteristic != ring.characteristic:
        raise ValueError(
            "an exact claim must target the backend characteristic"
        )
    return {
        "backend_arithmetic": "exact",
        "backend_characteristic": ring.characteristic,
        "claim_assurance": claim,
        "target_characteristic": target_characteristic,
        "characteristic_zero_lift": (
            "not_needed"
            if claim == "exact" and ring.characteristic == 0
            else (
                "not_claimed"
                if claim == "modular" and target_characteristic == 0
                else "outside_scope"
            )
        ),
        "note": policy.note,
    }


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

    schema = "support-saturation-certificate.v3"

    def __init__(self, options: CompilerOptions | None = None) -> None:
        self.options = options or CompilerOptions()

    def compile_problem(
        self, problem: SupportSaturationProblem
    ) -> dict[str, Any]:
        """Compile the shared, role-aware support-saturation input schema."""

        filtration = (
            NormalFiltration(
                ideal=problem.completion_ideal,
                orders=problem.jet_orders,
                strategy=problem.jet_strategy,
                transition_annihilator=problem.transition_annihilator,
            )
            if problem.completion_ideal
            else None
        )
        result = self.compile(
            problem.presentation,
            problem.support_ideal,
            distinguished_class=problem.distinguished_class,
            filtration=filtration,
            assurance=problem.assurance,
        )
        problem_data = problem.to_mapping()
        canonical_problem = json.dumps(
            problem_data, sort_keys=True, separators=(",", ":")
        )
        result["problem_sha256"] = hashlib.sha256(
            canonical_problem.encode()
        ).hexdigest()
        result["problem"] = problem_data
        result["certificate_state"] = _certificate_state(
            problem.presentation.ring, problem.assurance
        )
        return result

    def compile(
        self,
        presentation: ModulePresentation,
        boundary_ideal: Sequence[str],
        distinguished_class: Sequence[str] | None = None,
        filtration: NormalFiltration | None = None,
        assurance: CertificateAssurance | None = None,
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
            "assurance": asdict(
                assurance or CertificateAssurance()
            ),
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
        result["certificate_state"] = _certificate_state(
            presentation.ring, assurance
        )
        return result

    def compile_distinguished_jets(
        self,
        presentation: ModulePresentation,
        boundary_ideal: Sequence[str],
        distinguished_class: Sequence[str],
        filtration: NormalFiltration,
        assurance: CertificateAssurance | None = None,
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
            "assurance": asdict(
                assurance or CertificateAssurance()
            ),
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
                    "least_annihilating_exponent": (
                        integer(f"JET_{order}_BOUNDARY_EXPONENT")
                        if integer(f"JET_{order}_BOUNDARY_EXPONENT") >= 0
                        else None
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
            exponents = [
                jet["boundary_annihilation_exponent"] for jet in jets
            ]
            uniform_test = {
                "status": (
                    "certified_on_requested_finite_tower"
                    if all(exponent >= 0 for exponent in exponents)
                    else "search_bound_exhausted_on_at_least_one_jet"
                ),
                "requested_orders": list(filtration.orders),
                "least_common_exponent": (
                    max(exponents)
                    if all(exponent >= 0 for exponent in exponents)
                    else None
                ),
                "all_order_uniform_bound_certified": False,
                "warning": (
                    "a finite jet prefix is not an all-order "
                    "uniform-exponent certificate"
                ),
            }
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
            uniform_test = {
                "status": "distinguished_class_only",
                "requested_orders": list(filtration.orders),
                "least_common_exponent": None,
                "all_order_uniform_bound_certified": False,
            }
        return {
            "schema": "support-saturation-finite-jet-certificate.v1",
            "input_sha256": input_hash,
            "input": input_data,
            "certificate_state": _certificate_state(
                presentation.ring, assurance
            ),
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
                "completion_ideal": list(filtration.ideal),
                "convention": "N_n=N+m^n F",
                "scope": scope,
                "torsion_exponent_search_bound": (
                    self.options.torsion_exponent_bound
                ),
                "uniform_exponent_test": uniform_test,
                "jets": jets,
                "transitions": transitions,
            },
        }

    def compile_distinguished_support_witness(
        self,
        presentation: ModulePresentation,
        boundary_ideal: Sequence[str],
        distinguished_class: Sequence[str],
        assurance: CertificateAssurance | None = None,
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
            'LIB "primdec.lib";',
            setup,
            "module N=InputPresentation;",
            f"ideal Boundary=std({_ideal(boundary)});",
            f"vector Distinguished={_module((distinguished,))}[1];",
            'print("@@CLASS_ZERO="+string(reduce(Distinguished,N)==0));',
            (
                "module DistinguishedKernel=std("
                "modulo(module(Distinguished),N));"
            ),
            (
                "ideal DistinguishedAnnihilator=std("
                "ideal(DistinguishedKernel));"
            ),
            (
                "ideal DistinguishedAnnihilatorRadical=std("
                "radical(DistinguishedAnnihilator));"
            ),
            _section_printer(
                "DISTINGUISHED_ANNIHILATOR",
                "DistinguishedAnnihilator",
            ),
            _section_printer(
                "DISTINGUISHED_ANNIHILATOR_RADICAL",
                "DistinguishedAnnihilatorRadical",
            ),
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
            "assurance": asdict(
                assurance or CertificateAssurance()
            ),
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
            "certificate_state": _certificate_state(
                presentation.ring, assurance
            ),
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
                "least_annihilating_exponent": exponent,
                "annihilator": parsed["sections"][
                    "DISTINGUISHED_ANNIHILATOR"
                ],
                "annihilator_radical": parsed["sections"][
                    "DISTINGUISHED_ANNIHILATOR_RADICAL"
                ],
                "boundary_power_tests": tests,
            },
            "associated_primes": {
                "status": (
                    "existence_of_boundary_containing_prime_certified; "
                    "primes_not_enumerated"
                ),
                "primes": None,
            },
            "associated_prime_candidates": {
                "status": (
                    "annihilator_radical_computed_but_prime_components_"
                    "not_enumerated"
                ),
                "primes": None,
                "support_radical": parsed["sections"][
                    "DISTINGUISHED_ANNIHILATOR_RADICAL"
                ],
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
                "ideal H0AnnihilatorRadical="
                "std(radical(H0Annihilator)); "
                + _section_printer("H0_RELATIONS", "H0Relations")
                + _section_printer("H0_ANNIHILATOR", "H0Annihilator")
                + _section_printer(
                    "H0_ANNIHILATOR_RADICAL",
                    "H0AnnihilatorRadical",
                )
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
                    (
                        "ideal DistinguishedAnnihilatorRadical=std("
                        "radical(DistinguishedAnnihilator));"
                    ),
                    _section_printer(
                        "DISTINGUISHED_REMAINDER",
                        "module(DistinguishedRemainder)",
                    ),
                    _section_printer(
                        "DISTINGUISHED_ANNIHILATOR",
                        "DistinguishedAnnihilator",
                    ),
                    _section_printer(
                        "DISTINGUISHED_ANNIHILATOR_RADICAL",
                        "DistinguishedAnnihilatorRadical",
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
            "ideals": {
                "support_ideal": input_data["boundary_ideal"],
                "completion_ideal": (
                    list(filtration.ideal) if filtration is not None else []
                ),
            },
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
                "annihilator_radical": sections.get(
                    "H0_ANNIHILATOR_RADICAL", []
                ),
                "least_annihilating_exponent": (
                    integer("H0_BOUNDARY_EXPONENT")
                    if integer("H0_BOUNDARY_EXPONENT") >= 0
                    else None
                ),
                "associated_primes": h0_associated,
            },
            "associated_primes": {
                "status": associated_status,
                "primes": associated,
                "boundary_containing_primes": boundary_containing,
            },
            "associated_prime_candidates": {
                "status": (
                    "exact_boundary_containing_associated_primes"
                    if h0_associated is not None
                    else (
                        "empty_for_zero_local_cohomology"
                        if saturation_equal
                        else (
                            "annihilator_radical_computed_but_prime_"
                            "components_not_enumerated"
                        )
                    )
                ),
                "primes": (
                    h0_associated
                    if h0_associated is not None
                    else ([] if saturation_equal else None)
                ),
                "support_radical": sections.get(
                    "H0_ANNIHILATOR_RADICAL", []
                ),
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
                "annihilator_radical": sections[
                    "DISTINGUISHED_ANNIHILATOR_RADICAL"
                ],
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
                        "least_annihilating_exponent": (
                            integer(f"JET_{order}_BOUNDARY_EXPONENT")
                            if integer(
                                f"JET_{order}_BOUNDARY_EXPONENT"
                            )
                            >= 0
                            else None
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
                    "uniform_exponent_test": {
                        "status": (
                            "certified_on_requested_finite_tower"
                            if all(
                                exponent >= 0 for exponent in exponents
                            )
                            else "search_bound_exhausted_on_at_least_one_jet"
                        ),
                        "requested_orders": list(filtration.orders),
                        "least_common_exponent": (
                            max(exponents)
                            if all(
                                exponent >= 0 for exponent in exponents
                            )
                            else None
                        ),
                        "all_order_uniform_bound_certified": False,
                        "warning": (
                            "a finite jet prefix is not an all-order "
                            "uniform-exponent certificate"
                        ),
                    },
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
                "completion_ideal": list(filtration.ideal),
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
