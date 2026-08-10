"""Exact leading profiles for logarithmic differentials at surface nodes.

The local models used here are deliberately small.  They record only data
that survive modulo the maximal ideal: a transverse coefficient order and a
residue/contact index.  This is enough to certify logarithmic etaleness over
a target node, but not enough to reconstruct a determinant curve over a
smooth target-boundary point.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd


@dataclass(frozen=True)
class LogNodeProfile:
    """Certified consequences of one two-dimensional monomial leading term."""

    name: str
    target_boundary_rank: int
    transverse_order: int
    contact_index: int
    source_blowups: int
    resolved_exponent_matrix: tuple[tuple[int, int], tuple[int, int]]
    residual_determinant: int | None
    fitting_one_is_unit: bool
    cokernel_model: str
    normalization_defect: str


@dataclass(frozen=True)
class ToricFanProfile:
    """Common regular refinement for one rank-two toric exponent map."""

    exponent_map: tuple[tuple[int, int], tuple[int, int]]
    original_source_rays: tuple[tuple[int, int], ...]
    target_rays: tuple[tuple[int, int], ...]
    target_preimage_rays: tuple[tuple[int, int], ...]
    refined_source_rays: tuple[tuple[int, int], ...]
    exponent_matrices: tuple[
        tuple[tuple[int, int], tuple[int, int]], ...
    ]
    determinants: tuple[int, ...]


def _validate(contact_index: int, characteristic: int) -> None:
    if contact_index <= 0:
        raise ValueError("the contact index must be positive")
    if characteristic < 0:
        raise ValueError("the characteristic cannot be negative")
    if characteristic and contact_index % characteristic == 0:
        raise ValueError("the leading log matrix is not invertible in this characteristic")


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def primitive(ray: tuple[int, int]) -> tuple[int, int]:
    divisor = gcd(abs(ray[0]), abs(ray[1]))
    if divisor == 0:
        raise ValueError("the zero vector is not a ray")
    result = (ray[0] // divisor, ray[1] // divisor)
    if result[0] < 0 or (result[0] == 0 and result[1] < 0):
        result = (-result[0], -result[1])
    return result


def _extended_gcd(first: int, second: int) -> tuple[int, int, int]:
    if second == 0:
        sign = 1 if first >= 0 else -1
        return abs(first), sign, 0
    common, left, right = _extended_gcd(second, first % second)
    return common, right, left - (first // second) * right


def regular_cone_rays(
    left: tuple[int, int], right: tuple[int, int]
) -> tuple[tuple[int, int], ...]:
    """Return the minimal left-anchored regular subdivision of a cone."""

    left = primitive(left)
    right = primitive(right)
    cone_determinant = determinant(left, right)
    if cone_determinant <= 0:
        raise ValueError("cone rays must be positively ordered")
    if cone_determinant == 1:
        return (left, right)

    common, bezout_x, bezout_y = _extended_gcd(left[0], left[1])
    if common != 1:
        raise AssertionError("a primitive ray lost coprimality")
    # left_x*z_y-left_y*z_x=1.
    candidate = (-bezout_y, bezout_x)
    remainder = determinant(candidate, right) % cone_determinant
    if remainder == 0:
        raise AssertionError("primitive cone endpoints produced no interior ray")
    shift = (remainder - determinant(candidate, right)) // cone_determinant
    middle = (
        candidate[0] + shift * left[0],
        candidate[1] + shift * left[1],
    )
    if determinant(left, middle) != 1:
        raise AssertionError("the first subdivision cone is not regular")
    tail = regular_cone_rays(middle, right)
    return (left, *tail)


def _apply_matrix(
    matrix: tuple[tuple[int, int], tuple[int, int]],
    ray: tuple[int, int],
) -> tuple[int, int]:
    return (
        matrix[0][0] * ray[0] + matrix[0][1] * ray[1],
        matrix[1][0] * ray[0] + matrix[1][1] * ray[1],
    )


def _inside_sector(
    ray: tuple[int, int],
    left: tuple[int, int],
    right: tuple[int, int],
) -> bool:
    return determinant(left, ray) >= 0 and determinant(ray, right) >= 0


def _ray_slope(ray: tuple[int, int]) -> tuple[int, Fraction]:
    """Sort first-quadrant rays by slope without floating-point arithmetic."""

    if ray[0] == 0:
        return (1, Fraction(0))
    return (0, Fraction(ray[1], ray[0]))


def compile_toric_fan_profile(
    *,
    exponent_map: tuple[tuple[int, int], tuple[int, int]],
    source_rays: tuple[tuple[int, int], ...],
    target_rays: tuple[tuple[int, int], ...],
) -> ToricFanProfile:
    """Refine a regular source sector until it maps cone-wise to a target fan."""

    if any(
        determinant(source_rays[index], source_rays[index + 1]) != 1
        for index in range(len(source_rays) - 1)
    ):
        raise ValueError("the input source fan must be regular")
    if any(
        determinant(target_rays[index], target_rays[index + 1]) != 1
        for index in range(len(target_rays) - 1)
    ):
        raise ValueError("the target fan must be regular")
    map_determinant = determinant(exponent_map[0], exponent_map[1])
    if map_determinant <= 0:
        raise ValueError("the exponent map must preserve orientation")

    # A target ray t has rational inverse direction adj(M)*t.  Add its
    # primitive preimage whenever that ray lies in the declared source sector.
    adjugate = (
        (exponent_map[1][1], -exponent_map[0][1]),
        (-exponent_map[1][0], exponent_map[0][0]),
    )
    preimages: set[tuple[int, int]] = set()
    for target in target_rays:
        numerator = _apply_matrix(adjugate, target)
        candidate = primitive(numerator)
        if _inside_sector(candidate, source_rays[0], source_rays[-1]):
            preimages.add(candidate)

    combined = set(source_rays) | preimages
    ordered = sorted(combined, key=_ray_slope)
    refined: list[tuple[int, int]] = [ordered[0]]
    for right in ordered[1:]:
        refined.extend(regular_cone_rays(refined[-1], right)[1:])

    matrices: list[tuple[tuple[int, int], tuple[int, int]]] = []
    determinants: list[int] = []
    for source_left, source_right in zip(refined, refined[1:]):
        image_left = _apply_matrix(exponent_map, source_left)
        image_right = _apply_matrix(exponent_map, source_right)
        target_cone = next(
            (
                (target_left, target_right)
                for target_left, target_right in zip(target_rays, target_rays[1:])
                if _inside_sector(image_left, target_left, target_right)
                and _inside_sector(image_right, target_left, target_right)
            ),
            None,
        )
        if target_cone is None:
            raise ValueError("the refined source cone crosses a target fan ray")
        target_left, target_right = target_cone
        # For a regular target basis (l,r), coordinates of z are
        # (det(z,r),det(l,z)).
        left_coordinates = (
            determinant(image_left, target_right),
            determinant(target_left, image_left),
        )
        right_coordinates = (
            determinant(image_right, target_right),
            determinant(target_left, image_right),
        )
        matrix = (
            (left_coordinates[0], right_coordinates[0]),
            (left_coordinates[1], right_coordinates[1]),
        )
        matrix_determinant = determinant(matrix[0], matrix[1])
        if matrix_determinant != map_determinant:
            raise AssertionError("the local exponent determinant changed under refinement")
        matrices.append(matrix)
        determinants.append(matrix_determinant)

    return ToricFanProfile(
        exponent_map=exponent_map,
        original_source_rays=source_rays,
        target_rays=target_rays,
        target_preimage_rays=tuple(sorted(preimages, key=_ray_slope)),
        refined_source_rays=tuple(refined),
        exponent_matrices=tuple(matrices),
        determinants=tuple(determinants),
    )


def target_node_profile(
    name: str,
    *,
    transverse_order: int,
    residue_index: int,
    characteristic: int = 0,
) -> LogNodeProfile:
    """Profile ``(pi,xi)=(tau*w^p,w^e)`` over a target SNC node.

    If ``p`` is negative, successive blowups in the terminal direction replace
    ``tau`` by ``u*w**(-p)``.  The resulting exponent matrix is triangular
    with determinant ``e``.  In characteristic zero (or tame
    characteristic), the logarithmic differential is therefore invertible
    at the resolved node, independently of higher-order terms.
    """

    _validate(residue_index, characteristic)
    source_blowups = max(0, -transverse_order)
    resolved_order = transverse_order + source_blowups
    matrix = ((1, resolved_order), (0, residue_index))
    determinant = residue_index
    return LogNodeProfile(
        name=name,
        target_boundary_rank=2,
        transverse_order=transverse_order,
        contact_index=residue_index,
        source_blowups=source_blowups,
        resolved_exponent_matrix=matrix,
        residual_determinant=determinant,
        fitting_one_is_unit=True,
        cokernel_model="0",
        normalization_defect="zero",
    )


def smooth_target_boundary_profile(
    name: str,
    *,
    transverse_order: int,
    contact_index: int,
    characteristic: int = 0,
) -> LogNodeProfile:
    """Leading profile ``(pi,z)=(tau*w^p,w^e)`` over a smooth boundary.

    The target uses ``dlog(pi), dz`` whereas the source node uses
    ``dlog(tau), dlog(w)``.  The leading matrix has a unit entry and
    determinant ``e*w**e``.  Thus its cokernel is cyclic, but transverse
    higher terms can change the reduced determinant curve; its normalization
    defect is intentionally left undetermined.
    """

    _validate(contact_index, characteristic)
    source_blowups = max(0, -transverse_order)
    resolved_order = transverse_order + source_blowups
    matrix = ((1, resolved_order), (0, contact_index))
    return LogNodeProfile(
        name=name,
        target_boundary_rank=1,
        transverse_order=transverse_order,
        contact_index=contact_index,
        source_blowups=source_blowups,
        resolved_exponent_matrix=matrix,
        residual_determinant=None,
        fitting_one_is_unit=True,
        cokernel_model=f"cyclic with leading determinant w^{contact_index}",
        normalization_defect="not determined by the leading terminal row",
    )


def close_smooth_profile_by_boundary_support(
    profile: LogNodeProfile,
    *,
    terminal_generic_order: int,
) -> LogNodeProfile:
    """Close a smooth-target profile when the determinant is boundary-supported.

    In a regular local source node ``R`` with boundary ``tau*w=0``, a
    determinant whose divisor is boundary-supported is a unit times
    ``tau**a*w**b``.  If its generic order on ``tau=0`` is zero, then ``a=0``;
    the leading contact calculation forces ``b`` to be the contact index.
    A unit entry then reduces the presentation to ``diag(1,w**b)``.
    """

    if profile.target_boundary_rank != 1:
        raise ValueError("boundary-support closure applies to a smooth target boundary")
    if not profile.fitting_one_is_unit:
        raise ValueError("a cyclic presentation requires unit Fitt_1")
    if terminal_generic_order != 0:
        raise ValueError("the terminal component would remain in the determinant divisor")
    return LogNodeProfile(
        name=profile.name,
        target_boundary_rank=profile.target_boundary_rank,
        transverse_order=profile.transverse_order,
        contact_index=profile.contact_index,
        source_blowups=profile.source_blowups,
        resolved_exponent_matrix=profile.resolved_exponent_matrix,
        residual_determinant=None,
        fitting_one_is_unit=True,
        cokernel_model=f"R/(w^{profile.contact_index})",
        normalization_defect="zero on the smooth reduced support w=0",
    )
