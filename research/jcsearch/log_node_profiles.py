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


@dataclass(frozen=True)
class CyclicSNCMatchingProfile:
    """Branchwise splitting of ``R/(u^a*v^b)`` at an SNC node."""

    left_multiplicity: int
    right_multiplicity: int
    cokernel_model: str
    fitting_zero: str
    fitting_one: str
    branchwise_model: str
    matching_quotient: str
    matching_length: int
    finite_support_torsion: str


@dataclass(frozen=True)
class CyclicBoundaryCharge:
    """Intersection/Chern charge of a thickened SNC boundary divisor."""

    multiplicities: tuple[int, ...]
    self_intersections: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]
    node_matching_length: int
    doubled_charge: int
    charge: Fraction
    cyclic_cokernel_ch2: Fraction


@dataclass(frozen=True)
class CyclicCokernelTwistProfile:
    """Chern data of a cyclic rank-two cokernel from its kernel line."""

    divisor_square: Fraction
    kernel_degree: Fraction
    cokernel_line_degree: Fraction
    cokernel_ch2: Fraction


@dataclass(frozen=True)
class ContractedCyclicCokernelProfile:
    """Chern data when the kernel line defines a Gauss map to ``P^1``."""

    divisor_square: Fraction
    gauss_degree: Fraction
    kernel_degree: Fraction
    cokernel_line_degree: Fraction
    cokernel_ch2: Fraction


@dataclass(frozen=True)
class TangentialKernelTrivializationProfile:
    """A fixed target covector trivializing a cyclic logarithmic kernel."""

    divisor_square: Fraction
    determinant_multiplicities: tuple[int, ...]
    tangential_pullback_orders: tuple[int, ...]
    tangential_excess_orders: tuple[int, ...]
    kernel_degree: Fraction
    gauss_degree: Fraction
    cokernel_ch2: Fraction


@dataclass(frozen=True)
class LogarithmicCh2Budget:
    """Global Chern numbers of a boundary-supported logarithmic complex."""

    source_log_square: Fraction
    target_log_square: Fraction
    source_log_c2: Fraction
    target_log_c2: Fraction
    geometric_degree: int
    ch2: Fraction


@dataclass(frozen=True)
class LogarithmicCh2ModelChange:
    """Change of the global budget under additional boundary blowups."""

    geometric_degree: int
    source_smooth_blowups: int
    target_smooth_blowups: int
    source_node_blowups: int
    target_node_blowups: int
    ch2_change: Fraction


@dataclass(frozen=True)
class ResidualPointBudget:
    """Signed residual after certified divisorial contributions are removed.

    This is only an effective point length when a separate exact filtration
    realizes it by a finite-length sheaf.  The arithmetic alone deliberately
    makes no positivity claim.
    """

    global_ch2: Fraction
    divisorial_ch2: tuple[Fraction, ...]
    residual_ch2: Fraction


@dataclass(frozen=True)
class ImmersedAffineRowProfile:
    """Chern data of a cyclic affine-boundary row over an immersed curve.

    The reduced kernel is the pullback of the logarithmic conormal line of
    the target normalization.  The full kernel degree on the Cartier
    thickening is multiplied by the transverse index.
    """

    target_log_curve_intersection: Fraction
    normalization_log_cotangent_degree: Fraction
    residue_degree: int
    transverse_index: int
    source_self_intersection: int
    reduced_kernel_degree: Fraction
    thickened_kernel_degree: Fraction
    determinant_divisor_square: Fraction
    cokernel_ch2: Fraction


def _validate(contact_index: int, characteristic: int) -> None:
    if contact_index <= 0:
        raise ValueError("the contact index must be positive")
    if characteristic < 0:
        raise ValueError("the characteristic cannot be negative")
    if characteristic and contact_index % characteristic == 0:
        raise ValueError("the leading log matrix is not invertible in this characteristic")


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def logarithmic_ch2_budget(
    source_log_square: int | Fraction,
    target_log_square: int | Fraction,
    geometric_degree: int,
    *,
    source_log_c2: int | Fraction = 1,
    target_log_c2: int | Fraction = 1,
) -> LogarithmicCh2Budget:
    """Return ``deg ch_2(Omega_X(log)-f^*Omega_Y(log))``.

    The general surface formula is

    ``(L_X^2-d*L_Y^2-2*e_X+2*d*e_Y)/2``,

    where ``e_X`` and ``e_Y`` are the integrated top logarithmic Chern
    classes.  For SNC completions of ``A^2`` both are one.
    """

    if geometric_degree <= 0:
        raise ValueError("the geometric degree must be positive")
    source_log_square = Fraction(source_log_square)
    target_log_square = Fraction(target_log_square)
    source_log_c2 = Fraction(source_log_c2)
    target_log_c2 = Fraction(target_log_c2)
    ch2 = (
        source_log_square
        - geometric_degree * target_log_square
        - 2 * source_log_c2
        + 2 * geometric_degree * target_log_c2
    ) / 2
    return LogarithmicCh2Budget(
        source_log_square=source_log_square,
        target_log_square=target_log_square,
        source_log_c2=source_log_c2,
        target_log_c2=target_log_c2,
        geometric_degree=geometric_degree,
        ch2=ch2,
    )


def logarithmic_ch2_from_ramification(
    pullback_target_dot_ramification: int | Fraction,
    ramification_square: int | Fraction,
    geometric_degree: int,
    *,
    source_log_c2: int | Fraction = 1,
    target_log_c2: int | Fraction = 1,
) -> Fraction:
    """Evaluate the same budget from ``L_X=f^*L_Y+R_log``."""

    if geometric_degree <= 0:
        raise ValueError("the geometric degree must be positive")
    return (
        Fraction(pullback_target_dot_ramification)
        + Fraction(ramification_square) / 2
        - Fraction(source_log_c2)
        + geometric_degree * Fraction(target_log_c2)
    )


def logarithmic_ch2_model_change(
    geometric_degree: int,
    *,
    source_smooth_blowups: int = 0,
    target_smooth_blowups: int = 0,
    source_node_blowups: int = 0,
    target_node_blowups: int = 0,
) -> LogarithmicCh2ModelChange:
    """Track the budget change under ordinary boundary blowups.

    A boundary-node blowup is log crepant.  A smooth-boundary blowup lowers
    ``(K+D)^2`` by one.  Thus additional source and target smooth blowups
    change the budget by ``(d*s_Y-s_X)/2``.
    """

    if geometric_degree <= 0:
        raise ValueError("the geometric degree must be positive")
    counts = (
        source_smooth_blowups,
        target_smooth_blowups,
        source_node_blowups,
        target_node_blowups,
    )
    if any(count < 0 for count in counts):
        raise ValueError("boundary blowup counts must be nonnegative")
    change = Fraction(
        geometric_degree * target_smooth_blowups - source_smooth_blowups,
        2,
    )
    return LogarithmicCh2ModelChange(
        geometric_degree=geometric_degree,
        source_smooth_blowups=source_smooth_blowups,
        target_smooth_blowups=target_smooth_blowups,
        source_node_blowups=source_node_blowups,
        target_node_blowups=target_node_blowups,
        ch2_change=change,
    )


def residual_point_budget(
    global_ch2: int | Fraction,
    divisorial_ch2: tuple[int | Fraction, ...],
) -> ResidualPointBudget:
    """Subtract certified divisorial packets without assuming effectivity."""

    global_ch2 = Fraction(global_ch2)
    contributions = tuple(Fraction(value) for value in divisorial_ch2)
    return ResidualPointBudget(
        global_ch2=global_ch2,
        divisorial_ch2=contributions,
        residual_ch2=global_ch2 - sum(contributions, Fraction(0)),
    )


def immersed_affine_row_profile(
    target_log_curve_intersection: int | Fraction,
    normalization_log_cotangent_degree: int | Fraction,
    residue_degree: int,
    transverse_index: int,
    source_self_intersection: int,
) -> ImmersedAffineRowProfile:
    """Return the cyclic ``ch_2`` of an immersed affine-boundary row.

    Let ``nu:C^nu->(Y,D_Y)`` be a logarithmic immersion and let a reduced
    source boundary component ``E`` map to ``C^nu`` with residue degree
    ``f``.  The conormal sequence gives

    ``deg_E K_red = f*(L_Y.C-deg Omega_Cnu(log S))``.

    If the transverse ramification index is ``e``, the determinant packet is
    ``D=eE``.  When ``Fitt_1`` is the unit ideal on that packet and the kernel
    line on ``D`` reduces to ``K_red``, the Cartier filtration gives
    ``deg_D K=e*deg_E K_red`` and GRR gives

    ``ch_2 = deg_D K + D^2/2``.

    Point-supported failures of immersion or cyclicity must be entered in a
    separate exact filtration; this helper deliberately does not absorb
    them.
    """

    if residue_degree <= 0:
        raise ValueError("the residue degree must be positive")
    if transverse_index <= 0:
        raise ValueError("the transverse index must be positive")
    if source_self_intersection >= 0:
        raise ValueError("a new exceptional boundary component must be negative")
    target_degree = Fraction(target_log_curve_intersection)
    normalization_degree = Fraction(normalization_log_cotangent_degree)
    reduced_kernel_degree = residue_degree * (
        target_degree - normalization_degree
    )
    thickened_kernel_degree = transverse_index * reduced_kernel_degree
    divisor_square = Fraction(
        transverse_index * transverse_index * source_self_intersection
    )
    cokernel_ch2 = thickened_kernel_degree + divisor_square / 2
    return ImmersedAffineRowProfile(
        target_log_curve_intersection=target_degree,
        normalization_log_cotangent_degree=normalization_degree,
        residue_degree=residue_degree,
        transverse_index=transverse_index,
        source_self_intersection=source_self_intersection,
        reduced_kernel_degree=reduced_kernel_degree,
        thickened_kernel_degree=thickened_kernel_degree,
        determinant_divisor_square=divisor_square,
        cokernel_ch2=cokernel_ch2,
    )


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


def cyclic_snc_matching_profile(
    left_multiplicity: int,
    right_multiplicity: int,
    *,
    left_parameter: str = "u",
    right_parameter: str = "v",
) -> CyclicSNCMatchingProfile:
    """Return the canonical branch-splitting quotient of a cyclic node.

    In a two-dimensional regular local ring, the principal ideals
    ``(u**a)`` and ``(v**b)`` have intersection ``(u**a*v**b)``.  The
    resulting Mayer--Vietoris sequence has finite quotient
    ``R/(u**a,v**b)`` of length ``a*b``.  The cyclic hypersurface module is
    Cohen--Macaulay, so it has no finite-support submodule; the quotient is a
    degree-one matching defect, not degree-zero torsion.
    """

    if left_multiplicity <= 0 or right_multiplicity <= 0:
        raise ValueError("SNC branch multiplicities must be positive")
    left_power = f"{left_parameter}^{left_multiplicity}"
    right_power = f"{right_parameter}^{right_multiplicity}"
    product = f"{left_power}*{right_power}"
    return CyclicSNCMatchingProfile(
        left_multiplicity=left_multiplicity,
        right_multiplicity=right_multiplicity,
        cokernel_model=f"R/({product})",
        fitting_zero=f"({product})",
        fitting_one="R",
        branchwise_model=f"R/({left_power}) direct-sum R/({right_power})",
        matching_quotient=f"R/({left_power},{right_power})",
        matching_length=left_multiplicity * right_multiplicity,
        finite_support_torsion="zero",
    )


def cyclic_boundary_charge(
    multiplicities: tuple[int, ...],
    self_intersections: tuple[int, ...],
    edges: tuple[tuple[int, int], ...],
) -> CyclicBoundaryCharge:
    """Compute ``D^2/2`` from component and node data for ``D=sum m_i*C_i``.

    For an SNC divisor without triple intersections,

    ``D^2/2 = sum_i m_i^2*C_i^2/2 + sum_{i--j} m_i*m_j``.

    The second sum is precisely the total length of the canonical branchwise
    splitting quotients of the cyclic modules ``R/(u^m_i*v^m_j)``.  For the
    untwisted cyclic cokernel ``O_D``, ``ch_2(O_D)=-D^2/2``.
    """

    if len(multiplicities) != len(self_intersections):
        raise ValueError("multiplicity and self-intersection vectors must agree")
    if any(multiplicity <= 0 for multiplicity in multiplicities):
        raise ValueError("boundary multiplicities must be positive")
    normalized_edges: list[tuple[int, int]] = []
    for left, right in edges:
        if left == right:
            raise ValueError("a boundary edge cannot be a loop")
        if not (0 <= left < len(multiplicities) and 0 <= right < len(multiplicities)):
            raise ValueError("a boundary edge index is out of range")
        normalized_edges.append(tuple(sorted((left, right))))
    if len(set(normalized_edges)) != len(normalized_edges):
        raise ValueError("boundary edges must be distinct")

    node_length = sum(
        multiplicities[left] * multiplicities[right]
        for left, right in normalized_edges
    )
    doubled_charge = sum(
        multiplicity * multiplicity * self_intersection
        for multiplicity, self_intersection in zip(
            multiplicities, self_intersections
        )
    ) + 2 * node_length
    charge = Fraction(doubled_charge, 2)
    return CyclicBoundaryCharge(
        multiplicities=multiplicities,
        self_intersections=self_intersections,
        edges=tuple(normalized_edges),
        node_matching_length=node_length,
        doubled_charge=doubled_charge,
        charge=charge,
        cyclic_cokernel_ch2=-charge,
    )


def cyclic_cokernel_twist_profile(
    divisor_square: int | Fraction,
    kernel_degree: int | Fraction,
) -> CyclicCokernelTwistProfile:
    """Use ``L=K tensor O_D(D)`` for a cyclic rank-two cokernel ``i_*L``."""

    divisor_square = Fraction(divisor_square)
    kernel_degree = Fraction(kernel_degree)
    cokernel_line_degree = kernel_degree + divisor_square
    cokernel_ch2 = cokernel_line_degree - divisor_square / 2
    return CyclicCokernelTwistProfile(
        divisor_square=divisor_square,
        kernel_degree=kernel_degree,
        cokernel_line_degree=cokernel_line_degree,
        cokernel_ch2=cokernel_ch2,
    )


def contracted_cyclic_cokernel_profile(
    divisor_square: int | Fraction,
    gauss_degree: int | Fraction,
) -> ContractedCyclicCokernelProfile:
    """Use ``K=gamma^*O(-1)`` on a determinant packet contracted to a point.

    The degree is that of the globally generated line
    ``gamma^*O(1)=K^(-1)`` on the possibly reducible or nonreduced Cartier
    curve.  It must therefore be nonnegative.
    """

    divisor_square = Fraction(divisor_square)
    gauss_degree = Fraction(gauss_degree)
    if gauss_degree < 0:
        raise ValueError("the kernel Gauss degree must be nonnegative")
    twist = cyclic_cokernel_twist_profile(divisor_square, -gauss_degree)
    return ContractedCyclicCokernelProfile(
        divisor_square=divisor_square,
        gauss_degree=gauss_degree,
        kernel_degree=twist.kernel_degree,
        cokernel_line_degree=twist.cokernel_line_degree,
        cokernel_ch2=twist.cokernel_ch2,
    )


def tangential_kernel_trivialization_profile(
    divisor_square: int | Fraction,
    determinant_multiplicities: tuple[int, ...],
    tangential_pullback_orders: tuple[int, ...],
) -> TangentialKernelTrivializationProfile:
    """Certify a constant kernel when ``div(f^*z)`` contains ``D``.

    The entries record orders along the irreducible boundary components of
    ``D``.  If the pullback order of one fixed target tangential coordinate
    is at least the determinant multiplicity on every component, then its
    logarithmic differential is divisible by the full Cartier ideal of
    ``D``.  The nonzero target covector ``dz`` trivializes the restricted
    kernel, so its Gauss degree is zero.
    """

    if len(determinant_multiplicities) != len(tangential_pullback_orders):
        raise ValueError("determinant multiplicities and pullback orders must agree")
    if not determinant_multiplicities:
        raise ValueError("the determinant packet must have a component")
    if any(multiplicity <= 0 for multiplicity in determinant_multiplicities):
        raise ValueError("determinant multiplicities must be positive")
    excess = tuple(
        order - multiplicity
        for order, multiplicity in zip(
            tangential_pullback_orders, determinant_multiplicities
        )
    )
    if any(value < 0 for value in excess):
        raise ValueError("the tangential pullback divisor does not contain D")
    contracted = contracted_cyclic_cokernel_profile(divisor_square, 0)
    return TangentialKernelTrivializationProfile(
        divisor_square=contracted.divisor_square,
        determinant_multiplicities=determinant_multiplicities,
        tangential_pullback_orders=tangential_pullback_orders,
        tangential_excess_orders=excess,
        kernel_degree=contracted.kernel_degree,
        gauss_degree=contracted.gauss_degree,
        cokernel_ch2=contracted.cokernel_ch2,
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
