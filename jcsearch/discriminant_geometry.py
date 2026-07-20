"""Dual-curve geometry of the discriminant normalization."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import factorial

import sympy as sp


@dataclass(frozen=True)
class UniversalPrimitive:
    """A universal polynomial graph and the coefficient open set it uses."""

    degree: int
    variable: sp.Symbol
    polynomial: sp.Expr
    coefficients: tuple[sp.Symbol, ...]
    open_factor: sp.Expr
    admissible: bool


@dataclass(frozen=True)
class SaturatedIncidence:
    """Equations for an incidence stratum with one Rabinowitsch gate."""

    name: str
    degree: int
    equations: tuple[sp.Expr, ...]
    marked_points: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    saturation_factor: sp.Expr
    gate: sp.Symbol
    contact_partition: tuple[int, ...]

    @property
    def saturated_equations(self):
        """Equations defining the saturation in an extended polynomial ring."""
        return self.equations + (
            1 - self.gate * self.saturation_factor,
        )

    @property
    def elimination_variables(self):
        """Variables removed to obtain the closed coefficient-space locus."""
        return self.marked_points + (self.gate,)

    @property
    def expected_incidence_dimension(self):
        """Contact-factor dimension before projection to coefficient space."""
        return contact_incidence_dimension(self.degree, self.contact_partition)

    def eliminate(self):
        """Compute generators of the coefficient-space elimination ideal."""
        return incidence_elimination_generators(self)


@dataclass(frozen=True)
class ContactPartitionIncidence:
    """Exact root and coefficient data for an arbitrary contact partition."""

    degree: int
    multiplicities: tuple[int, ...]
    full_contact: bool
    variable: sp.Symbol
    marked_roots: tuple[sp.Symbol, ...]
    symmetry_blocks: tuple[tuple[int, ...], ...]
    symmetry_order: int
    quotient_coordinates: tuple[sp.Symbol, ...]
    root_to_quotient: tuple[tuple[sp.Symbol, sp.Expr], ...]
    residual_coefficients: tuple[sp.Symbol, ...]
    contact_polynomial: sp.Expr
    residual_polynomial: sp.Expr
    M: sp.Expr
    phi: sp.Expr
    scale_denominator: sp.Expr
    normalization: sp.Expr
    scale: sp.Expr
    slope: sp.Expr
    intercept: sp.Expr
    primitive: sp.Expr
    distinct_root_factor: sp.Expr
    residual_nonvanishing_factor: sp.Expr
    weighted_admissibility_factor: sp.Expr
    coefficient_parameters: tuple[sp.Symbol, ...]
    normalized_universal_primitive: sp.Expr
    coefficient_space_elimination_ideal: SaturatedIncidence

    @property
    def expected_seed_dimension(self):
        """Expected dimension in the normalized seed coefficient space."""
        internal = len(self.quotient_coordinates) + len(self.residual_coefficients)
        return internal - 1

    @property
    def Phi(self):
        """Alias matching the conventional notation Phi_lambda."""
        return self.phi

    @property
    def a(self):
        """Scale in ``H-sW+t=a*M``."""
        return self.scale

    @property
    def s(self):
        """Slope parameter in the inverse pencil."""
        return self.slope

    @property
    def t(self):
        """Intercept parameter in the inverse pencil."""
        return self.intercept

    @property
    def H(self):
        """Normalized admissible primitive."""
        return self.primitive


@dataclass(frozen=True)
class MultipleOmissionIncidence:
    """A common normalized seed with several full-contact omitted values."""

    degree: int
    partitions: tuple[tuple[int, ...], ...]
    contacts: tuple[ContactPartitionIncidence, ...]
    coefficient_parameters: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    marked_points: tuple[sp.Symbol, ...]
    saturation_factor: sp.Expr
    gate: sp.Symbol
    omission_coordinates: tuple[tuple[sp.Expr, sp.Expr], ...]
    pairwise_value_differences: tuple[tuple[sp.Expr, sp.Expr], ...]
    distinct_value_gates: tuple[tuple[sp.Symbol, sp.Symbol], ...]

    @property
    def saturated_equations(self):
        """Equations for simultaneous admissible full-contact incidences."""
        return self.equations + (1 - self.gate * self.saturation_factor,)

    @property
    def common_value_ideal(self):
        """Pairwise ideals whose vanishing means the omitted values coincide."""
        return self.pairwise_value_differences

    @property
    def distinct_omitted_value_equations(self):
        """Rabinowitsch equations forcing every pair of values to differ."""
        gates = tuple(
            1 - gate_s * delta_s - gate_t * delta_t
            for (gate_s, gate_t), (delta_s, delta_t) in zip(
                self.distinct_value_gates, self.pairwise_value_differences
            )
        )
        return self.saturated_equations + gates


@dataclass(frozen=True)
class TwoOmissionIncidence:
    """Uniform affine-difference presentation for two contact partitions."""

    degree: int
    partitions: tuple[tuple[int, ...], tuple[int, ...]]
    variable: sp.Symbol
    left: ContactPartitionIncidence
    right: ContactPartitionIncidence
    left_scale: sp.Symbol
    right_scale: sp.Symbol
    affine_slope: sp.Symbol
    affine_intercept: sp.Symbol
    affine_difference: sp.Expr
    equations: tuple[sp.Expr, ...]
    collision_discriminant: sp.Expr
    saturation_factor: sp.Expr
    shared_seed_incidence: MultipleOmissionIncidence

    @property
    def high_coefficient_equations(self):
        """Equations saying the difference has degree at most one."""
        polynomial = sp.Poly(self.affine_difference, self.variable)
        return tuple(
            polynomial.coeff_monomial(self.variable**power)
            for power in range(2, self.degree + 1)
        )


@dataclass(frozen=True)
class MaximalContactPhi:
    """The quotient-coordinate hypersurface for a maximal 2/3 partition."""

    double_root_count: int
    triple_root_count: int
    degree: int
    variable: sp.Symbol
    double_root_polynomial: sp.Expr
    triple_root_polynomial: sp.Expr
    double_coordinates: tuple[sp.Symbol, ...]
    triple_coordinates: tuple[sp.Symbol, ...]
    M: sp.Expr
    phi: sp.Expr

    @property
    def partition(self):
        """The corresponding maximal full-contact partition."""
        return (3,) * self.triple_root_count + (2,) * self.double_root_count


def discriminant_param(H, W):
    """Return the tangent-line coordinates (s,t) of the graph of H."""
    derivative = sp.diff(H, W)
    return sp.expand(derivative), sp.expand(W * derivative - H)


def cusp_polynomial(H, W):
    """Return the ramification polynomial of the discriminant normalization."""
    return sp.expand(sp.diff(H, W, 2))


def bitangent_equations(H, W, r, u):
    """Return divided equations for two graph points sharing a tangent line."""
    derivative = sp.diff(H, W)
    slope_equation = sp.cancel(
        (derivative.subs(W, r) - derivative.subs(W, u)) / (r - u)
    )
    intercept_equation = sp.cancel(
        (
            r * derivative.subs(W, r)
            - H.subs(W, r)
            - u * derivative.subs(W, u)
            + H.subs(W, u)
        )
        / (r - u)
    )
    return sp.expand(slope_equation), sp.expand(intercept_equation)


def symmetric_bitangent_equations(H, W, pair_sum, pair_product):
    """Return bitangent equations on unordered pairs of normalization points."""
    r, u = sp.symbols("_bitangent_r _bitangent_u")
    equations = bitangent_equations(H, W, r, u)
    result = []
    for equation in equations:
        symmetric, remainder, mapping = sp.symmetrize(
            equation, [r, u], formal=True
        )
        assert remainder == 0
        result.append(
            sp.expand(
                symmetric.subs(
                    {
                        mapping[0][0]: pair_sum,
                        mapping[1][0]: pair_product,
                    }
                )
            )
        )
    return tuple(result)


def ordinary_cusp_determinant(H, W):
    """Return det(nu''(W),nu'''(W)) modulo the cusp equation H''=0."""
    slope, intercept = discriminant_param(H, W)
    second = sp.Matrix([sp.diff(slope, W, 2), sp.diff(intercept, W, 2)])
    third = sp.Matrix([sp.diff(slope, W, 3), sp.diff(intercept, W, 3)])
    determinant = sp.expand(sp.det(sp.Matrix.hstack(second, third)))
    remainder = sp.rem(determinant, cusp_polynomial(H, W), W)
    return sp.factor(remainder)


def contact_incidence_dimension(degree: int, multiplicities):
    """Dimension bound for a common-tangent contact incidence modulo lines.

    A tangent line with distinct contact points of multiplicities ``m_i`` has
    ``H-line = product((W-r_i)**m_i) * Q``.  The returned dimension counts the
    contact points and the coefficients of Q.  ``None`` means that the total
    contact exceeds the polynomial degree, so the incidence is empty.
    """
    multiplicities = tuple(int(value) for value in multiplicities)
    if degree < 2 or not multiplicities or any(value < 2 for value in multiplicities):
        raise ValueError("contact multiplicities must all be at least two")
    total_contact = sum(multiplicities)
    if total_contact > degree:
        return None
    quotient_coefficients = degree - total_contact + 1
    return len(multiplicities) + quotient_coefficients


def collision_partition_blocks(finer, coarser):
    """Return blocks witnessing ``finer <= coarser`` under root collision.

    The relation means that every part of ``coarser`` is the sum of one block
    of parts of ``finer``.  ``None`` is returned when no such merging exists.
    Parts and blocks use deterministic decreasing order.
    """
    finer = tuple(sorted((int(value) for value in finer), reverse=True))
    coarser = tuple(sorted((int(value) for value in coarser), reverse=True))
    if (
        not finer
        or not coarser
        or any(value < 2 for value in finer + coarser)
        or sum(finer) != sum(coarser)
        or len(finer) < len(coarser)
    ):
        return None

    def search(remaining, target_index):
        if target_index == len(coarser):
            return () if not remaining else None
        target = coarser[target_index]
        for size in range(1, len(remaining) + 1):
            for block in combinations(remaining, size):
                if sum(finer[index] for index in block) != target:
                    continue
                block_set = set(block)
                tail = search(
                    tuple(index for index in remaining if index not in block_set),
                    target_index + 1,
                )
                if tail is not None:
                    return (tuple(block),) + tail
        return None

    return search(tuple(range(len(finer))), 0)


def collision_precedes(finer, coarser):
    """Return whether the coarser partition is obtained by merging parts."""
    return collision_partition_blocks(finer, coarser) is not None


def maximal_two_three_partitions(degree: int):
    """Return types indexing maximal strata: partitions using only 2 and 3.

    These are the minimal elements for ``collision_precedes`` because that
    order points from a finer partition toward a coarser collision.
    """
    degree = int(degree)
    if degree < 2:
        return ()
    partitions = []
    for threes in range(degree // 3 + 1):
        remainder = degree - 3 * threes
        if remainder % 2 == 0:
            partitions.append((3,) * threes + (2,) * (remainder // 2))
    return tuple(sorted(partitions, reverse=True))


def contact_semigroup_atoms(minimum_multiplicity=2):
    """Atoms of the additive contact semigroup ``{m >= r}``.

    An allowed multiplicity is indecomposable into two allowed summands
    exactly when it lies in ``r,...,2*r-1``.
    """
    minimum_multiplicity = int(minimum_multiplicity)
    if minimum_multiplicity < 2:
        raise ValueError("the contact threshold must be at least two")
    return tuple(range(minimum_multiplicity, 2 * minimum_multiplicity))


def contact_atom_refinement(partition, minimum_multiplicity=2):
    """Split every contact part into atoms above the chosen threshold."""
    minimum_multiplicity = int(minimum_multiplicity)
    atoms = set(contact_semigroup_atoms(minimum_multiplicity))
    partition = tuple(int(value) for value in partition)
    if not partition or any(value < minimum_multiplicity for value in partition):
        raise ValueError("every part must meet the contact threshold")
    refinement = []
    for value in partition:
        while value >= 2 * minimum_multiplicity:
            refinement.append(minimum_multiplicity)
            value -= minimum_multiplicity
        refinement.append(value)
    assert set(refinement) <= atoms
    return tuple(sorted(refinement, reverse=True))


def multiplicity_excess(partition, baseline=2):
    """Return ``sum(m_i-baseline)`` for a full-contact partition."""
    baseline = int(baseline)
    partition = tuple(int(value) for value in partition)
    if baseline < 2 or not partition or any(value < baseline for value in partition):
        raise ValueError("every part must meet the multiplicity baseline")
    return sum(value - baseline for value in partition)


def maximal_two_three_phi(double_root_count, triple_root_count, prefix="maxphi"):
    """Build ``Phi`` for ``M=Q_2**2*Q_3**3`` in quotient coordinates."""
    double_root_count = int(double_root_count)
    triple_root_count = int(triple_root_count)
    if double_root_count < 0 or triple_root_count < 0:
        raise ValueError("root counts must be nonnegative")
    if double_root_count + triple_root_count == 0:
        raise ValueError("at least one contact root is required")
    W = sp.symbols(f"_{prefix}_W")
    double_coordinates = tuple(
        sp.symbols(f"_{prefix}_q0:{double_root_count}")
    )
    triple_coordinates = tuple(
        sp.symbols(f"_{prefix}_r0:{triple_root_count}")
    )
    double_root_polynomial = sp.expand(
        W**double_root_count
        + sum(
            coefficient * W**index
            for index, coefficient in enumerate(double_coordinates)
        )
    )
    triple_root_polynomial = sp.expand(
        W**triple_root_count
        + sum(
            coefficient * W**index
            for index, coefficient in enumerate(triple_coordinates)
        )
    )
    M = sp.expand(double_root_polynomial**2 * triple_root_polynomial**3)
    phi = sp.expand(
        M.subs(W, 1) - M.subs(W, 0) - sp.diff(M, W).subs(W, 0)
    )
    return MaximalContactPhi(
        double_root_count=double_root_count,
        triple_root_count=triple_root_count,
        degree=2 * double_root_count + 3 * triple_root_count,
        variable=W,
        double_root_polynomial=double_root_polynomial,
        triple_root_polynomial=triple_root_polynomial,
        double_coordinates=double_coordinates,
        triple_coordinates=triple_coordinates,
        M=M,
        phi=phi,
    )


def affine_difference_mason_defect(left_partition, right_partition):
    """Return the strict Mason obstruction ``n-l(left)-l(right)``.

    A positive value rules out two coprime monic full-contact polynomials with
    these partitions whose difference has degree at most one.
    """
    left = tuple(int(value) for value in left_partition)
    right = tuple(int(value) for value in right_partition)
    if (
        not left
        or not right
        or any(value < 2 for value in left + right)
        or sum(left) != sum(right)
    ):
        raise ValueError("partitions must have the same degree and parts at least two")
    return sum(left) - len(left) - len(right)


def normalized_seed_space_dimension(degree: int):
    """Dimension of the normalized admissible seed space A_n."""
    if degree < 3:
        raise ValueError("inverse degree must be at least three")
    return degree - 3


def exceptional_stratum_dimension(multiplicities):
    """Dimension of a nonempty full-contact stratum E_lambda."""
    multiplicities = tuple(int(value) for value in multiplicities)
    if not multiplicities or any(value < 2 for value in multiplicities):
        raise ValueError("partition parts must all be at least two")
    return len(multiplicities) - 1


def exceptional_stratum_codimension(degree: int, multiplicities):
    """Codimension of E_lambda in the normalized seed space A_n."""
    multiplicities = tuple(int(value) for value in multiplicities)
    if sum(multiplicities) != degree:
        raise ValueError("a full-contact partition must sum to the degree")
    return normalized_seed_space_dimension(degree) - exceptional_stratum_dimension(
        multiplicities
    )


def weighted_power_sums(multiplicities, roots, count=None):
    """Return p_j=sum_i lambda_i*r_i**j for the weighted root multiset."""
    multiplicities = tuple(int(value) for value in multiplicities)
    roots = tuple(roots)
    if len(multiplicities) != len(roots):
        raise ValueError("one multiplicity is required for every root")
    if count is None:
        count = len(roots)
    return tuple(
        sp.expand(sum(weight * root**power for weight, root in zip(multiplicities, roots)))
        for power in range(1, int(count) + 1)
    )


def weighted_newton_top_coefficients(multiplicities, roots, count=None):
    """Recover the top monic coefficients from weighted Newton sums."""
    powers = weighted_power_sums(multiplicities, roots, count)
    coefficients = [sp.Integer(1)]
    for index, power_sum in enumerate(powers, start=1):
        coefficient = -(
            power_sum
            + sum(
                coefficients[offset] * powers[index - offset - 1]
                for offset in range(1, index)
            )
        ) / index
        coefficients.append(sp.expand(coefficient))
    return tuple(coefficients[1:])


def weighted_vandermonde_determinant(multiplicities, roots):
    """Closed Jacobian determinant for the top coefficient map.

    For ``M=prod(W-r_i)**lambda_i`` and top monic coefficients
    ``(c_1,...,c_l)``, this is ``det(d c_j / d r_i)``.
    """
    multiplicities = tuple(int(value) for value in multiplicities)
    roots = tuple(roots)
    if len(multiplicities) != len(roots):
        raise ValueError("one multiplicity is required for every root")
    vandermonde = sp.prod(
        roots[right] - roots[left]
        for left in range(len(roots))
        for right in range(left + 1, len(roots))
    )
    return sp.expand((-1) ** len(roots) * sp.prod(multiplicities) * vandermonde)


def contact_partition_incidence(degree, multiplicities, full_contact=True, _prefix=None):
    """Return the exact incidence attached to an arbitrary contact partition.

    Equal multiplicities are quotiented before elimination by replacing their
    marked roots with the elementary coefficients of the corresponding monic
    root polynomial.  The normalization is ``c=1``, so the returned primitive
    satisfies ``H'(1)=-1``.
    """
    degree = int(degree)
    multiplicities = tuple(sorted((int(value) for value in multiplicities), reverse=True))
    if degree < 3 or not multiplicities or any(value < 2 for value in multiplicities):
        raise ValueError("contact multiplicities must be nonempty and at least two")
    total_contact = sum(multiplicities)
    if total_contact > degree:
        raise ValueError("total contact cannot exceed the polynomial degree")
    if full_contact and total_contact != degree:
        raise ValueError("a full-contact partition must sum to the degree")

    partition_tag = "_".join(str(value) for value in multiplicities)
    prefix = _prefix or f"_cp{degree}_{partition_tag}"
    W = sp.symbols(f"{prefix}_W")
    marked_roots = tuple(sp.symbols(f"{prefix}_r0:{len(multiplicities)}"))

    contact_polynomial = sp.Integer(1)
    distinct_factors = []
    quotient_coordinates = []
    root_to_quotient = []
    symmetry_blocks = []
    block_polynomials = []
    symmetry_order = 1
    offset = 0
    for multiplicity in sorted(set(multiplicities), reverse=True):
        count = multiplicities.count(multiplicity)
        indices = tuple(range(offset, offset + count))
        offset += count
        symmetry_blocks.append(indices)
        symmetry_order *= factorial(count)
        block_coordinates = tuple(
            sp.symbols(f"{prefix}_e{multiplicity}_{index}")
            for index in range(1, count + 1)
        )
        quotient_coordinates.extend(block_coordinates)
        block_polynomial = W**count + sum(
            (-1) ** index * coordinate * W ** (count - index)
            for index, coordinate in enumerate(block_coordinates, start=1)
        )
        block_polynomial = sp.expand(block_polynomial)
        block_polynomials.append(block_polynomial)
        contact_polynomial *= block_polynomial**multiplicity

        raw_block = sp.Poly(
            sp.prod(W - marked_roots[index] for index in indices), W
        )
        for index, coordinate in enumerate(block_coordinates, start=1):
            elementary = (-1) ** index * raw_block.coeff_monomial(W ** (count - index))
            root_to_quotient.append((coordinate, sp.expand(elementary)))
        if count > 1:
            distinct_factors.append(sp.discriminant(block_polynomial, W))

    for left_index, left in enumerate(block_polynomials):
        for right in block_polynomials[left_index + 1 :]:
            distinct_factors.append(sp.resultant(left, right, W))

    residual_degree = degree - total_contact
    if full_contact:
        residual_coefficients = ()
        residual_polynomial = sp.Integer(1)
    else:
        residual_coefficients = tuple(sp.symbols(f"{prefix}_q0:{residual_degree}"))
        residual_polynomial = W**residual_degree + sum(
            coefficient * W**index
            for index, coefficient in enumerate(residual_coefficients)
        )
        residual_polynomial = sp.expand(residual_polynomial)

    contact_polynomial = sp.expand(contact_polynomial)
    M = sp.expand(contact_polynomial * residual_polynomial)
    derivative = sp.diff(M, W)
    derivative_at_zero = derivative.subs(W, 0)
    scale_denominator = sp.expand(derivative.subs(W, 1) - derivative_at_zero)
    phi = sp.expand(M.subs(W, 1) - M.subs(W, 0) - derivative_at_zero)
    normalization = sp.Integer(1)
    scale = sp.cancel(-normalization / scale_denominator)
    slope = sp.cancel(normalization * derivative_at_zero / scale_denominator)
    intercept = sp.cancel(-normalization * M.subs(W, 0) / scale_denominator)
    primitive = sp.cancel(scale * M + slope * W - intercept)

    distinct_root_factor = sp.factor(sp.prod(distinct_factors))
    if residual_degree:
        residual_nonvanishing_factor = sp.factor(
            sp.prod(
                sp.resultant(block_polynomial, residual_polynomial, W)
                for block_polynomial in block_polynomials
            )
        )
    else:
        residual_nonvanishing_factor = sp.Integer(1)
    weighted_admissibility_factor = sp.factor(
        distinct_root_factor
        * residual_nonvanishing_factor
        * scale_denominator
        * (sp.diff(M, W, 2).subs(W, 1) - 2 * scale_denominator)
    )

    coefficient_parameters = tuple(
        sp.symbols(f"{prefix}_h{index}") for index in range(3, degree)
    )
    top_coefficient = sp.cancel(
        (
            -normalization
            - sum(
                (index - 2) * coefficient
                for index, coefficient in zip(
                    range(3, degree), coefficient_parameters
                )
            )
        )
        / (degree - 2)
    )
    all_coefficients = coefficient_parameters + (top_coefficient,)
    normalized_universal_primitive = sp.expand(
        sum(
            coefficient * (W**index - W**2)
            for index, coefficient in zip(range(3, degree + 1), all_coefficients)
        )
    )
    coefficient_equations = (phi,) + tuple(
        sp.expand(
            scale_denominator * coefficient
            + normalization * sp.Poly(M, W).coeff_monomial(W**index)
        )
        for index, coefficient in zip(range(3, degree), coefficient_parameters)
    )
    gate = sp.symbols(f"{prefix}_coefficient_gate")
    coefficient_space_elimination_ideal = SaturatedIncidence(
        name=f"contact_partition_{partition_tag}",
        degree=degree,
        equations=coefficient_equations,
        marked_points=tuple(quotient_coordinates) + residual_coefficients,
        parameters=coefficient_parameters,
        saturation_factor=weighted_admissibility_factor,
        gate=gate,
        contact_partition=multiplicities,
    )

    return ContactPartitionIncidence(
        degree=degree,
        multiplicities=multiplicities,
        full_contact=bool(full_contact),
        variable=W,
        marked_roots=marked_roots,
        symmetry_blocks=tuple(symmetry_blocks),
        symmetry_order=symmetry_order,
        quotient_coordinates=tuple(quotient_coordinates),
        root_to_quotient=tuple(root_to_quotient),
        residual_coefficients=residual_coefficients,
        contact_polynomial=contact_polynomial,
        residual_polynomial=residual_polynomial,
        M=M,
        phi=phi,
        scale_denominator=scale_denominator,
        normalization=normalization,
        scale=scale,
        slope=slope,
        intercept=intercept,
        primitive=primitive,
        distinct_root_factor=distinct_root_factor,
        residual_nonvanishing_factor=residual_nonvanishing_factor,
        weighted_admissibility_factor=weighted_admissibility_factor,
        coefficient_parameters=coefficient_parameters,
        normalized_universal_primitive=normalized_universal_primitive,
        coefficient_space_elimination_ideal=coefficient_space_elimination_ideal,
    )


def multiple_omission_incidence(degree, partitions):
    """Return the incidence for one seed with several omitted values.

    Each partition supplies a factorization
    ``H-s_i*W+t_i=a_i*M_i``.  The coefficient equations identify the
    normalized primitives, while ``pairwise_value_differences`` records the
    cleared numerators of ``s_i-s_j`` and ``t_i-t_j``.  Vanishing of both
    numerators means one common omitted value; the extra Rabinowitsch gates in
    ``distinct_omitted_value_equations`` instead force genuinely different
    target values.
    """
    degree = int(degree)
    partitions = tuple(
        tuple(sorted((int(value) for value in partition), reverse=True))
        for partition in partitions
    )
    if len(partitions) < 2:
        raise ValueError("at least two contact partitions are required")
    contacts = tuple(
        contact_partition_incidence(
            degree,
            partition,
            full_contact=True,
            _prefix=f"_multi{degree}_{index}_{'_'.join(map(str, partition))}",
        )
        for index, partition in enumerate(partitions)
    )
    coefficient_parameters = tuple(
        sp.symbols(f"_multi{degree}_h{index}") for index in range(3, degree)
    )
    equations = []
    marked_points = []
    saturation_factors = []
    omission_coordinates = []
    for contact in contacts:
        parameter_substitution = dict(
            zip(contact.coefficient_parameters, coefficient_parameters)
        )
        equations.extend(
            sp.expand(equation.subs(parameter_substitution))
            for equation in contact.coefficient_space_elimination_ideal.equations
        )
        marked_points.extend(contact.quotient_coordinates)
        saturation_factors.append(contact.weighted_admissibility_factor)
        omission_coordinates.append((contact.s, contact.t))

    pairwise_value_differences = []
    distinct_value_gates = []
    for left in range(len(contacts)):
        for right in range(left + 1, len(contacts)):
            first = contacts[left]
            second = contacts[right]
            first_zero = first.M.subs(first.variable, 0)
            second_zero = second.M.subs(second.variable, 0)
            first_prime_zero = sp.diff(first.M, first.variable).subs(first.variable, 0)
            second_prime_zero = sp.diff(second.M, second.variable).subs(second.variable, 0)
            delta_s = sp.expand(
                first_prime_zero * second.scale_denominator
                - second_prime_zero * first.scale_denominator
            )
            delta_t = sp.expand(
                -first_zero * second.scale_denominator
                + second_zero * first.scale_denominator
            )
            pairwise_value_differences.append((delta_s, delta_t))
            distinct_value_gates.append(
                sp.symbols(
                    f"_multi{degree}_value_gate_{left}_{right}_s "
                    f"_multi{degree}_value_gate_{left}_{right}_t"
                )
            )

    return MultipleOmissionIncidence(
        degree=degree,
        partitions=partitions,
        contacts=contacts,
        coefficient_parameters=coefficient_parameters,
        equations=tuple(equations),
        marked_points=tuple(marked_points),
        saturation_factor=sp.factor(sp.prod(saturation_factors)),
        gate=sp.symbols(f"_multi{degree}_gate"),
        omission_coordinates=tuple(omission_coordinates),
        pairwise_value_differences=tuple(pairwise_value_differences),
        distinct_value_gates=tuple(distinct_value_gates),
    )


def two_omission_incidence(degree, left_partition, right_partition):
    """Build ``a*M_left-b*M_right=alpha*W+beta`` uniformly.

    The returned high-coefficient equations are the compact root-space form
    of the shared-seed incidence.  The full ``multiple_omission_incidence`` is
    retained alongside it to impose endpoint normalization and distinguish
    common from genuinely different omitted target values.
    """
    degree = int(degree)
    shared = multiple_omission_incidence(
        degree, (tuple(left_partition), tuple(right_partition))
    )
    left, right = shared.contacts
    W = sp.symbols(f"_two{degree}_W")
    left_polynomial = left.M.subs(left.variable, W)
    right_polynomial = right.M.subs(right.variable, W)
    left_scale, right_scale = sp.symbols(f"_two{degree}_a _two{degree}_b")
    affine_slope, affine_intercept = sp.symbols(
        f"_two{degree}_alpha _two{degree}_beta"
    )
    affine_difference = sp.expand(
        left_scale * left_polynomial
        - right_scale * right_polynomial
        - affine_slope * W
        - affine_intercept
    )
    polynomial = sp.Poly(affine_difference, W)
    equations = tuple(
        sp.expand(polynomial.coeff_monomial(W**power))
        for power in range(degree + 1)
    )
    collision_discriminant = sp.factor(
        left.distinct_root_factor * right.distinct_root_factor
    )
    return TwoOmissionIncidence(
        degree=degree,
        partitions=(left.multiplicities, right.multiplicities),
        variable=W,
        left=left,
        right=right,
        left_scale=left_scale,
        right_scale=right_scale,
        affine_slope=affine_slope,
        affine_intercept=affine_intercept,
        affine_difference=affine_difference,
        equations=equations,
        collision_discriminant=collision_discriminant,
        saturation_factor=sp.factor(
            left_scale
            * right_scale
            * collision_discriminant
            * left.weighted_admissibility_factor
            * right.weighted_admissibility_factor
        ),
        shared_seed_incidence=shared,
    )


def tangent_chord_normalization(G, W, alpha, beta):
    """Normalize a tangent chord of G to the weighted endpoints zero and one.

    If the tangent at ``alpha`` also meets the graph at ``beta``, the result H
    satisfies H(0)=H'(0)=H(1)=0.  Its dual parameterization differs from that
    of G only by an affine source reparameterization and an invertible affine
    target transformation.
    """
    difference = sp.sympify(beta) - sp.sympify(alpha)
    if difference == 0:
        raise ValueError("the tangent-chord endpoints must be distinct")
    shifted = sp.sympify(alpha) + difference * W
    tangent = G.subs(W, alpha) + sp.diff(G, W).subs(W, alpha) * difference * W
    return sp.expand(G.subs(W, shifted) - tangent)


def universal_primitive(degree: int, W, prefix="h", admissible=True):
    """Return the universal degree-n primitive modulo affine-linear terms.

    In the admissible chart, ``H(0)=H'(0)=H(1)=0`` is imposed identically:
    ``H=sum(h_k*(W**k-W**2), k=3,...,n)``.  The open factor excludes degree
    drop, ``H'(1)=0``, and the forbidden weighted value
    ``H''(1)/(-H'(1))=-2``.
    """
    if degree < 3:
        raise ValueError("inverse degree must be at least three")
    start = 3 if admissible else 2
    coefficients = tuple(
        sp.symbols(f"{prefix}{index}") for index in range(start, degree + 1)
    )
    if admissible:
        polynomial = sum(
            coefficient * (W**index - W**2)
            for index, coefficient in zip(range(3, degree + 1), coefficients)
        )
    else:
        polynomial = sum(
            coefficient * W**index
            for index, coefficient in zip(range(2, degree + 1), coefficients)
        )
    polynomial = sp.expand(polynomial)
    leading = coefficients[-1]
    if admissible:
        first_at_one = sp.diff(polynomial, W).subs(W, 1)
        second_at_one = sp.diff(polynomial, W, 2).subs(W, 1)
        open_factor = leading * first_at_one * (second_at_one - 2 * first_at_one)
    else:
        open_factor = leading
    return UniversalPrimitive(
        degree=degree,
        variable=W,
        polynomial=polynomial,
        coefficients=coefficients,
        open_factor=sp.factor(open_factor),
        admissible=bool(admissible),
    )


def universal_discriminant_incidences(model: UniversalPrimitive, prefix="inc"):
    """Build compatibility wrappers for the original named incidences.

    New full- and partial-contact work should use
    :func:`contact_partition_incidence`.  This dictionary retains the ordinary
    ordered-bitangent incidence and the first three named bad strata used by
    the all-degree theorem. Saturation is represented by a Rabinowitsch
    equation ``1-gate*factor``.
    """
    H = model.polynomial
    W = model.variable
    r, u, v = sp.symbols(f"{prefix}_r {prefix}_u {prefix}_v")
    derivative2 = sp.diff(H, W, 2)
    derivative3 = sp.diff(H, W, 3)

    def at(expression, point):
        return sp.expand(expression.subs(W, point))

    bitangent_ru = bitangent_equations(H, W, r, u)
    bitangent_rv = bitangent_equations(H, W, r, v)

    def incidence(name, equations, points, saturation, contacts):
        gate = sp.symbols(f"{prefix}_{name}_gate")
        return SaturatedIncidence(
            name=name,
            degree=model.degree,
            equations=tuple(sp.expand(equation) for equation in equations),
            marked_points=points,
            parameters=model.coefficients,
            saturation_factor=model.open_factor * saturation,
            gate=gate,
            contact_partition=contacts,
        )

    return {
        "ordinary_bitangent": incidence(
            "ordinary_bitangent",
            bitangent_ru,
            (r, u),
            (r - u) * at(derivative2, r) * at(derivative2, u),
            (2, 2),
        ),
        "higher_cusp": incidence(
            "higher_cusp",
            (at(derivative2, r), at(derivative3, r)),
            (r,),
            sp.Integer(1),
            (4,),
        ),
        "cusp_branch": incidence(
            "cusp_branch",
            (at(derivative2, r),) + bitangent_ru,
            (r, u),
            (r - u) * at(derivative3, r) * at(derivative2, u),
            (3, 2),
        ),
        "tritangent": incidence(
            "tritangent",
            bitangent_ru + bitangent_rv,
            (r, u, v),
            (r - u)
            * (r - v)
            * (u - v)
            * at(derivative2, r)
            * at(derivative2, u)
            * at(derivative2, v),
            (2, 2, 2),
        ),
    }


def incidence_elimination_generators(incidence: SaturatedIncidence):
    """Eliminate marked points and the saturation gate from an incidence."""
    variables = incidence.elimination_variables + incidence.parameters
    basis = sp.groebner(incidence.saturated_equations, *variables, order="lex")
    parameter_set = set(incidence.parameters)
    return tuple(
        sp.factor(polynomial.as_expr())
        for polynomial in basis.polys
        if polynomial.as_expr().free_symbols <= parameter_set
    )


def deterministic_generic_primitive(degree: int, W):
    """Return the rational admissible audit seed used in degrees 3 and above."""
    if degree < 3:
        raise ValueError("inverse degree must be at least three")
    extra = sum((index + 1) * W**index for index in range(degree - 2))
    return sp.expand(W**2 * (1 - W) * extra)


def partition_dual_geometry(partition):
    """Annotate a full root-multiplicity partition by its dual-curve geometry."""
    multiplicities = tuple(sorted((int(value) for value in partition), reverse=True))
    nontrivial = tuple(value for value in multiplicities if value > 1)
    if len(nontrivial) == 1 and nontrivial[0] == sum(multiplicities):
        if nontrivial[0] == 3:
            return "ordinary cusp; maximally ramified"
        return "maximally ramified higher cusp"
    if nontrivial == (3, 2):
        return "ordinary cusp branch meeting another tangent branch"
    if nontrivial == (2, 2, 2):
        return "tritangent line / triple normalization point"
    if nontrivial == (3,):
        return "ordinary cusp"
    if nontrivial == (2, 2):
        return "ordinary node"
    if nontrivial == (4,):
        return "higher cusp"
    return "higher or multiple dual-curve singularity"
