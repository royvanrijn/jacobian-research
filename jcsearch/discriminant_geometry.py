"""Dual-curve geometry of the discriminant normalization."""

from __future__ import annotations

from dataclasses import dataclass
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


def contact_partition_incidence(degree, multiplicities, full_contact=True):
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
    prefix = f"_cp{degree}_{partition_tag}"
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
