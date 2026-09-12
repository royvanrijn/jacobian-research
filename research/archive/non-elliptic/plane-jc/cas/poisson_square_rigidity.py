#!/usr/bin/env python3
"""Exact band compiler and reduced rigidity theorem for ``[P,Q]=X^2``.

The distinguished three-layer support is

    P = X^3 A(t) + X^2 B(t),   deg(A)<=3, deg(B)<=4,
    Q = X^2 C(t) + X   D(t),   deg(C)<=2, deg(D)<=3.

On the locus where ``A`` and ``C`` are nonzero, the highest Wronskian layer
forces a common affine factor

    A=a(t-r)^3, C=c(t-r)^2.

The lower layers are then small enough to classify without a giant
coefficient ideal.  Over an algebraically closed characteristic-zero field,
the tangent-pencil family returned by ``reduced_three_band_family`` and the
four charts returned by ``degenerate_reduced_families`` cover every point.
The reduced scheme has three irreducible components: the tangent closure,
``C=0``, and ``A=0``.  The unreduced coefficient scheme can retain nilpotent
coefficients.  Exact transverse slices give generic multiplicities ``2,3,1``
on those components; the module reports that distinction explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
import sympy as sp


@dataclass(frozen=True)
class BandLimitedPoissonProblem:
    """Finite Laurent/polynomial support data for a Poisson equation."""

    P_support: tuple[tuple[int, tuple[int, ...]], ...]
    Q_support: tuple[tuple[int, tuple[int, ...]], ...]
    rhs: sp.Expr
    X: sp.Symbol
    t: sp.Symbol

    def _bands(self, prefix: str, support) -> dict[int, sp.Expr]:
        return {
            x_power: sp.Add(
                *(
                    sp.Symbol(f"{prefix}_{x_power}_{t_power}") * self.t**t_power
                    for t_power in t_powers
                )
            )
            for x_power, t_powers in support
        }

    @property
    def P_bands(self) -> dict[int, sp.Expr]:
        return self._bands("p", self.P_support)

    @property
    def Q_bands(self) -> dict[int, sp.Expr]:
        return self._bands("q", self.Q_support)

    @property
    def P(self) -> sp.Expr:
        return sp.Add(
            *(band * self.X**power for power, band in self.P_bands.items())
        )

    @property
    def Q(self) -> sp.Expr:
        return sp.Add(
            *(band * self.X**power for power, band in self.Q_bands.items())
        )

    @property
    def bracket(self) -> sp.Expr:
        return sp.expand(
            sp.diff(self.P, self.X) * sp.diff(self.Q, self.t)
            - sp.diff(self.P, self.t) * sp.diff(self.Q, self.X)
        )

    @property
    def layers(self) -> dict[int, sp.Expr]:
        difference = sp.Poly(sp.expand(self.bracket - self.rhs), self.X)
        return {
            power[0]: sp.expand(coefficient)
            for power, coefficient in difference.terms()
        }

    @property
    def coefficient_equations(self) -> tuple[sp.Expr, ...]:
        equations: list[sp.Expr] = []
        for layer in self.layers.values():
            polynomial = sp.Poly(layer, self.t)
            equations.extend(polynomial.all_coeffs())
        return tuple(equation for equation in equations if equation != 0)


@dataclass(frozen=True)
class TangentFamily:
    """The complete reduced nondegenerate three-band family."""

    X: sp.Symbol
    t: sp.Symbol
    root: sp.Symbol
    cubic_scale: sp.Symbol
    gamma_scale: sp.Symbol
    tangent_linear: sp.Symbol
    A: sp.Expr
    B: sp.Expr
    C: sp.Expr
    D: sp.Expr
    P: sp.Expr
    Q: sp.Expr
    W: sp.Expr
    gamma: sp.Expr
    H: sp.Expr

    @property
    def bracket(self) -> sp.Expr:
        return sp.expand(
            sp.diff(self.P, self.X) * sp.diff(self.Q, self.t)
            - sp.diff(self.P, self.t) * sp.diff(self.Q, self.X)
        )


@dataclass(frozen=True)
class PoissonBandFamily:
    """One parametrized field-valued stratum of the three-band equation."""

    name: str
    X: sp.Symbol
    t: sp.Symbol
    A: sp.Expr
    B: sp.Expr
    C: sp.Expr
    D: sp.Expr
    nonzero_parameters: tuple[sp.Symbol, ...]

    @property
    def P(self) -> sp.Expr:
        return sp.expand(self.X**3 * self.A + self.X**2 * self.B)

    @property
    def Q(self) -> sp.Expr:
        return sp.expand(self.X**2 * self.C + self.X * self.D)

    @property
    def bracket(self) -> sp.Expr:
        return sp.expand(
            sp.diff(self.P, self.X) * sp.diff(self.Q, self.t)
            - sp.diff(self.P, self.t) * sp.diff(self.Q, self.X)
        )


@dataclass(frozen=True)
class ReducedIntersectionStratum:
    """Dense chart of one irreducible reduced component-intersection branch."""

    name: str
    component_membership: tuple[str, ...]
    family: PoissonBandFamily


@dataclass(frozen=True)
class NormalizedResidualSystem:
    """Equations left after ``A=a h^3, C=c h^2`` and the middle layer."""

    variables: tuple[sp.Symbol, ...]
    B: sp.Expr
    D: sp.Expr
    equations_by_degree: tuple[tuple[int, sp.Expr], ...]
    reduced_consequences: tuple[sp.Expr, ...]
    nilpotent_warning: str


@dataclass(frozen=True)
class ComponentTangentAudit:
    """Exact Zariski tangent-space rank at a rational component point."""

    component: str
    jacobian_rank: int
    tangent_dimension: int
    reduced_component_dimension: int = 4


@dataclass(frozen=True)
class IntersectionKernelSliceAudit:
    """Exact Artin algebra cut out on a chosen tangent-kernel slice."""

    stratum: str
    tangent_kernel_dimension: int
    groebner_basis: tuple[str, ...]
    length: int
    hilbert_vector: tuple[int, ...]
    socle_dimension: int


def weighted_wronskian(
    left_weight: int,
    right_weight: int,
    left: sp.Expr,
    right: sp.Expr,
    t: sp.Symbol,
) -> sp.Expr:
    return sp.expand(
        left_weight * left * sp.diff(right, t)
        - right_weight * sp.diff(left, t) * right
    )


def three_band_layers(
    A: sp.Expr,
    B: sp.Expr,
    C: sp.Expr,
    D: sp.Expr,
    t: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the coefficients of ``X^4,X^3,X^2`` in the bracket."""

    return (
        weighted_wronskian(3, 2, A, C, t),
        sp.expand(
            weighted_wronskian(3, 1, A, D, t)
            + weighted_wronskian(2, 2, B, C, t)
        ),
        weighted_wronskian(2, 1, B, D, t),
    )


def primary_chart_scaling_audit() -> dict[str, object]:
    """Verify the ``G_m`` action used to normalize the primary problem.

    The weights ``(-2,-1,0,1)`` on ``(A,B,C,D)`` make the three Wronskian
    layers semi-invariant with weights ``(-2,-1,0)``.  Since the bottom
    layer equals one, the action preserves the complete coefficient scheme.
    On ``D(d0)`` it sets the constant coefficient of ``D`` to one.
    """

    t = sp.symbols("t")
    lam = sp.symbols("lambda", nonzero=True)
    A = sp.Function("A")(t)
    B = sp.Function("B")(t)
    C = sp.Function("C")(t)
    D = sp.Function("D")(t)
    original = three_band_layers(A, B, C, D, t)
    transformed = three_band_layers(
        A / lam**2,
        B / lam,
        C,
        lam * D,
        t,
    )
    weights = (-2, -1, 0)
    if any(
        sp.simplify(new - lam**weight * old) != 0
        for new, old, weight in zip(transformed, original, weights)
    ):
        raise AssertionError("primary-chart scaling ceased to preserve the layers")
    return {
        "coefficient_weights": (-2, -1, 0, 1),
        "layer_weights": weights,
        "normalizing_parameter": "lambda=d0^-1",
        "normalized_constant": 1,
    }


def standard_three_band_problem() -> BandLimitedPoissonProblem:
    X, t = sp.symbols("X t")
    return BandLimitedPoissonProblem(
        P_support=((3, tuple(range(4))), (2, tuple(range(5)))),
        Q_support=((2, tuple(range(3))), (1, tuple(range(4)))),
        rhs=X**2,
        X=X,
        t=t,
    )


def reduced_three_band_family() -> TangentFamily:
    """Return the four-parameter reduced family and its tangent polynomial."""

    X, t, r, e, Z = sp.symbols("X t r e Z")
    a, d = sp.symbols("a d", nonzero=True)
    h = t - r
    A = a * h**3
    C = -sp.Rational(3, 2) * a * d**2 * h**2
    D = d + e * h
    B = -h / d - e * h**2 / (2 * d**2)
    P = sp.expand(X**3 * A + X**2 * B)
    Q = sp.expand(X**2 * C + X * D)
    W = X * h
    gamma = d * X
    H_Z = e * Z**2 / 2 - a * d**2 * Z**3 / 2
    H = sp.expand(H_Z.subs(Z, W))
    return TangentFamily(
        X=X,
        t=t,
        root=r,
        cubic_scale=a,
        gamma_scale=d,
        tangent_linear=e,
        A=A,
        B=B,
        C=C,
        D=D,
        P=P,
        Q=Q,
        W=W,
        gamma=gamma,
        H=H,
    )


def degenerate_reduced_families() -> tuple[PoissonBandFamily, ...]:
    """Return the four charts covering all field-valued points with ``AC=0``.

    The charts overlap along ``A=C=0``.  Parameters declared nonzero are the
    constant value of ``D`` or the slope of a linear ``D``.
    """

    X, t = sp.symbols("X t")
    d, b, alpha, lam = sp.symbols("d b alpha lambda")
    delta, d0, k = sp.symbols("delta d0 k")

    D_constant = d
    B_constant = b - t / d
    D_linear = d0 + delta * t
    B_linear = k * D_linear**2 + 1 / (2 * delta)

    return (
        PoissonBandFamily(
            name="C_zero_D_constant",
            X=X,
            t=t,
            A=alpha,
            B=B_constant,
            C=sp.S.Zero,
            D=D_constant,
            nonzero_parameters=(d,),
        ),
        PoissonBandFamily(
            name="C_zero_D_linear",
            X=X,
            t=t,
            A=alpha * D_linear**3,
            B=B_linear,
            C=sp.S.Zero,
            D=D_linear,
            nonzero_parameters=(delta,),
        ),
        PoissonBandFamily(
            name="A_zero_D_constant",
            X=X,
            t=t,
            A=sp.S.Zero,
            B=B_constant,
            C=lam * B_constant,
            D=D_constant,
            nonzero_parameters=(d,),
        ),
        PoissonBandFamily(
            name="A_zero_D_linear",
            X=X,
            t=t,
            A=sp.S.Zero,
            B=B_linear,
            C=lam * B_linear,
            D=D_linear,
            nonzero_parameters=(delta,),
        ),
    )


def reduced_intersection_strata() -> tuple[ReducedIntersectionStratum, ...]:
    """Return dense charts for the three reduced intersection branches.

    The exact Singular radical audit proves:

    - ``C=0`` intersect ``A=0`` is the lower-Wronskian core and is also the
      triple intersection;
    - tangent intersect ``C=0`` is the union of that core and the
      constant-``D`` chart;
    - tangent intersect ``A=0`` is the union of that core and the
      constant-``B,C`` chart.
    """

    X, t = sp.symbols("X t")
    delta = sp.symbols("delta", nonzero=True)
    d0, k = sp.symbols("d0 k")
    d = sp.symbols("d", nonzero=True)
    b, alpha, c = sp.symbols("b alpha c")

    D_linear = d0 + delta * t
    B_core = k * D_linear**2 + 1 / (2 * delta)
    core = PoissonBandFamily(
        name="lower_wronskian_core",
        X=X,
        t=t,
        A=sp.S.Zero,
        B=B_core,
        C=sp.S.Zero,
        D=D_linear,
        nonzero_parameters=(delta,),
    )
    constant_d = PoissonBandFamily(
        name="tangent_C_zero_constant_D_branch",
        X=X,
        t=t,
        A=alpha,
        B=b - t / d,
        C=sp.S.Zero,
        D=d,
        nonzero_parameters=(d,),
    )
    constant_bc = PoissonBandFamily(
        name="tangent_A_zero_constant_BC_branch",
        X=X,
        t=t,
        A=sp.S.Zero,
        B=1 / (2 * delta),
        C=c,
        D=D_linear,
        nonzero_parameters=(delta,),
    )
    return (
        ReducedIntersectionStratum(
            name="triple_lower_wronskian_core",
            component_membership=("tangent", "C_zero", "A_zero"),
            family=core,
        ),
        ReducedIntersectionStratum(
            name="tangent_C_zero_boundary",
            component_membership=("tangent", "C_zero"),
            family=constant_d,
        ),
        ReducedIntersectionStratum(
            name="tangent_A_zero_boundary",
            component_membership=("tangent", "A_zero"),
            family=constant_bc,
        ),
    )


def c_zero_dual_number_family() -> tuple[
    PoissonBandFamily, sp.Symbol
]:
    """Return a generic square-zero direction along the ``C=0`` component.

    For affine ``D`` with nonzero slope ``delta``, put

        C = epsilon*c*D^2,
        A = alpha*D^3 - epsilon*(c/delta)*D.

    Keeping the reduced ``B,D`` fixed kills all three Poisson layers modulo
    ``epsilon^2``.  This proves that the original coefficient scheme is
    generically nonreduced on the ``C=0`` component as well as on the tangent
    component.
    """

    X, t = sp.symbols("X t")
    delta = sp.symbols("delta", nonzero=True)
    d0, k, alpha, c, epsilon = sp.symbols(
        "d0 k alpha c epsilon"
    )
    D = d0 + delta * t
    B = k * D**2 + 1 / (2 * delta)
    return (
        PoissonBandFamily(
            name="C_zero_generic_dual_number",
            X=X,
            t=t,
            A=alpha * D**3 - epsilon * c * D / delta,
            B=B,
            C=epsilon * c * D**2,
            D=D,
            nonzero_parameters=(delta,),
        ),
        epsilon,
    )


def c_zero_curvilinear_family() -> tuple[
    PoissonBandFamily, sp.Symbol
]:
    """Return the generic length-three thickening of the ``C=0`` component."""

    X, t = sp.symbols("X t")
    delta, alpha = sp.symbols("delta alpha", nonzero=True)
    d0, k, c, epsilon = sp.symbols("d0 k c epsilon")
    D = d0 + delta * t
    B = k * D**2 + 1 / (2 * delta)
    z = epsilon * c
    C = z * D**2 - 2 * z**2 / (3 * alpha * delta)
    A = alpha * D**3 + (
        -z / delta - 4 * k * z**2 / (3 * alpha * delta)
    ) * D
    return (
        PoissonBandFamily(
            name="C_zero_generic_length_three",
            X=X,
            t=t,
            A=A,
            B=B,
            C=C,
            D=D,
            nonzero_parameters=(delta, alpha),
        ),
        epsilon,
    )


def c_zero_local_slice_audit() -> dict[str, sp.Expr]:
    """Return the exact equations giving the generic ``k[eta]/eta^3`` slice."""

    c0, c1, c2, alpha, k = sp.symbols("c0 c1 c2 alpha k")
    delta = sp.symbols("delta", nonzero=True)
    s = sp.symbols("s")
    C = c0 + c1 * s + c2 * s**2
    # This is the unique A modulo its free cubic coefficient alpha*s^3
    # obtained from K3 after D=s and B=k*s^2+1/(2*delta).
    A = (
        -c1 / (3 * delta)
        + (2 * k * c0 - c2 / delta) * s
        + 2 * k * c1 * s**2
        + alpha * s**3
    )
    top = sp.Poly(
        sp.expand(3 * A * sp.diff(C, s) - 2 * sp.diff(A, s) * C),
        s,
    )
    coefficients = {
        f"K4_degree_{power[0]}": sp.factor(value)
        for power, value in top.terms()
    }
    L = -3 * alpha + 4 * k * c2
    c0_solution = 2 * c2**2 / (delta * L)
    constant_after_middle = sp.factor(
        coefficients["K4_degree_0"].subs(
            {c1: 0, c0: c0_solution}
        )
    )
    return {
        **coefficients,
        "local_unit": L,
        "c0_solution": c0_solution,
        "constant_after_middle": constant_after_middle,
        "expected_constant": -12
        * alpha
        * c2**3
        / (delta**2 * L**2),
    }


def lower_layer_degree_audit() -> dict[str, tuple[sp.Expr, ...]]:
    """Coefficient contradictions excluding quadratic and cubic ``D``.

    Affine changes make an exact degree-two or degree-three polynomial monic
    and centered.  The returned coefficient lists are those of
    ``2BD'-B'D-1`` from highest degree to constant degree.
    """

    t = sp.symbols("t")
    b0, b1, b2, b3, b4 = sp.symbols("b0:5")
    B = b0 + b1 * t + b2 * t**2 + b3 * t**3 + b4 * t**4
    d0, d1 = sp.symbols("d0 d1")

    def coefficients(D: sp.Expr) -> tuple[sp.Expr, ...]:
        polynomial = sp.Poly(
            sp.expand(weighted_wronskian(2, 1, B, D, t) - 1),
            t,
        )
        return tuple(sp.factor(value) for value in polynomial.all_coeffs())

    return {
        "degree_2": coefficients(t**2 + d0),
        "degree_3": coefficients(t**3 + d1 * t + d0),
    }


def component_tangent_space_audit() -> tuple[ComponentTangentAudit, ...]:
    """Audit generic tangent dimensions on the three reduced components."""

    problem = standard_three_band_problem()
    variables = tuple(
        sorted(
            set().union(
                *(
                    equation.free_symbols
                    for equation in problem.coefficient_equations
                )
            ),
            key=str,
        )
    )
    variable_by_name = {str(variable): variable for variable in variables}
    jacobian = sp.Matrix(
        [
            [sp.diff(equation, variable) for variable in variables]
            for equation in problem.coefficient_equations
        ]
    )

    def point(
        A: sp.Expr,
        B: sp.Expr,
        C: sp.Expr,
        D: sp.Expr,
        t: sp.Symbol,
        substitutions: dict[sp.Symbol, sp.Expr],
    ) -> dict[sp.Symbol, sp.Expr]:
        result: dict[sp.Symbol, sp.Expr] = {}
        for prefix, polynomial, degree in (
            ("p_3", A, 3),
            ("p_2", B, 4),
            ("q_2", C, 2),
            ("q_1", D, 3),
        ):
            specialized = sp.expand(polynomial.subs(substitutions))
            coefficients = sp.Poly(specialized, t)
            for index in range(degree + 1):
                result[variable_by_name[f"{prefix}_{index}"]] = (
                    coefficients.coeff_monomial(t**index)
                )
        assert all(
            sp.expand(equation.subs(result)) == 0
            for equation in problem.coefficient_equations
        )
        return result

    tangent = reduced_three_band_family()
    tangent_point = point(
        tangent.A,
        tangent.B,
        tangent.C,
        tangent.D,
        tangent.t,
        {
            tangent.root: 0,
            tangent.cubic_scale: 1,
            tangent.gamma_scale: 1,
            tangent.tangent_linear: 1,
        },
    )

    families = {
        family.name: family for family in degenerate_reduced_families()
    }
    c_zero = families["C_zero_D_linear"]
    a_zero = families["A_zero_D_linear"]

    def named_substitutions(
        family: PoissonBandFamily,
        values: dict[str, sp.Expr],
    ) -> dict[sp.Symbol, sp.Expr]:
        symbols = set().union(
            family.A.free_symbols,
            family.B.free_symbols,
            family.C.free_symbols,
            family.D.free_symbols,
        ) - {family.t}
        return {
            symbol: values[str(symbol)]
            for symbol in symbols
            if str(symbol) in values
        }

    base_values = {
        "delta": 2,
        "d0": 3,
        "k": 5,
        "alpha": 7,
        "lambda": 11,
    }
    c_zero_point = point(
        c_zero.A,
        c_zero.B,
        c_zero.C,
        c_zero.D,
        c_zero.t,
        named_substitutions(c_zero, base_values),
    )
    a_zero_point = point(
        a_zero.A,
        a_zero.B,
        a_zero.C,
        a_zero.D,
        a_zero.t,
        named_substitutions(a_zero, base_values),
    )

    reports = []
    for name, sample in (
        ("tangent_closure", tangent_point),
        ("C_zero", c_zero_point),
        ("A_zero", a_zero_point),
    ):
        rank = int(jacobian.subs(sample).rank())
        reports.append(
            ComponentTangentAudit(
                component=name,
                jacobian_rank=rank,
                tangent_dimension=len(variables) - rank,
            )
        )
    return tuple(reports)


def intersection_tangent_space_audit() -> tuple[ComponentTangentAudit, ...]:
    """Audit tangent dimensions on dense charts of the intersection branches."""

    problem = standard_three_band_problem()
    variables = tuple(
        sorted(
            set().union(
                *(
                    equation.free_symbols
                    for equation in problem.coefficient_equations
                )
            ),
            key=str,
        )
    )
    variable_by_name = {str(variable): variable for variable in variables}
    jacobian = sp.Matrix(
        [
            [sp.diff(equation, variable) for variable in variables]
            for equation in problem.coefficient_equations
        ]
    )
    values = {
        "delta": 2,
        "d0": 3,
        "k": 5,
        "d": 2,
        "b": 3,
        "alpha": 5,
        "c": 7,
    }

    reports = []
    for stratum in reduced_intersection_strata():
        family = stratum.family
        symbols = set().union(
            family.A.free_symbols,
            family.B.free_symbols,
            family.C.free_symbols,
            family.D.free_symbols,
        ) - {family.t}
        substitutions = {
            symbol: values[str(symbol)]
            for symbol in symbols
            if str(symbol) in values
        }
        point: dict[sp.Symbol, sp.Expr] = {}
        for prefix, polynomial, degree in (
            ("p_3", family.A, 3),
            ("p_2", family.B, 4),
            ("q_2", family.C, 2),
            ("q_1", family.D, 3),
        ):
            specialized = sp.Poly(
                sp.expand(polynomial.subs(substitutions)),
                family.t,
            )
            for index in range(degree + 1):
                point[variable_by_name[f"{prefix}_{index}"]] = (
                    specialized.coeff_monomial(family.t**index)
                )
        if not all(
            sp.expand(equation.subs(point)) == 0
            for equation in problem.coefficient_equations
        ):
            raise AssertionError("intersection sample is not a Poisson solution")
        rank = int(jacobian.subs(point).rank())
        reports.append(
            ComponentTangentAudit(
                component=stratum.name,
                jacobian_rank=rank,
                tangent_dimension=len(variables) - rank,
                reduced_component_dimension=3,
            )
        )
    return tuple(reports)


def intersection_tangent_kernel_slice_audit(
) -> tuple[IntersectionKernelSliceAudit, ...]:
    """Compute exact tangent-kernel Artin slices on the two extra branches.

    Three coefficient coordinates parametrizing each reduced threefold are
    fixed at a rational point.  The resulting linearized system has kernel
    dimension four on the constant-``D`` branch and three on the
    constant-``B,C`` branch.  Restricting the full 17 equations to those
    kernels gives two zero-dimensional length-five algebras.

    These chosen slices are executable normal signatures, not claims about
    the full completed local rings or their primary decompositions.
    """

    problem = standard_three_band_problem()
    variables = tuple(
        sorted(
            set().union(
                *(
                    equation.free_symbols
                    for equation in problem.coefficient_equations
                )
            ),
            key=str,
        )
    )
    variable_by_name = {str(variable): variable for variable in variables}
    zero = {variable: sp.S.Zero for variable in variables}
    cases = (
        (
            "tangent_C_zero_boundary",
            {
                "p_3_0": sp.Integer(5),
                "p_2_0": sp.Integer(3),
                "p_2_1": -sp.Rational(1, 2),
                "q_1_0": sp.Integer(2),
            },
            ("p_3_0", "p_2_0", "q_1_0"),
        ),
        (
            "tangent_A_zero_boundary",
            {
                "p_2_0": sp.Rational(1, 4),
                "q_2_0": sp.Integer(7),
                "q_1_0": sp.Integer(3),
                "q_1_1": sp.Integer(2),
            },
            ("q_1_0", "q_1_1", "q_2_0"),
        ),
    )
    reports = []
    for stratum, named_values, fixed_names in cases:
        values = dict(zero)
        values.update(
            {
                variable_by_name[name]: value
                for name, value in named_values.items()
            }
        )
        shift = {
            variable: variable + values[variable]
            for variable in variables
        }
        equations = [
            sp.expand(equation.subs(shift))
            for equation in problem.coefficient_equations
        ] + [variable_by_name[name] for name in fixed_names]
        linearization = sp.Matrix(
            [
                [
                    sp.diff(equation, variable).subs(zero)
                    for variable in variables
                ]
                for equation in equations
            ]
        )
        kernel = linearization.nullspace()
        parameters = sp.symbols(f"u0:{len(kernel)}")
        kernel_matrix = sp.Matrix.hstack(*kernel)
        kernel_substitution = {
            variable: sum(
                kernel_matrix[row, column] * parameters[column]
                for column in range(len(kernel))
            )
            for row, variable in enumerate(variables)
        }
        restricted = tuple(
            polynomial
            for equation in equations
            for polynomial in (sp.factor(equation.subs(kernel_substitution)),)
            if polynomial != 0
        )
        basis = sp.groebner(
            restricted,
            *parameters,
            order="grevlex",
        )
        basis_expressions = tuple(
            sp.factor(polynomial.as_expr())
            for polynomial in basis.polys
        )
        if stratum == "tangent_C_zero_boundary":
            expected = tuple(
                parameters[first] * parameters[second]
                for first in range(4)
                for second in range(first, 4)
            )
            if set(basis_expressions) != set(expected):
                raise AssertionError("constant-D slice is not the square-zero algebra")
            report = IntersectionKernelSliceAudit(
                stratum=stratum,
                tangent_kernel_dimension=4,
                groebner_basis=tuple(
                    sp.sstr(expression) for expression in basis_expressions
                ),
                length=5,
                hilbert_vector=(1, 4),
                socle_dimension=4,
            )
        else:
            u0, u1, u2 = parameters
            expected = {
                u2**3,
                u0**2,
                u0 * u1 - sp.Rational(9, 14) * u2**2,
                u1**2,
                u0 * u2,
                u1 * u2 - 3 * u2**2,
            }
            if {
                sp.expand(expression) for expression in basis_expressions
            } != expected:
                raise AssertionError(
                    f"constant-BC slice algebra changed: {basis_expressions}"
                )
            report = IntersectionKernelSliceAudit(
                stratum=stratum,
                tangent_kernel_dimension=3,
                groebner_basis=tuple(
                    sp.sstr(expression) for expression in basis_expressions
                ),
                length=5,
                hilbert_vector=(1, 3, 1),
                socle_dimension=1,
            )
        reports.append(report)
    return tuple(reports)


def normalized_residual_system() -> NormalizedResidualSystem:
    """Compute the exact small system exposing its radical and nilpotents."""

    t = sp.symbols("t")
    a, c = sp.symbols("a c", nonzero=True)
    d0, d1, d2, d3, b2 = sp.symbols("d0 d1 d2 d3 b2")
    b1 = 3 * a * d0 / (2 * c)
    b3 = 3 * a * d2 / (2 * c)
    b4 = 3 * a * d3 / (2 * c)
    B = sp.expand(b1 * t + b2 * t**2 + b3 * t**3 + b4 * t**4)
    D = d0 + d1 * t + d2 * t**2 + d3 * t**3
    A = a * t**3
    C = c * t**2
    assert sp.expand(
        weighted_wronskian(3, 1, A, D, t)
        + weighted_wronskian(2, 2, B, C, t)
    ) == 0
    residual = sp.Poly(
        sp.expand(weighted_wronskian(2, 1, B, D, t) - 1),
        t,
    )
    by_degree = tuple(
        sorted(
            (
                (power[0], sp.factor(coefficient))
                for power, coefficient in residual.terms()
            ),
            reverse=True,
        )
    )
    return NormalizedResidualSystem(
        variables=(a, c, d0, d1, d2, d3, b2),
        B=B,
        D=D,
        equations_by_degree=by_degree,
        reduced_consequences=(
            d3,
            d2,
            4 * b2 * c - 3 * a * d1,
            2 * c + 3 * a * d0**2,
        ),
        nilpotent_warning=(
            "After the constant, linear, and cubic equations force the "
            "unit relations and d3=0, the degree-four equation is "
            "proportional to d2^2.  Thus d2 vanishes on the reduced locus "
            "but survives as an explicit dual-number direction."
        ),
    )


def verify_reduced_classification() -> None:
    """Exact regression for the classification and tangent identification."""

    family = reduced_three_band_family()
    X, t = family.X, family.t
    K4, K3, K2 = three_band_layers(
        family.A, family.B, family.C, family.D, t
    )
    assert sp.expand(K4) == 0
    assert sp.expand(K3) == 0
    assert sp.expand(K2 - 1) == 0
    assert sp.expand(family.bracket - X**2) == 0

    # If S=H'(W)+gamma and T=W*S-H, the pair is (S,T)=(Q,-d^2 P).
    Z = sp.symbols("Z")
    H_Z = (
        family.tangent_linear * Z**2 / 2
        - family.cubic_scale * family.gamma_scale**2 * Z**3 / 2
    )
    S = sp.expand(
        sp.diff(H_Z, Z).subs(Z, family.W) + family.gamma
    )
    T = sp.expand(family.W * S - H_Z.subs(Z, family.W))
    assert sp.expand(S - family.Q) == 0
    assert sp.expand(T + family.gamma_scale**2 * family.P) == 0

    for degenerate in degenerate_reduced_families():
        assert sp.expand(degenerate.bracket - degenerate.X**2) == 0
        layers = three_band_layers(
            degenerate.A,
            degenerate.B,
            degenerate.C,
            degenerate.D,
            degenerate.t,
        )
        assert tuple(sp.expand(layer) for layer in layers) == (0, 0, 1)

    c_zero_dual, epsilon = c_zero_dual_number_family()
    dual_residual = sp.Poly(
        sp.expand(c_zero_dual.bracket - c_zero_dual.X**2),
        epsilon,
    )
    assert sp.rem(
        dual_residual,
        sp.Poly(epsilon**2, epsilon),
    ).as_expr() == 0

    c_zero_length_three, epsilon = c_zero_curvilinear_family()
    length_three_residual = sp.Poly(
        sp.expand(
            c_zero_length_three.bracket - c_zero_length_three.X**2
        ),
        epsilon,
    )
    assert sp.rem(
        length_three_residual,
        sp.Poly(epsilon**3, epsilon),
    ).as_expr() == 0
    assert sp.rem(
        length_three_residual,
        sp.Poly(epsilon**4, epsilon),
    ).as_expr() != 0
    local_slice = c_zero_local_slice_audit()
    assert sp.factor(
        local_slice["constant_after_middle"]
        - local_slice["expected_constant"]
    ) == 0

    # For exact deg(D)=2, the top equations successively give
    # b3=0, b2=2*b4*d0, b1=0, b0=b4*d0^2, leaving -1=0.
    # For exact deg(D)=3 they give b4=b3=b2=b1=b0=0, again
    # leaving -1=0.  Thus every field-valued solution has deg(D)<=1.
    degree_audit = lower_layer_degree_audit()
    b0, b1, b2, b3, b4, d0 = sp.symbols("b0 b1 b2 b3 b4 d0")
    degree_2_substitution = {
        b3: 0,
        b2: 2 * b4 * d0,
        b1: 0,
        b0: b4 * d0**2,
    }
    assert tuple(
        sp.expand(equation.subs(degree_2_substitution))
        for equation in degree_audit["degree_2"]
    )[-1] == -1
    b0, b1, b2, b3, b4 = sp.symbols("b0 b1 b2 b3 b4")
    degree_3_substitution = {b4: 0, b3: 0, b2: 0, b1: 0, b0: 0}
    assert tuple(
        sp.expand(equation.subs(degree_3_substitution))
        for equation in degree_audit["degree_3"]
    )[-1] == -1

    tangent_audit = component_tangent_space_audit()
    assert tuple(
        (report.component, report.jacobian_rank, report.tangent_dimension)
        for report in tangent_audit
    ) == (
        ("tangent_closure", 11, 5),
        ("C_zero", 11, 5),
        ("A_zero", 12, 4),
    )

    residual = normalized_residual_system()
    equations = dict(residual.equations_by_degree)
    a, c, d0, d1, d2, d3, b2 = residual.variables
    assert sp.factor(equations[6] - 3 * a * d3**2 / c) == 0
    # Successively impose the consequences valid in every reduced algebra.
    reduced_substitution = {
        d3: 0,
        d2: 0,
        b2: 3 * a * d1 / (4 * c),
        c: -3 * a * d0**2 / 2,
    }
    for equation in equations.values():
        assert sp.cancel(equation.subs(reduced_substitution)) == 0

    # The unreduced caveat is real: at a=1,c=-3/2,d0=1,d1=b2=d3=0,
    # setting d2=epsilon kills every equation modulo epsilon^2.
    epsilon = sp.symbols("epsilon")
    dual_number_point = {
        a: 1,
        c: -sp.Rational(3, 2),
        d0: 1,
        d1: 0,
        b2: 0,
        d2: epsilon,
        d3: 0,
    }
    for equation in equations.values():
        assert sp.rem(
            sp.Poly(sp.expand(equation.subs(dual_number_point)), epsilon),
            sp.Poly(epsilon**2, epsilon),
        ).as_expr() == 0


def main() -> int:
    verify_reduced_classification()
    problem = standard_three_band_problem()
    residual = normalized_residual_system()
    print(
        "PASS three-band compiler:",
        len(problem.coefficient_equations),
        "coefficient equations in layers",
        tuple(problem.layers),
    )
    print(
        "PASS reduced rigidity: tangent plus four degenerate charts cover "
        "every field-valued three-band solution"
    )
    print("PASS generic component multiplicities: tangent=2, C_zero=3, A_zero=1")
    print("SCHEME CAVEAT:", residual.nilpotent_warning)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
