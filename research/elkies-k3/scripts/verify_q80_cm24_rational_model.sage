#!/usr/bin/env sage
"""Verify the rational q=80 CM24 boundary surface and its tangent cone.

The surface was recovered from the unique ordinary lift of the fully marked
mod-7 short-basis seed.  This script checks its fibers exactly and computes
the two four-variable polynomial-section schemes used by the CM24 basis:
D1/D2 on the D5 spinor chart, and D3 on the D5 identity chart with its I4/I2
node incidences built in.

The Weierstrass surface is defined over QQ.  The three selected marked
sections are defined over the compositum of QQ(sqrt(-3)) and QQ(sqrt(-6)).
The final characteristic-zero calculation certifies two rational tangent-cone
lines at this enhanced CM point.  It does *not* construct or algebraize either
one-parameter rank-19 branch, so this file is an exact boundary certificate,
not a generic-family verifier.
"""

from sage.all import *


polynomials = PolynomialRing(QQ, "T")
T = polynomials.gen()
A = T**2*(-3+QQ(9)/4*T-QQ(9)/4*T**2+QQ(9)/4*T**3)
B = T**3*(2-QQ(315)/32*T+9*T**2-QQ(9)/16*T**3-QQ(27)/32*T**5)
discriminant = 4*A**3+27*B**2
factorization = discriminant.factor()
rho = QQ(-3)
node = QQ(-33)/2
assert discriminant.valuation(T) == 7
assert discriminant.valuation(T-1) == 4
assert discriminant.valuation(T-rho) == 2
assert gcd(
    PolynomialRing(QQ, "x").gen()**3+A(rho)*PolynomialRing(QQ, "x").gen()+B(rho),
    3*PolynomialRing(QQ, "x").gen()**2+A(rho),
).monic() == PolynomialRing(QQ, "x").gen()-node
print(f"Q80CM24QQ|A={A}|B={B}|Delta={factorization}", flush=True)


quadratic.<sqrt_minus_six> = QuadraticField(-6)


def section_scheme(profile):
    ring = PolynomialRing(QQ, names=("a", "b", "c", "d"), order="degrevlex")
    a, b, c, d = ring.gens()
    local_polynomials = PolynomialRing(ring, "t")
    t = local_polynomials.gen()
    if profile == "spinor":
        X = t*(1+a*t)
        Y = t**2*(b+c*t+d*t**2)
    else:
        # X passes through (d,w) at T=1,-3; solve its two upper coefficients
        # from the free constant a.  Y vanishes at both nodes.
        matrix_values = matrix(QQ, [[1, 1], [-3, 9]])
        x1, x2 = matrix_values.solve_right(vector(ring, [QQ(-1)/2-a, node-a]))
        X = a+x1*t+x2*t**2
        Y = (t-1)*(t+3)*(b+c*t+d*t**2)
    identity = Y**2-X**3-local_polynomials(A)*X-local_polynomials(B)
    equations = tuple(ring(value) for value in identity.list() if value)
    ideal = ring.ideal(equations)
    groebner = ideal.groebner_basis()
    assert ideal.dimension() == 0
    quotient_degree = ideal.vector_space_dimension()
    print(
        f"Q80CM24QQ|profile={profile}|equations={len(equations)}|"
        f"basis={len(groebner)}|quotient_degree={quotient_degree}|groebner={tuple(groebner)}",
        flush=True,
    )
    solutions = ideal.variety(ring=quadratic)
    print(
        f"Q80CM24QQ|profile={profile}|solutions={len(solutions)}|"
        f"coordinates={tuple(tuple(solution[variable] for variable in ring.gens()) for solution in solutions)}",
        flush=True,
    )
    lex_ring = PolynomialRing(QQ, names=("b", "c", "d", "a"), order="lex")
    lex_ideal = lex_ring.ideal(
        [
            lex_ring(equation.subs(dict(zip(ring.gens(), (lex_ring("a"), lex_ring("b"), lex_ring("c"), lex_ring("d"))))))
            for equation in equations
        ]
    )
    lex_basis = lex_ideal.groebner_basis()
    a_generator = lex_ring("a")
    a_eliminants = [
        polynomial for polynomial in lex_basis
        if polynomial.variables() == (a_generator,)
    ]
    print(
        f"Q80CM24QQ|profile={profile}|a_eliminants={tuple(poly.factor() for poly in a_eliminants)}",
        flush=True,
    )
    return ring, X, Y, ideal, groebner, solutions


section_scheme("spinor")
section_scheme("identity")

# Select the three short sections singled out by the ordinary mod-7 marking.
# The corresponding residue embedding has sqrt(-3)=2 and sqrt(-6)=-1.
quadratic_three.<sqrt_minus_three> = QuadraticField(-3)
composite, embed_three, embed_six, _ = quadratic_three.composite_fields(
    quadratic, both_maps=True
)[0]
s3 = embed_three(sqrt_minus_three)
s6 = embed_six(sqrt_minus_six)
delta = QQ(3)/8*s6
function_polynomials = PolynomialRing(composite, "t")
t = function_polynomials.gen()
function_field = function_polynomials.fraction_field()
A_function = function_field(A)
B_function = function_field(B)

D1 = (
    function_field(t*(1-QQ(3)/2*t)),
    function_field(t**2*(-3*delta+2*delta*t+delta*t**2)),
)
a2 = QQ(3)/4*(1-s3)
D2 = (
    function_field(t*(1+a2*t)),
    function_field(t**2*(-3*delta+(a2*s6/2)*t-delta*t**2)),
)
a3 = QQ(9)/4*(s3-1)
d3 = -delta
c3 = a3/s6
b3 = -3*a3*c3/(a3-QQ(9)/2)
x31, x32 = matrix(composite, [[1, 1], [-3, 9]]).solve_right(
    vector(composite, [QQ(-1)/2-a3, node-a3])
)
D3 = (
    function_field(a3+x31*t+x32*t**2),
    function_field((t-1)*(t+3)*(b3+c3*t+d3*t**2)),
)
for label, point in (("D1", D1), ("D2", D2), ("D3", D3)):
    assert point[1]**2 == point[0]**3+A_function*point[0]+B_function


def add_points(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and y_left == -y_right:
        return None
    if x_left == x_right:
        slope = (3*x_left**2+A_function)/(2*y_left)
    else:
        slope = (y_right-y_left)/(x_right-x_left)
    x_sum = slope**2-x_left-x_right
    y_sum = -y_left+slope*(x_left-x_sum)
    return function_field(x_sum), function_field(y_sum)


def negate(point):
    return None if point is None else (point[0], -point[1])


def multiply_point(multiplier, point):
    if multiplier < 0:
        return multiply_point(-multiplier, negate(point))
    answer = None
    addend = point
    while multiplier:
        if multiplier & 1:
            answer = add_points(answer, addend)
        multiplier >>= 1
        if multiplier:
            addend = add_points(addend, addend)
    return answer


G1 = D1
G2 = add_points(add_points(D1, negate(D2)), negate(D3))
G3 = add_points(multiply_point(-4, D1), negate(D3))


def homogeneous_coordinates(point):
    x_coordinate, y_coordinate = point
    x_denominator = function_polynomials(x_coordinate.denominator())
    roots = x_denominator.sqrt(all=True)
    assert roots
    Z = next(
        root for root in roots
        if function_polynomials(y_coordinate.denominator()) == root**3
    )
    return (
        function_polynomials(x_coordinate.numerator()),
        function_polynomials(y_coordinate.numerator()),
        Z,
    )


def pair_intersection(left, right):
    XP, YP, ZP = homogeneous_coordinates(left)
    XQ, YQ, ZQ = homogeneous_coordinates(right)
    D = XP*ZQ**2-XQ*ZP**2
    S = YP*ZQ**3+YQ*ZP**3
    H = D*ZP*ZQ
    N = S**2-D**2*(XP*ZQ**2+XQ*ZP**2)
    first = gcd(H, N)
    cancellation = gcd(H, N//first).monic()
    H_reduced = H//cancellation
    N_reduced = N//cancellation**2
    finite = H_reduced.degree()
    excess = N_reduced.degree()-2*finite-4
    assert excess <= 0 or excess % 2 == 0
    return finite+(excess//2 if excess > 0 else 0)


coordinate_basis = (composite.one(), s3, s6, s3*s6)
basis_matrix = matrix(QQ, [vector(value) for value in coordinate_basis]).transpose()


def biquadratic_coordinates(value):
    coordinates = basis_matrix.solve_right(vector(composite(value)))
    assert sum(coefficient*basis for coefficient, basis in zip(coordinates, coordinate_basis)) == value
    return tuple(coordinates)


def polynomial_table(value):
    value = function_polynomials(value)
    return tuple(biquadratic_coordinates(value[index]) for index in range(value.degree()+1))


for label, point, expected_pole in (("G1", G1, 0), ("G2", G2, 0), ("G3", G3, 1)):
    assert point[1]**2 == point[0]**3+A_function*point[0]+B_function
    X_homogeneous, Y_homogeneous, Z_homogeneous = homogeneous_coordinates(point)
    assert Z_homogeneous.degree() == expected_pole
    print(
        f"Q80CM24QQ|section={label}|X={polynomial_table(X_homogeneous)}|"
        f"Y={polynomial_table(Y_homogeneous)}|Z={polynomial_table(Z_homogeneous)}",
        flush=True,
    )
pairs = (
    pair_intersection(G1, G2),
    pair_intersection(G1, G3),
    pair_intersection(G2, G3),
)
assert pairs == (2, 3, 1)
print(f"Q80CM24QQ|generic_basis_pairs={pairs}|status=PASS", flush=True)

# Linearize the *generic* pole profile (0,0,3) at the CM point.  Inflate the
# visible pole-one G3 by an arbitrary monic quadratic Q; its two coefficients
# are cancellation-gauge directions at the special point.  The remaining
# kernel direction must move the surface and is the desired smooth rank-19
# deformation tangent.
Q_cancel = t**2+1
X3_visible, Y3_visible, Z3_visible = homogeneous_coordinates(G3)
X3_full = X3_visible*Q_cancel**2
Y3_full = Y3_visible*Q_cancel**3
Z3_full = Z3_visible*Q_cancel
assert (X3_full.degree(), Y3_full.degree(), Z3_full.degree()) == (10, 15, 3)
h3 = Y3_full[0]/(X3_full[0]*Z3_full[0])
l3 = Y3_full[15]/X3_full[10]
assert X3_full[0] == h3**2*Z3_full[0]**2
assert Y3_full[0] == h3**3*Z3_full[0]**3
assert X3_full[10] == l3**2 and Y3_full[15] == l3**3

deformation_names = (
    ("d", "p", "q", "r", "b1", "b2", "b3", "b4", "e")
    + ("c1", "u10", "u11")
    + ("k2", "l2", "x21", "x22", "x23", "y21", "y22", "y23", "y24", "y25")
    + ("z30", "z31", "z32", "h3", "l3")
    + tuple(f"x3{index}" for index in range(1, 10))
    + tuple(f"y3{index}" for index in range(1, 15))
)
assert len(deformation_names) == 50
d_value = QQ(-1)/2
A_base = function_polynomials(A)
B_base = function_polynomials(B)
Delta_variation_factor_A = 12*A_base**2
Delta_variation_factor_B = 54*B_base
X1_base, Y1_base = map(function_polynomials, G1)
X2_base, Y2_base = map(function_polynomials, G2)


def i4_cone_linearization(x_coordinate, y_coordinate, dx, dy, dA, dd):
    x_coordinate = function_field(x_coordinate)
    y_coordinate = function_field(y_coordinate)
    dx = function_field(dx)
    dy = function_field(dy)
    eta = 6*d_value*y_coordinate.derivative()(1)
    xi = 6*d_value*x_coordinate.derivative()(1)+A_function.derivative()(1)
    delta_eta = 6*dd*y_coordinate.derivative()(1)+6*d_value*dy.derivative()(1)
    delta_xi = (
        6*dd*x_coordinate.derivative()(1)+6*d_value*dx.derivative()(1)
        +function_field(dA).derivative()(1)
    )
    return 2*eta*delta_eta-3*dd*xi**2-6*d_value*xi*delta_xi


def zero_variation():
    return {
        "dd": composite.zero(),
        "dA": function_polynomials.zero(),
        "dB": function_polynomials.zero(),
        "dX1": function_polynomials.zero(),
        "dY1": function_polynomials.zero(),
        "dX2": function_polynomials.zero(),
        "dY2": function_polynomials.zero(),
        "dX3": function_polynomials.zero(),
        "dY3": function_polynomials.zero(),
        "dZ3": function_polynomials.zero(),
    }


variations = []
for name in deformation_names:
    variation = zero_variation()
    if name == "d":
        variation["dd"] = 1
    elif name in ("p", "q", "r"):
        variation["dA"] = t**({"p": 3, "q": 4, "r": 5}[name])
    elif name in ("b1", "b2", "b3", "b4", "e"):
        variation["dB"] = t**({"b1": 4, "b2": 5, "b3": 6, "b4": 7, "e": 8}[name])
    elif name == "c1":
        variation["dX1"] = t**2
    elif name == "u10":
        variation["dY1"] = t**2*(t-1)
    elif name == "u11":
        variation["dY1"] = t**3*(t-1)
    elif name == "k2":
        k2 = Y2_base[0]/X2_base[0]
        variation["dX2"] = 2*k2
        variation["dY2"] = 3*k2**2
    elif name == "l2":
        l2 = Y2_base[6]/X2_base[4]
        variation["dX2"] = 2*l2*t**4
        variation["dY2"] = 3*l2**2*t**6
    elif name.startswith("x2"):
        variation["dX2"] = t**int(name[2:])
    elif name.startswith("y2"):
        variation["dY2"] = t**int(name[2:])
    elif name.startswith("z3"):
        variation["dZ3"] = t**int(name[2:])
    elif name == "h3":
        variation["dX3"] = 2*h3*Z3_full[0]**2
        variation["dY3"] = 3*h3**2*Z3_full[0]**3
    elif name == "l3":
        variation["dX3"] = 2*l3*t**10
        variation["dY3"] = 3*l3**2*t**15
    elif name.startswith("x3"):
        variation["dX3"] = t**int(name[2:])
    elif name.startswith("y3"):
        variation["dY3"] = t**int(name[2:])
    else:
        raise AssertionError(name)
    variations.append(variation)


def section_linearization(X, Y, Z, variation):
    dA, dB = variation["dA"], variation["dB"]
    dX, dY, dZ = variation["dX3"], variation["dY3"], variation["dZ3"]
    return (
        2*Y*dY-3*X**2*dX-dA*X*Z**4-A_base*dX*Z**4
        -4*A_base*X*Z**3*dZ-dB*Z**6-6*B_base*Z**5*dZ
    )


columns = []
for variation in variations:
    dA, dB, dd = variation["dA"], variation["dB"], variation["dd"]
    ambient_column = [
        dA(1)+6*d_value*dd,
        dB(1)-6*d_value**2*dd,
    ]
    delta_discriminant = Delta_variation_factor_A*dA+Delta_variation_factor_B*dB
    ambient_column += [delta_discriminant.derivative(order)(1) for order in range(1, 4)]
    dF1 = (
        2*Y1_base*variation["dY1"]-3*X1_base**2*variation["dX1"]
        -dA*X1_base-A_base*variation["dX1"]-dB
    )
    dF2 = (
        2*Y2_base*variation["dY2"]-3*X2_base**2*variation["dX2"]
        -dA*X2_base-A_base*variation["dX2"]-dB
    )
    dF3 = section_linearization(X3_full, Y3_full, Z3_full, variation)
    column = ambient_column
    column += [dF1[index] for index in range(9)]
    column += [dF2[index] for index in range(13)]
    column += [dF3[index] for index in range(31)]
    column += [function_polynomials(variation["dX1"])(1)-dd]
    column += [
        function_polynomials(variation["dX3"])(1)-dd*Z3_full(1)**2
        -2*d_value*Z3_full(1)*function_polynomials(variation["dZ3"])(1),
        function_polynomials(variation["dY3"])(1),
    ]
    dx3_affine = (
        variation["dX3"]*Z3_full-2*X3_full*variation["dZ3"]
    )/Z3_full**3
    dy3_affine = (
        variation["dY3"]*Z3_full-3*Y3_full*variation["dZ3"]
    )/Z3_full**4
    column += [
        i4_cone_linearization(
            X1_base, Y1_base, variation["dX1"], variation["dY1"], dA, dd
        ),
        i4_cone_linearization(
            function_field(X3_full/Z3_full**2),
            function_field(Y3_full/Z3_full**3),
            dx3_affine,
            dy3_affine,
            dA,
            dd,
        ),
    ]
    columns.append(vector(composite, column))
deformation_jacobian = matrix(composite, columns).transpose()
deformation_kernel = deformation_jacobian.right_kernel().basis_matrix()
surface_projection = deformation_kernel.matrix_from_columns(range(9))
print(
    f"Q80CM24TANGENT|stage=raw|variables=50|equations={deformation_jacobian.nrows()}|"
    f"rank={deformation_jacobian.rank()}|kernel={deformation_kernel.nrows()}|"
    f"surface_projection={surface_projection.rank()}",
    flush=True,
)
assert deformation_kernel.nrows() == 8
assert surface_projection.rank() == 2
surface_tangent_rows = []
for row in surface_projection.row_space().basis_matrix().rows():
    pivot = next(value for value in row if value)
    surface_tangent_rows.append(tuple(value/pivot for value in row))
print(
    f"Q80CM24TANGENT|variables=50|equations={deformation_jacobian.nrows()}|"
    f"rank={deformation_jacobian.rank()}|kernel=8|pure_section_kernel=6|"
    f"surface_tangent_plane=2|surface_tangents="
    f"{tuple(tuple(biquadratic_coordinates(value) for value in row) for row in surface_tangent_rows)}|"
    "status=PASS",
    flush=True,
)

# Resolve the six vertical tangent directions by returning to the minimal
# pole-one chart that is special to the CM24 point.  This is the exact
# characteristic-zero analogue of ``verify_q80_rank19_deformation_gf7.sage``:
# after imposing the D5/E6 component zeros for G1, the marked section scheme
# has only the two surface-moving tangent directions.  Computing its degree-2
# Kuranishi map therefore gives the intrinsic tangent cone of the two
# rank-19 branches, rather than a cone polluted by cancellation gauges from
# the artificially inflated pole-three presentation above.
cone_names = (
    ("D", "P", "Q", "E")
    + tuple(f"X1_{index}" for index in range(5))
    + tuple(f"Y1_{index}" for index in range(7))
    + tuple(f"X2_{index}" for index in range(5))
    + tuple(f"Y2_{index}" for index in range(7))
    + ("LAM",)
    + tuple(f"N_{index}" for index in range(7))
    + tuple(f"M_{index}" for index in range(10))
)


class Jet2:
    """A tiny exact K[epsilon]/(epsilon^3) implementation."""

    __slots__ = ("coefficients",)
    precision = 3
    base_ring = composite

    @classmethod
    def coerce_base(cls, value):
        try:
            if value.parent() is cls.base_ring:
                return value
        except AttributeError:
            pass
        return cls.base_ring(value)

    def __init__(self, constant=0, linear=0, quadratic=0):
        values = (constant, linear, quadratic)
        self.coefficients = tuple(
            self.coerce_base(values[index]) if index < len(values) else self.base_ring.zero()
            for index in range(self.precision)
        )

    @classmethod
    def from_coefficients(cls, coefficients):
        result = cls()
        values = list(coefficients)
        result.coefficients = tuple(
            cls.coerce_base(values[index]) if index < len(values) else cls.base_ring.zero()
            for index in range(cls.precision)
        )
        return result

    @property
    def constant(self):
        return self.coefficients[0]

    @property
    def linear(self):
        return self.coefficients[1]

    @property
    def quadratic(self):
        return self.coefficients[2]

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Jet2) else Jet2(value)

    def __add__(self, other):
        other = self.coerce(other)
        return self.from_coefficients(
            [left+right for left, right in zip(self.coefficients, other.coefficients)]
        )

    __radd__ = __add__

    def __neg__(self):
        return self.from_coefficients([-value for value in self.coefficients])

    def __sub__(self, other):
        return self+(-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other)-self

    def __mul__(self, other):
        other = self.coerce(other)
        return self.from_coefficients(
            [
                sum(
                    (self.coefficients[index]*other.coefficients[degree-index]
                     for index in range(degree+1)),
                    self.base_ring.zero(),
                )
                for degree in range(self.precision)
            ]
        )

    __rmul__ = __mul__

    def inverse(self):
        assert self.constant
        coefficients = [1/self.constant]
        for degree in range(1, self.precision):
            coefficients.append(
                -sum(
                    self.coefficients[index]*coefficients[degree-index]
                    for index in range(1, degree+1)
                )/self.constant
            )
        return self.from_coefficients(coefficients)

    def __truediv__(self, other):
        return self*self.coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.coerce(other)*self.inverse()

    def __pow__(self, exponent):
        assert exponent >= 0
        result = Jet2(1)
        factor = self
        while exponent:
            if exponent & 1:
                result *= factor
            exponent >>= 1
            if exponent:
                factor *= factor
        return result


def jet_poly_add(left, right):
    length = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else Jet2())
        +(right[index] if index < len(right) else Jet2())
        for index in range(length)
    ]


def jet_poly_neg(value):
    return [-coefficient for coefficient in value]


def jet_poly_mul(left, right, truncate=None):
    length = len(left)+len(right)-1
    if truncate is not None:
        length = min(length, truncate)
    result = [Jet2() for _ in range(length)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i+j < length:
                result[i+j] = result[i+j]+left_value*right_value
    return result


def jet_poly_pow(value, exponent, truncate=None):
    result = [Jet2(1)]
    for _ in range(exponent):
        result = jet_poly_mul(result, value, truncate=truncate)
    return result


def jet_poly_scale(value, scalar):
    scalar = Jet2.coerce(scalar)
    return [scalar*coefficient for coefficient in value]


def jet_poly_subtract(first, *others):
    result = first
    for other in others:
        result = jet_poly_add(result, jet_poly_neg(other))
    return result


def jet_poly_value_at_one(value):
    return sum(value, Jet2())


def jet_poly_derivative_at_one(value):
    return sum((Jet2(index)*value[index] for index in range(1, len(value))), Jet2())


jet_matrix_inverse = matrix(
    QQ, 4, 4, lambda row, column: binomial(4+column, row)
).inverse()


def evaluate_cone_system(values):
    D, P, Q, E = values[:4]
    remaining = values[4:]
    X1 = list(remaining[0:5])
    Y1 = list(remaining[5:12])
    X2 = list(remaining[12:17])
    Y2 = list(remaining[17:24])
    LAM = remaining[24]
    N = list(remaining[25:32])
    M = list(remaining[32:42])
    R = -3*D**2+3-P-Q
    A_jet = [Jet2(), Jet2(), Jet2(-3), P, Q, R]
    A_at_one = [
        sum(
            (Jet2(binomial(index, order))*A_jet[index] for index in range(order, 6)),
            Jet2(),
        )
        for order in range(4)
    ]
    denominator = -3*D**2
    u_series = list(A_at_one)
    u_series[0] = u_series[0]+3*D**2
    u_series = [coefficient/denominator for coefficient in u_series]
    branch_series = [Jet2(1), Jet2(), Jet2(), Jet2()]
    branch_series = jet_poly_add(
        branch_series, jet_poly_scale(u_series, QQ(3)/2)
    )
    branch_series = jet_poly_add(
        branch_series,
        jet_poly_scale(jet_poly_pow(u_series, 2, truncate=4), QQ(3)/8),
    )
    branch_series = jet_poly_add(
        branch_series,
        jet_poly_scale(jet_poly_pow(u_series, 3, truncate=4), -QQ(1)/16),
    )
    branch_series = jet_poly_scale(branch_series, 2*D**3)
    fixed_series = [
        Jet2(2*binomial(3, order))+E*binomial(8, order)
        for order in range(4)
    ]
    difference = [branch_series[index]-fixed_series[index] for index in range(4)]
    middle_B = [
        sum((Jet2(jet_matrix_inverse[row, column])*difference[column] for column in range(4)), Jet2())
        for row in range(4)
    ]
    B_jet = [Jet2(), Jet2(), Jet2(), Jet2(2)]+middle_B+[E]
    H = [-LAM, Jet2(1)]
    residuals = (
        jet_poly_subtract(
            jet_poly_pow(Y1, 2), jet_poly_pow(X1, 3), jet_poly_mul(A_jet, X1), B_jet
        ),
        jet_poly_subtract(
            jet_poly_pow(Y2, 2), jet_poly_pow(X2, 3), jet_poly_mul(A_jet, X2), B_jet
        ),
        jet_poly_subtract(
            jet_poly_pow(M, 2),
            jet_poly_pow(N, 3),
            jet_poly_mul(jet_poly_mul(A_jet, N), jet_poly_pow(H, 4)),
            jet_poly_mul(B_jet, jet_poly_pow(H, 6)),
        ),
    )
    equations = []
    for residual, expected_degree in zip(residuals, (12, 12, 18)):
        residual += [Jet2()]*(expected_degree+1-len(residual))
        equations.extend(residual[:expected_degree+1])
    X1_one = jet_poly_value_at_one(X1)
    Y1_one = jet_poly_value_at_one(Y1)
    N_one = jet_poly_value_at_one(N)
    M_one = jet_poly_value_at_one(M)
    A_prime = jet_poly_derivative_at_one(A_jet)
    X1_prime = jet_poly_derivative_at_one(X1)
    Y1_prime = jet_poly_derivative_at_one(Y1)
    H_one = 1-LAM
    M_affine_prime = jet_poly_derivative_at_one(M)/H_one**3-3*M_one/H_one**4
    N_affine_prime = jet_poly_derivative_at_one(N)/H_one**2-2*N_one/H_one**3
    equations.extend(
        (
            X1_one-D,
            Y1_one,
            N_one-D*H_one**2,
            M_one,
            X1[1]-1,
            (6*D*Y1_prime)**2-3*D*(6*D*X1_prime+A_prime)**2,
            (6*D*M_affine_prime)**2-3*D*(6*D*N_affine_prime+A_prime)**2,
        )
    )
    assert len(equations) == 52
    return equations


def padded_coefficients(polynomial, length):
    polynomial = function_polynomials(polynomial)
    return tuple(polynomial[index] for index in range(length))


X1_seed, Y1_seed, _ = homogeneous_coordinates(G1)
X2_seed, Y2_seed, _ = homogeneous_coordinates(G2)
N_seed, M_seed, H_seed = homogeneous_coordinates(G3)
assert H_seed[1] == 1
cone_seed_values = (
    (d_value, QQ(9)/4, -QQ(9)/4, -QQ(27)/32)
    + padded_coefficients(X1_seed, 5)
    + padded_coefficients(Y1_seed, 7)
    + padded_coefficients(X2_seed, 5)
    + padded_coefficients(Y2_seed, 7)
    + (-H_seed[0],)
    + padded_coefficients(N_seed, 7)
    + padded_coefficients(M_seed, 10)
)
assert len(cone_seed_values) == len(cone_names) == 46
assert not any(
    value.constant for value in evaluate_cone_system(list(map(Jet2, cone_seed_values)))
)

cone_jacobian_columns = []
for column in range(46):
    inputs = [Jet2(value, 1 if index == column else 0) for index, value in enumerate(cone_seed_values)]
    cone_jacobian_columns.append(
        vector(composite, [value.linear for value in evaluate_cone_system(inputs)])
    )
cone_jacobian = matrix(composite, cone_jacobian_columns).transpose()
component_zero_names = {
    "X1_0", "Y1_0", "Y1_1", "X1_3", "X1_4", "Y1_5", "Y1_6",
}
cone_active_columns = [
    index for index, name in enumerate(cone_names) if name not in component_zero_names
]
cone_active_names = tuple(cone_names[index] for index in cone_active_columns)
cone_resolved_jacobian = cone_jacobian.matrix_from_columns(cone_active_columns)
cone_tangent = cone_resolved_jacobian.right_kernel_matrix()
assert len(cone_active_names) == 39
assert cone_resolved_jacobian.rank() == 37
assert cone_tangent.nrows() == 2


def directional_second(direction):
    inputs = []
    active_lookup = {column: index for index, column in enumerate(cone_active_columns)}
    for column, seed_value in enumerate(cone_seed_values):
        linear = direction[active_lookup[column]] if column in active_lookup else 0
        inputs.append(Jet2(seed_value, linear))
    return vector(composite, [value.quadratic for value in evaluate_cone_system(inputs)])


quadratic_00 = directional_second(cone_tangent[0])
quadratic_11 = directional_second(cone_tangent[1])
quadratic_sum = directional_second(cone_tangent[0]+cone_tangent[1])
quadratic_01 = quadratic_sum-quadratic_00-quadratic_11
cone_left_kernel = cone_resolved_jacobian.left_kernel_matrix()
obstruction_00 = cone_left_kernel*quadratic_00
obstruction_01 = cone_left_kernel*quadratic_01
obstruction_11 = cone_left_kernel*quadratic_11
quadratic_ring = PolynomialRing(composite, names=("tau0", "tau1"))
tau0, tau1 = quadratic_ring.gens()
quadratic_polynomials = []
for row in range(cone_left_kernel.nrows()):
    polynomial = (
        obstruction_00[row]*tau0**2
        +obstruction_01[row]*tau0*tau1
        +obstruction_11[row]*tau1**2
    )
    if polynomial and polynomial not in quadratic_polynomials:
        quadratic_polynomials.append(polynomial)
quadratic_basis = tuple(quadratic_ring.ideal(quadratic_polynomials).groebner_basis())
assert len(quadratic_basis) == 1
quadratic_generator = quadratic_basis[0]
quadratic_factorization = quadratic_generator.factor()
assert sum(exponent for _, exponent in quadratic_factorization) == 2
surface_tangent_coordinates = tuple(
    tuple(cone_tangent[row, cone_active_names.index(name)] for name in ("D", "P", "Q", "E"))
    for row in range(2)
)
print(
    f"Q80CM24CONE|variables=39|equations=52|"
    f"rank=37|tangent_dimension=2|surface_tangents="
    f"{tuple(tuple(biquadratic_coordinates(value) for value in row) for row in surface_tangent_coordinates)}|"
    f"equation={quadratic_generator}|factorization={quadratic_factorization}|status=PASS",
    flush=True,
)
print(
    "Q80CM24QQ|scope=boundary_surface_sections_and_tangent_cone"
    "|generic_one_parameter_branch=NOT_CONSTRUCTED",
    flush=True,
)
print("Q80CM24QQ|status=PASS", flush=True)
