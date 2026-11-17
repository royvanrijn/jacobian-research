"""Exact controls for the unpublished Mathieu manuscript (SymPy required).

Run from the repository root with::

    .venv/bin/python unpublished/mathieu-reductive-homogeneous-spaces/verify.py

Polynomial identities are checked exactly. Moment and congruence checks use
only the finite orders listed in the output; they do not prove the all-order
claims or the classification of reductive homogeneous spaces in main.tex.
"""

import json
from fractions import Fraction
from math import comb, factorial

import sympy as S


X, Y, Z, W = S.symbols("X Y Z W")


def require(condition, message):
    """Keep verification active even when Python is invoked with -O."""
    if not condition:
        raise AssertionError(message)


def binary_form(n, x=X, y=Y):
    if n == 4:
        return x**2 * y**2
    if n == 8:
        return x**8 + 14*x**4*y**4 + y**8
    if n == 12:
        return x*y*(x**10 + 11*x**5*y**5 - y**10)
    raise ValueError(f"Unsupported binary degree: {n}")


def check_klein_identities():
    result = {}
    for n in (8, 12):
        f = binary_form(n)
        if n == 8:
            g = X**2*Y**2*(X**4-Y**4)**2
            t = X**12-33*X**8*Y**4-33*X**4*Y**8+Y**12
            stated_j = X*Y*(X**4-Y**4)*t/6
        else:
            g = X**20+Y**20-228*X**15*Y**5+494*X**10*Y**10+228*X**5*Y**15
            t = X**30+Y**30+522*(X**25*Y**5-X**5*Y**25)-10005*(X**20*Y**10+X**10*Y**20)
            stated_j = -t/12
        v = S.Poly(g, X, Y).total_degree()
        j = (S.diff(f, X)*S.diff(g, Y)-S.diff(f, Y)*S.diff(g, X))/(n*v)
        require(S.Poly(j-stated_j, X, Y).is_zero, f"Degree {n}: Jacobian")
        relation = j**2-f**3*g/36+3*g**3 if n == 8 else j**2-12*f**5-g**3/144
        require(S.Poly(relation, X, Y).is_zero, f"Degree {n}: Klein relation")
        result[str(n)] = "PASS"
    return result


def require_octahedral_invariance(form):
    transformed = form.xreplace({
        X: (X+S.I*Y)/S.sqrt(2), Y: (S.I*X+Y)/S.sqrt(2),
    })
    require(S.Poly(transformed-form, X, Y, extension=[S.sqrt(2), S.I]).is_zero,
            "Octahedral quarter-turn generator")


def check_invariance():
    # The diagonal generators have eigenvalues exp(+-pi*i/4) and
    # exp(+-pi*i/5): check their action by exponent differences.
    for n, order in ((8, 8), (12, 10)):
        require(all((i-j) % order == 0 for (i, j), _ in S.Poly(binary_form(n), X, Y).terms()),
                f"Degree {n}: diagonal generator")
    require_octahedral_invariance(binary_form(8))
    alpha = S.sqrt((5+S.sqrt(5))/10)
    beta = S.sqrt((5-S.sqrt(5))/10)
    require(S.simplify(alpha**2+beta**2) == 1, "Icosahedral generator determinant")
    # The SU(2) matrix is i*[[beta,alpha],[alpha,-beta]]. Since i**12=1,
    # its common scalar i can be removed in this degree-twelve substitution.
    difference = binary_form(12, beta*X+alpha*Y, alpha*X-beta*Y)-binary_form(12)
    require(S.Poly(difference, X, Y, extension=[alpha, beta]).is_zero,
            "Icosahedral non-diagonal generator")
    p4 = X**2*Y**2*(1+X*Z*W**2)**4*(1+Y*Z**2*W)**4
    q4 = Z**2*W**2
    weyl = {X: Y, Y: -X, Z: W, W: -Z}
    require(S.expand(p4.xreplace(weyl)-p4) == 0, "Normalizer witness Weyl invariance")
    require(S.expand(q4.xreplace(weyl)-q4) == 0, "Normalizer marker Weyl invariance")
    return "PASS: diagonal and non-diagonal generators for each exceptional form; normalizer Weyl element"


def univariate_coefficients(expression):
    return {i: int(c) for (i,), c in S.Poly(expression, X).terms() if c}


def haar_moment(n, m, marked):
    """Integrate the zero-weight component by the monomial beta formula.

    This does not use a Klein identity, a modular congruence or the radial
    closed form. Taylor coefficients are formed one top-row monomial at a
    time, avoiding the much larger expansion of the complete witness P**m.
    """
    b = n//2
    a = b-1
    upper_power, lower_power = a*m, b*m
    s = lower_power-b if marked else lower_power
    f = binary_form(n, X, 1)
    upper = univariate_coefficients(f**upper_power)
    lower = univariate_coefficients(f**lower_power)
    total_lower_degree = n*lower_power
    row_degree = n*upper_power+s
    answer = Fraction(0)
    for r in range(s+1):
        lower_coefficient = sum(
            coefficient*comb(k, r)*comb(total_lower_degree-k, s-r)*X**(k-r)
            for k, coefficient in lower.items()
            if k >= r and total_lower_degree-k >= s-r
        )
        lower_polynomial = S.Poly(lower_coefficient, X, domain=S.ZZ)
        if marked:
            lower_polynomial, remainder = lower_polynomial.div(S.Poly(f, X, domain=S.ZZ))
            require(remainder.is_zero, f"Degree {n}, m={m}: polynomial division by F")
        coefficients = univariate_coefficients(lower_polynomial.as_expr())
        for j, coefficient in upper.items():
            top_x = j+r
            bottom_z = row_degree-top_x
            if bottom_z not in coefficients:
                continue
            answer += Fraction(
                coefficient*coefficients[bottom_z]*(-1)**bottom_z
                *factorial(top_x)*factorial(bottom_z),
                factorial(row_degree+1),
            )
    return answer


def valuation(q, p):
    def integer_valuation(value):
        value = abs(value)
        count = 0
        while value and value % p == 0:
            value //= p
            count += 1
        return count
    require(q != 0, "Valuation of zero is not finite")
    return integer_valuation(q.numerator)-integer_valuation(q.denominator)


def coefficient_mod_p(n, p):
    """Extract [X**(p-1)]f(X,1)**((p-1)/b) by sparse multiplication."""
    base = univariate_coefficients(binary_form(n, X, 1))
    coefficients = {0: 1}
    for _ in range((p-1)//(n//2)):
        product = {}
        for i, left in coefficients.items():
            for j, right in base.items():
                if i+j <= p-1:
                    product[i+j] = (product.get(i+j, 0)+left*right) % p
        coefficients = product
    return coefficients.get(p-1, 0)


def hasse_closed_value(n, p):
    if n == 8:
        return comb((p-1)//2, (p-1)//4)*pow(-pow(12, -1, p), (p-1)//4, p) % p
    return (comb((p-1)//2, (p-1)//3)*pow(12, (p-1)//6, p)
            *pow(pow(144, -1, p), (p-1)//3, p)) % p


def check_exceptional_moments():
    result = {}
    for n in (8, 12):
        a, b = n//2-1, n//2
        c = b*(n-1)
        orders = {}
        for m in (1, 2):
            require(haar_moment(n, m, False) == 0, f"Degree {n}, m={m}: purity")
            marked = haar_moment(n, m, True)
            p = c*m-a
            require(factorial(p) % marked.denominator == 0,
                    f"Degree {n}, m={m}: factorial denominator bound")
            entry = {"pure": "0", "marked": str(marked)}
            if S.isprime(p):
                require(valuation(marked, p) == -1, f"Degree {n}, m={m}: p-valuation")
                scaled = marked*p
                residue = scaled.numerator*pow(scaled.denominator, -1, p) % p
                expected = comb(n*b*m, b*(m-1))*coefficient_mod_p(n, p) % p
                require(residue == expected != 0, f"Degree {n}, m={m}: Haar congruence")
                entry.update({"prime": p, "valuation": -1, "scaled_moment_mod_prime": residue})
            orders[str(m)] = entry
        result[str(n)] = orders
    return result


def check_binary_coefficients():
    result = {}
    for n in (8, 12):
        b = n//2
        a, c = b-1, b*(n-1)
        primes = sorted(set(S.primerange(7, 101)).intersection(
            p for p in range(7, 101) if (p-1) % b == 0
        ) | {c*m-a for m in range(1, 9) if S.isprime(c*m-a)})
        for p in primes:
            f = S.Poly(binary_form(n, X, 1), X, modulus=p)
            require(S.gcd(f, f.diff()).degree() == 0,
                    f"Degree {n}, p={p}: simple binary roots")
            require(coefficient_mod_p(n, p) == hasse_closed_value(n, p) != 0,
                    f"Degree {n}, p={p}: binary coefficient")
        result[str(n)] = primes
    return result


def integer_beta(left_power, right_power):
    """Integral of t**left_power * (1-t)**right_power on [0, 1]."""
    return Fraction(factorial(left_power)*factorial(right_power),
                    factorial(left_power+right_power+1))


def radial_phase_moment(m, marked):
    """Integrate the explicit SU(2) phase coefficient, independently of rows.

    For t=|x|**2 and u=conjugate(x*y), P=(xy)**2 times
    (1-t*u)**4 * (1+(1-t)*u)**4 and Q=u**2. Phase averaging
    selects u**(2*m) or u**(2*m-2); integrate each binomial summand
    against t**(2*m)*(1-t)**(2*m) by the beta integral.
    """
    phase_degree = 2*m-2 if marked else 2*m
    return sum(((-1)**r*comb(4*m, r)*comb(4*m, phase_degree-r)
                *integer_beta(2*m+r, 2*m+phase_degree-r)
                for r in range(phase_degree+1)), Fraction(0))


def radial_sum(m):
    length = 2*m-2
    return sum((Fraction((-1)**r*comb(length, r),
                         (4*m-r)*(4*m-r-1)*(2*m+r+1)*(2*m+r+2))
                for r in range(length+1)), Fraction(0))


def radial_double_integral(m):
    """Integrate S_m over 0<s<t<1 after s=r*t, then double.

    Expand only (1-r*t)*(1-t). Its four terms leave elementary
    t-integrals and two beta integrals in r. This does not use the
    denominator sum or its claimed simplified value.
    """
    length = 2*m-2
    first = integer_beta(2*m, length)
    second = integer_beta(2*m+1, length)
    return 2*(first*(Fraction(1, 6*m)-Fraction(1, 6*m+1))
              -second*(Fraction(1, 6*m+1)-Fraction(1, 6*m+2)))


def require_radial_moments(order, claimed_marked):
    """Compare a claimed value with row, phase and double-integral routes."""
    length = 2*order-2
    beta_sum = radial_sum(order)
    double_integral = radial_double_integral(order)
    evaluated_sum = (integer_beta(2*order, length)
                     / (6*order*(6*order+2)))
    require(beta_sum == double_integral == evaluated_sum > 0,
            f"Radial beta/double-integral evaluation m={order}")
    prefactor = Fraction(factorial(4*order)**2,
                         factorial(6*order-1)*factorial(length))
    direct = haar_moment(4, order, True)
    phase_marked = radial_phase_moment(order, True)
    require(prefactor*beta_sum == claimed_marked == direct == phase_marked > 0,
            f"Radial marked moment m={order}")
    require(haar_moment(4, order, False) == radial_phase_moment(order, False) == 0,
            f"Radial pure moment m={order}")
    return direct


def check_radial_identities():
    f = Z**2*W**2
    quotient = (Z+f*X)**2*(W+f*Y)**2/f
    require(S.cancel(quotient-(1+X*Z*W**2)**2*(1+Y*Z**2*W)**2) == 0,
            "Normalizer polynomial factorization")
    m = S.symbols("m", integer=True, positive=True)
    triangle_factor = 2*(1/(6*m*(6*m+1))
                         -(2*m+1)/(4*m*(6*m+1)*(6*m+2)))
    require(S.cancel(triangle_factor-1/(6*m*(6*m+2))) == 0,
            "Symbolic double-integral factor simplification")
    prefactor = S.factorial(4*m)**2/S.factorial(6*m-1)/S.factorial(2*m-2)
    evaluated_sum = S.factorial(2*m-2)*S.factorial(2*m)/(6*m*(6*m+2)*S.factorial(4*m-1))
    closed = 2*m/(3*m+1)*S.factorial(2*m)*S.factorial(4*m)/S.factorial(6*m)
    require(S.combsimp(prefactor*evaluated_sum-closed) == 0,
            "Symbolic factorial simplification of evaluated radial sum")
    moments = {}
    for order in range(1, 9):
        closed_value = Fraction(2*order*factorial(2*order)*factorial(4*order),
                                (3*order+1)*factorial(6*order))
        moments[str(order)] = str(require_radial_moments(order, closed_value))
    return {"polynomial_and_rational_identities": "PASS",
            "symbolic_factorial_simplification": "PASS",
            "independent_phase_beta_moments": "PASS: pure and marked, m=1,...,8",
            "beta_double_integral_evaluation": "PASS: m=1,...,8",
            "direct_moments": moments}


def check_negative_controls():
    def must_reject(check, description):
        try:
            check()
        except AssertionError:
            return "REJECTED"
        raise AssertionError(f"Negative control accepted: {description}")

    # Exercise the same invariance and independent-moment comparison paths
    # used above, with substantive mathematical inputs deliberately changed.
    return {
        "octahedral_middle_coefficient_14_to_13": must_reject(
            lambda: require_octahedral_invariance(binary_form(8)-X**4*Y**4),
            "corrupted octahedral form"),
        "radial_marked_m1_claim_1_over_31": must_reject(
            lambda: require_radial_moments(1, Fraction(1, 31)),
            "corrupted first marked moment"),
    }


def main():
    checks = {
        "Jacobian_and_Klein_identities": check_klein_identities(),
        "generator_invariance": check_invariance(),
        "exceptional_direct_haar_moments": check_exceptional_moments(),
        "binary_coefficient_prime_checks": check_binary_coefficients(),
        "radial_checks": check_radial_identities(),
        "negative_controls": check_negative_controls(),
    }
    print(json.dumps({
        "exact_checks": checks,
        "scope": ("Exact polynomial and rational identities; direct exceptional Haar moments "
                  "at m=1,2; binary coefficient checks at the listed primes; radial row/phase "
                  "moments and beta/double-integral evaluations at m=1,...,8; rejection of "
                  "the two listed corrupted inputs. The all-order proofs and reductive classification "
                  "are in main.tex and are not certified by these finite controls."),
        "formal_verification": False,
        "external_review": False,
    }, indent=2))


if __name__ == "__main__":
    main()
