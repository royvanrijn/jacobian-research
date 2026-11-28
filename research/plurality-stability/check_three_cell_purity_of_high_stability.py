"""Check the rational enclosures in the circular-purity comparison."""

from fractions import Fraction
from math import comb, factorial, isqrt


Q = Fraction
PI_LOWER = Q("3.14159")
PI_UPPER = Q("3.14160")


def require(condition, description):
    if not condition:
        raise ArithmeticError(description)


def sqrt_bounds(value):
    scale = 10**10
    lower = Q(isqrt(scale**2 * value.numerator // value.denominator), scale)
    return lower, lower + Q(1, scale)


def exp_negative_bounds(value):
    lower = sum((-value) ** j / factorial(j) for j in range(20))
    return lower, lower + value**20 / factorial(20)


def log_bounds(value):
    ratio = (value - 1) / (value + 1)
    lower = 2 * sum(ratio ** (2 * j + 1) / (2 * j + 1) for j in range(10))
    remainder = 2 * ratio**21 / (21 * (1 - ratio**2))
    return lower, lower + remainder


def gamma_cdf_bounds(value):
    root_lower, root_upper = sqrt_bounds(value)
    series_lower = sum(
        (-value) ** j / (factorial(j) * (j + Q(3, 2))) for j in range(20)
    )
    series_upper = series_lower + value**20 / (factorial(20) * Q(43, 2))
    lower = 2 * value * root_lower * series_lower / sqrt_bounds(PI_UPPER)[1]
    upper = 2 * value * root_upper * series_upper / sqrt_bounds(PI_LOWER)[0]
    return lower, upper


def tilt_density_bounds(correlation, radius_squared):
    complement = 1 - correlation**2
    kappa = correlation**2 / complement
    root_lower, root_upper = sqrt_bounds(radius_squared)
    complement_lower, complement_upper = sqrt_bounds(complement)
    exponential_lower, exponential_upper = exp_negative_bounds(kappa * radius_squared)
    lower = (
        2 * root_lower * exponential_lower
        / (sqrt_bounds(PI_UPPER)[1] * complement * complement_upper)
    )
    upper = (
        2 * root_upper * exponential_upper
        / (sqrt_bounds(PI_LOWER)[0] * complement * complement_lower)
    )
    return lower, upper


def arcsine_sum(value):
    return sum(
        Q(comb(2 * j, j), 4**j) * value ** (2 * j + 1) / (2 * j + 1) for j in range(7)
    )


def check_high_stability_floor():
    for value, low, high in [
        (Q(3, 2), ".4054", ".4055"),
        (Q(3), "1.0986", "1.0987"),
        (Q(1000, 357), "1.030", "1.031"),
    ]:
        lower, upper = log_bounds(value)
        require(Q(low) < lower < upper < Q(high), "logarithm enclosure")
    radial_cap = (
        4 * Q(".4055") ** 2 + (Q("1.0987") - 2 * Q(".4054")) ** 2 + Q("1.0987") ** 2
    ) / 9
    require(radial_cap < Q(".217"), "radial degree-two cap")
    gamma_lower, gamma_upper = gamma_cdf_bounds(Q("1.030"))
    require(Q(".4399") < gamma_lower < gamma_upper < Q(".4401"), "gamma tail in PH6")
    require(16 / (9 * PI_LOWER) * (1 - Q(".4399")) ** 2 < Q(".178"), "first-moment cap")

    left, right = Q(".63"), Q(".68")
    upper_at_right = (
        Q(1, 3) + Q(2, 3) * right**4
        + Q(".455") * (right**2 - right**4) + Q(".178") * (right - right**2)
    )
    derivative_at_left = Q(".178") + 2 * (Q(".455") - Q(".178")) * left
    derivative_at_left += 4 * (Q(2, 3) - Q(".455")) * left**3
    require(upper_at_right < Q(".627716"), "polynomial upper bound")
    require(derivative_at_left > Q(".73"), "polynomial derivative")
    require(Q(".73") ** 2 < 1 - right**2 and Q(".94") ** 2 < 1 - right**2 / 4,
            "simplex derivative denominators")
    require(PI_UPPER / 2 + right / Q(".73") < Q("2.51"), "inverse cosine cap")
    require(
        3 / (2 * PI_LOWER**2) * (Q("2.51") / Q(".73") + PI_UPPER / (4 * Q(".94")))
        < Q(".66"),
        "simplex derivative upper bound",
    )
    half = right / 2
    require(arcsine_sum(right) > Q(".74769"), "first inverse sine lower bound")
    require(
        Q(".3469") < arcsine_sum(half)
        and arcsine_sum(half) + half**15 / (15 * (1 - half**2)) < Q(".347"),
        "second inverse sine enclosure",
    )
    simplex_lower = (
        Q(1, 3) + 3 / (4 * PI_UPPER**2) * (Q(".74769") + Q(".3469"))
        * (PI_LOWER + Q(".74769") - Q(".347"))
    )
    require(simplex_lower > Q(".6279"), "simplex value lower bound")


def check_tilt_transfer():
    for correlation, radius_squared, low, high in [
        (".63", ".163", ".8737", ".8738"),
        (".63", ".164", ".8758", ".8759"),
        (".68", ".164", "1.0067", "1.0069"),
        (".63", "1.679", "1.0339", "1.0341"),
        (".68", "1.679", ".8751", ".8753"),
        (".68", "1.680", ".8746", ".8748"),
    ]:
        lower, upper = tilt_density_bounds(Q(correlation), Q(radius_squared))
        require(Q(low) < lower < upper < Q(high), "tilted density crossing")
    for value, low, high in [(".164", ".8487", ".8488"), ("1.679", ".1865", ".1866")]:
        lower, upper = exp_negative_bounds(Q(value))
        require(Q(low) < lower < upper < Q(high), "exponential tail enclosure")
    for value, low, high in [
        (Q(".163") / (1 - Q(".63") ** 2), ".0900", ".0902"),
        (Q("1.680") / (1 - Q(".68") ** 2), ".8998", ".9000"),
    ]:
        lower, upper = gamma_cdf_bounds(value)
        require(Q(low) < lower < upper < Q(high), "gamma tail enclosure")
    deficit = Q(7, 8) * (1 - Q(".8487") + Q(".1866")) - Q(".0900") - (1 - Q(".9000"))
    require(deficit < Q(".105663"), "tilted mixing deficit")
    require(Q(7, 8) * Q(119, 500) - Q(2, 3) * deficit > Q(11, 80), "tilted mixing floor")


def scaled_fourier_coefficient(order):
    if order == 0:
        return 0
    cosine_mod_four = (1, 0, -1, 0)[order % 4]
    cosine_mod_three = 1 if order % 3 == 0 else -Q(1, 2)
    return Q(2 * (-1) ** order - 1 + 8 * cosine_mod_four - 9 * cosine_mod_three, order**2)


def check_mixed_remainder():
    moments = [
        sum(comb(n, k) * scaled_fourier_coefficient(abs(n - 2 * k)) for k in range(n + 1))
        / 2**n
        for n in range(1, 5)
    ]
    require(moments == [Q(3, 2), -Q(5, 16), Q(19, 24), -Q(53, 256)], "cosine moments")
    coefficient = Q(19, 144) - Q(265, 6144)
    require(Q(3, 2) - Q(5, 32) ** 2 / (4 * coefficient) > Q(143, 100),
            "small angular parameter reserve")
    exponential_lower = sum(Q(7, 2) ** j / factorial(j) for j in range(8))
    require((exponential_lower - 1) / 160 > Q(143, 100) / PI_LOWER**2,
            "large angular parameter reserve")
    require(PI_UPPER**2 * Q("1.005") * Q("1.245") < Q("12.5"), "derivative cap")
    require(Q(143, 100) / PI_UPPER**2 - Q("12.5") * Q(".005") > Q(".08"),
            "mixed reserve near the edge")
    require(
        4 / Q("1.005") ** 2 - 3 > Q(".92")
        and 1 + (1 + Q(".32")) / Q(".48") < 4
        and Q(".92") * (exponential_lower - 4) / 160 > Q(".08"),
        "large angular parameter reserve near the edge",
    )


if __name__ == "__main__":
    check_high_stability_floor()
    check_tilt_transfer()
    check_mixed_remainder()
    print("Verified the rational purity floors and the mixed-pair remainder bounds.")
