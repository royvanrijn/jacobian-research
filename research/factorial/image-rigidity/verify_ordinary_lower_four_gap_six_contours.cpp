#include <algorithm>
#include <cassert>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <vector>

#if defined(__FAST_MATH__) || defined(NDEBUG)
#error "Compile with assertions enabled and without fast-math"
#endif

namespace {

static_assert(std::numeric_limits<double>::is_iec559 && FLT_EVAL_METHOD == 0);
constexpr double infinity = std::numeric_limits<double>::infinity();

double downward(double value) { return std::nextafter(value, -infinity); }
double upward(double value) { return std::nextafter(value, infinity); }

struct Interval {
    double low, high;

    Interval(double point = 0) : low(point), high(point) {}
    Interval(double lower, double upper) : low(lower), high(upper) {}

    double width() const { return high - low; }
    double midpoint() const { return (low + high) / 2; }
};

Interval operator+(Interval a, Interval b) {
    return {downward(a.low + b.low), upward(a.high + b.high)};
}

Interval operator-(Interval a, Interval b) {
    return {downward(a.low - b.high), upward(a.high - b.low)};
}

Interval operator-(Interval a) { return {-a.high, -a.low}; }

Interval operator*(Interval a, Interval b) {
    const double products[] = {
        a.low * b.low, a.low * b.high, a.high * b.low, a.high * b.high
    };
    auto limits = std::minmax_element(products, products + 4);
    return {downward(*limits.first), upward(*limits.second)};
}

Interval inverse(Interval value) {
    assert(value.low > 0 || value.high < 0);
    return {downward(1 / value.high), upward(1 / value.low)};
}

Interval operator/(Interval a, Interval b) { return a * inverse(b); }

Interval square(Interval value) {
    const double left = value.low * value.low;
    const double right = value.high * value.high;
    const double lower = value.low > 0 ? downward(left)
        : value.high < 0 ? downward(right) : 0;
    return {lower, upward(std::max(left, right))};
}

Interval rational(int numerator, int denominator) {
    return Interval(numerator) / Interval(denominator);
}

Interval exponential(Interval argument) {
    // (1 + x/M)^M <= exp(x) <= (1 - x/M)^(-M) for |x| < M.
    const int exponent = 24;
    const Interval step = argument / Interval(1 << exponent);
    Interval lower = 1 + step;
    Interval upper_inverse = 1 - step;
    assert(lower.low > 0 && upper_inverse.low > 0);
    for (int i = 0; i < exponent; ++i) {
        lower = square(lower);
        upper_inverse = square(upper_inverse);
    }
    return {lower.low, inverse(upper_inverse).high};
}

struct ComplexBox {
    Interval real, imag;

    ComplexBox(Interval x = 0, Interval y = 0) : real(x), imag(y) {}
    explicit ComplexBox(double x) : real(x), imag(0) {}
};

ComplexBox operator+(ComplexBox a, ComplexBox b) {
    return {a.real + b.real, a.imag + b.imag};
}

ComplexBox operator-(ComplexBox a, ComplexBox b) {
    return {a.real - b.real, a.imag - b.imag};
}

ComplexBox operator*(ComplexBox a, ComplexBox b) {
    return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
}

ComplexBox operator*(ComplexBox value, Interval scale) {
    return {value.real * scale, value.imag * scale};
}

ComplexBox operator/(ComplexBox value, Interval scale) {
    return value * inverse(scale);
}

Interval norm_squared(ComplexBox value) {
    return square(value.real) + square(value.imag);
}

ComplexBox inverse(ComplexBox value) {
    return ComplexBox(value.real, -value.imag) / norm_squared(value);
}

ComplexBox operator/(ComplexBox a, ComplexBox b) { return a * inverse(b); }

ComplexBox power(ComplexBox base, int exponent) {
    ComplexBox result(1);
    for (; exponent > 0; exponent >>= 1) {
        if (exponent & 1) result = result * base;
        if (exponent > 1) base = base * base;
    }
    return result;
}

int binomial(int n, int k) {
    if (k > n) return 0;
    int result = 1;
    for (int i = 1; i <= k; ++i) result = result * (n - i + 1) / i;
    return result;
}

enum class Piece { local, horizontal, split_vertical, left_vertical, origin };
struct Box { Interval x, y, t; };

const Interval split = rational(63, 10);
const Interval local_radius = rational(2, 5);

bool near_horizontal_bound(ComplexBox saddle, ComplexBox saddle_inverse,
                           Interval t, Interval modulus_limit_squared) {
    // Expand B-e^t through order 10 where B=(s/r)^4 D/6 and s=r+t.
    // The constant and linear coefficients cancel exactly.
    ComplexBox inverse_power = saddle_inverse * saddle_inverse;
    ComplexBox sum(0);
    Interval t_power = 1;
    Interval exp_coefficient = rational(1, 2);
    for (int k = 2; k <= 10; ++k) {
        ComplexBox coefficient = (
            (ComplexBox(10) - saddle) * Interval(binomial(4, k))
            + (saddle - ComplexBox(4)) * Interval(binomial(10, k))
        ) * inverse_power / Interval(6);
        sum = sum + (coefficient - ComplexBox(exp_coefficient)) * t_power;
        t_power = t_power * t;
        inverse_power = inverse_power * saddle_inverse;
        exp_coefficient = exp_coefficient / Interval(k + 1);
    }
    const Interval exp_tail = Interval(-1, 1) / Interval(60000);
    sum.real = sum.real + exp_tail;
    const ComplexBox difference = sum * exponential(-t);
    const Interval correction = 2 * difference.real + square(t) * norm_squared(difference);
    return correction.high < 0
        && (1 + square(t) * correction).high < modulus_limit_squared.low;
}

bool certify(Box box, bool upper, Piece piece) {
    const auto [x, y, t] = box;
    const ComplexBox saddle(x, y);
    const Interval radius = upper ? rational(113, 100) : rational(55, 100);
    const Interval distance_squared = square(x - (upper ? 10 : 4)) + square(y);
    if (distance_squared.low > square(radius).high) return true;

    // Multiply by |r-4|^2 to include the lower disk's pole in the exclusion test.
    const Interval norm = norm_squared(saddle);
    const Interval norm_squared_again = square(norm);
    const Interval norm_sixth = square(norm_squared_again) * norm_squared_again;
    const Interval numerator = norm_sixth * norm_squared(ComplexBox(10) - saddle);
    const Interval parameter = upper ? 112000 : 88500;
    const Interval denominator = square(parameter) * norm_squared(saddle - ComplexBox(4));
    if (upper ? numerator.low > denominator.high : numerator.high < denominator.low) {
        return true;
    }

    const ComplexBox saddle_inverse = inverse(saddle);
    const ComplexBox saddle_minus_four = saddle - ComplexBox(4);
    if (piece == Piece::local) {
        const ComplexBox curvature = (saddle - ComplexBox(5)) * (saddle - ComplexBox(8))
            * saddle_inverse * saddle_inverse;
        const Interval curvature_lower = upper ? rational(42, 1000) : rational(74, 1000);
        const Interval curvature_upper = upper ? rational(14, 100) : rational(38, 100);
        if (curvature.real.low < curvature_lower.high ||
            norm_squared(curvature).high > square(curvature_upper).low) return false;

        const ComplexBox ratio = ComplexBox(1) + saddle_inverse * t;
        const ComplexBox numerator_part = saddle_minus_four * power(ratio, 6);
        const ComplexBox divisor = ComplexBox(10) - saddle + numerator_part;
        if (norm_squared(divisor).low <= 0) return false;
        const ComplexBox eta = numerator_part * Interval(6) / divisor;
        const ComplexBox third = (ComplexBox(8) + eta * Interval(20)
            - eta * eta * Interval(15) + eta * eta * eta * Interval(2))
            * power(inverse(saddle + ComplexBox(t)), 3);
        const Interval limit = upper ? rational(75, 1000) : rational(21, 100);
        return norm_squared(third).high < square(limit).low;
    }

    ComplexBox displacement;
    if (piece == Piece::horizontal) {
        const Interval real_position = x + t;
        if (real_position.high < (upper ? split.low : 2)
            || real_position.low > (upper ? 20 : split.high)) return true;
        displacement = ComplexBox(t);
    } else if (piece == Piece::split_vertical) {
        displacement = ComplexBox(split - x, (t - 1) * y);
    } else if (piece == Piece::left_vertical) {
        displacement = ComplexBox(2 - x, (t - 1) * y);
    } else {
        displacement = ComplexBox(t - x, -y);
    }

    const Interval limit_squared = square(rational(399, 400));
    if (piece == Piece::horizontal && t.low >= -2 && t.high <= 2
        && near_horizontal_bound(saddle, saddle_inverse, t, limit_squared)) return true;

    const ComplexBox ratio = ComplexBox(1) + displacement * saddle_inverse;
    const ComplexBox normalized = power(ratio, 4)
        * (ComplexBox(10) - saddle + saddle_minus_four * power(ratio, 6)) / Interval(6);
    const Interval modulus_squared =
        norm_squared(normalized) * exponential(-2 * displacement.real);
    return modulus_squared.high < limit_squared.low;
}

void certify_piece(bool upper, Piece piece, Interval path_parameter) {
    // The initial rectangles strictly contain the displayed half disks.
    const Interval x = upper ? Interval(8.86, 11.14) : Interval(3.44, 4.56);
    const Interval y = upper ? Interval(0, 1.14) : Interval(0, 0.56);
    std::vector<std::pair<Box, int>> remaining = {{{x, y, path_parameter}, 0}};
    int box_count = 0;
    int maximum_depth = 0;
    while (!remaining.empty()) {
        auto [box, depth] = remaining.back();
        remaining.pop_back();
        ++box_count;
        maximum_depth = std::max(maximum_depth, depth);
        if (certify(box, upper, piece)) continue;
        assert(depth < 60);

        Interval* largest = &box.x;
        if (box.y.width() > largest->width()) largest = &box.y;
        if (box.t.width() > largest->width()) largest = &box.t;
        const Interval original = *largest;
        const double midpoint = original.midpoint();
        assert(original.low < midpoint && midpoint < original.high);
        largest->high = midpoint;
        remaining.emplace_back(box, depth + 1);
        *largest = original;
        largest->low = midpoint;
        remaining.emplace_back(box, depth + 1);
    }
    std::cout << (upper ? "upper" : "lower") << ' ' << static_cast<int>(piece)
        << ' ' << box_count << ' ' << maximum_depth << '\n';
}

double magnitude_upper(ComplexBox value) {
    return upward(std::sqrt(norm_squared(value).high));
}

double magnitude_lower(ComplexBox value) {
    return std::max(0.0, downward(std::sqrt(std::max(0.0, norm_squared(value).low))));
}

ComplexBox stationary_box(ComplexBox parameter, bool upper) {
    const std::complex<double> w(parameter.real.midpoint(), parameter.imag.midpoint());
    std::complex<double> estimate(upper ? 10 : 4, 0);
    for (int step = 0; step < 50; ++step) {
        auto fifth = estimate * estimate * estimate * estimate * estimate;
        auto sixth = fifth * estimate;
        estimate -= ((estimate - 10.0) * sixth + w * (estimate - 4.0))
            / (7.0 * sixth - 60.0 * fifth + w);
    }
    assert(std::isfinite(estimate.real()) && std::isfinite(estimate.imag()));

    // Verify the proposed center by comparing the shifted polynomial with its linear term.
    const ComplexBox center(estimate.real(), estimate.imag());
    const Interval radius = rational(1, 100000);
    const Interval original_radius = upper ? rational(113, 100) : rational(55, 100);
    assert((Interval(magnitude_upper(center - ComplexBox(upper ? 10 : 4))) + radius).high
           < original_radius.low);
    const ComplexBox fifth = power(center, 5);
    const ComplexBox sixth = fifth * center;
    const ComplexBox constant = (center - ComplexBox(10)) * sixth
        + parameter * (center - ComplexBox(4));
    const ComplexBox linear = sixth * Interval(7) - fifth * Interval(60) + parameter;
    Interval error(magnitude_upper(constant));
    Interval radius_power = square(radius);
    for (int k = 2; k <= 7; ++k) {
        ComplexBox term = power(center, 7 - k) * Interval(binomial(7, k));
        if (k <= 6) term = term - power(center, 6 - k) * Interval(10 * binomial(6, k));
        error = error + Interval(magnitude_upper(term)) * radius_power;
        radius_power = radius_power * radius;
    }
    assert(error.high < (radius * Interval(magnitude_lower(linear))).low);
    return center + ComplexBox(Interval(-radius.high, radius.high),
                               Interval(-radius.high, radius.high));
}

Interval height_squared(ComplexBox saddle, ComplexBox parameter) {
    return norm_squared(power(saddle, 4) * (parameter + power(saddle, 6)))
        * exponential(-2 * saddle.real);
}

void certify_boundaries() {
    int samples = 0;
    const int subdivisions = 32;
    const Interval required_ratio_squared = square(rational(1035, 1000));
    for (int sign : {-1, 1}) {
        for (int j = 0; j <= subdivisions; ++j) {
            const Interval t = rational(j, subdivisions);
            const Interval divisor = 1 + square(t);
            const ComplexBox unit((1 - square(t)) * Interval(sign) / divisor, 2 * t / divisor);
            for (int radius : {88500, 112000}) {
                const ComplexBox parameter = unit * Interval(radius);
                const Interval upper = height_squared(stationary_box(parameter, true), parameter);
                const Interval lower = height_squared(stationary_box(parameter, false), parameter);
                if (radius == 88500) assert(upper.low > (required_ratio_squared * lower).high);
                else assert(lower.low > (required_ratio_squared * upper).high);
                ++samples;
            }
        }
    }
    std::cout << "boundary " << samples << '\n';
}

}  // namespace

int main() {
    assert(std::fegetround() == FE_TONEAREST);
    for (bool upper : {true, false}) {
        certify_piece(upper, Piece::local, Interval(-local_radius.high, local_radius.high));
        certify_piece(upper, Piece::horizontal, Interval(upper ? -5 : -3, -local_radius.low));
        certify_piece(upper, Piece::horizontal, Interval(local_radius.low, upper ? 12 : 3));
        certify_piece(upper, Piece::split_vertical, Interval(0, 1));
        if (!upper) {
            certify_piece(upper, Piece::left_vertical, Interval(0, 1));
            certify_piece(upper, Piece::origin, Interval(0, 2));
        }
    }
    certify_boundaries();
}

