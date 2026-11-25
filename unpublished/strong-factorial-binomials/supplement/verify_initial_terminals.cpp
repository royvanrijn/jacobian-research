// Exact finite certificate for the initial ordinary binomial terminals.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using std::int64_t;

constexpr int modulus = 1000000007;
constexpr int first_terminal = 3;
constexpr int last_terminal = 24;
constexpr int max_gap = 42;
constexpr int max_argument = last_terminal * (max_gap * (last_terminal - 1) - 1);

int multiply_mod(int a, int b) {
    return int64_t(a) * b % modulus;
}

int inverse_mod(int value) {
    int result = 1;
    for (int power = modulus - 2; power > 0; power /= 2) {
        if (power % 2 != 0) result = multiply_mod(result, value);
        value = multiply_mod(value, value);
    }
    return result;
}

struct Factorials {
    std::array<int, max_argument + 1> direct;
    std::array<int, max_argument + 1> inverse;

    Factorials() {
        direct[0] = 1;
        for (int j = 1; j <= max_argument; ++j) {
            direct[j] = multiply_mod(j, direct[j - 1]);
        }
        inverse[max_argument] = inverse_mod(direct[max_argument]);
        for (int j = max_argument; j >= 1; --j) {
            inverse[j - 1] = multiply_mod(j, inverse[j]);
        }
    }
};

std::vector<int> reciprocal(int m, int a, int d, const Factorials& factorials) {
    std::vector<int> values(m + 1);
    for (int k = 0; k <= m; ++k) {
        const int binomial = multiply_mod(
            factorials.direct[m],
            multiply_mod(factorials.inverse[k], factorials.inverse[m - k]));
        const int factorial_ratio = multiply_mod(
            factorials.direct[a * m + d * k], factorials.inverse[a * m]);
        values[m - k] = multiply_mod(binomial, factorial_ratio);
    }
    return values;
}

bool coprime_mod(std::vector<int> dividend, std::vector<int> divisor) {
    while (divisor.size() > 1) {
        const int inverse_lead = inverse_mod(divisor.back());
        const int degree = divisor.size() - 1;
        for (int i = int(dividend.size()) - 1; i >= degree; --i) {
            if (dividend[i] == 0) continue;
            const int scale = multiply_mod(dividend[i], inverse_lead);
            for (int j = 0; j <= degree; ++j) {
                int& entry = dividend[i - degree + j];
                entry -= multiply_mod(scale, divisor[j]);
                if (entry < 0) entry += modulus;
            }
        }
        while (!dividend.empty() && dividend.back() == 0) dividend.pop_back();
        if (dividend.empty()) return false;
        dividend.swap(divisor);
    }
    return divisor.size() == 1;
}

int main() {
    for (int p = 2; int64_t(p) * p <= modulus; ++p) {
        if (modulus % p == 0) return 1;
    }
    if (max_argument >= modulus) return 1;
    const Factorials factorials;
    int total = 0;
    int unresolved = 0;
    for (int N = first_terminal; N <= last_terminal; ++N) {
        int row = 0;
        for (int d = 1; d <= max_gap; ++d) {
            for (int a = 1; a < d * (N - 2); ++a) {
                if (d * d >= 20 * (a + d)) continue;
                ++row;
                if (!coprime_mod(reciprocal(N, a, d, factorials),
                                 reciprocal(N - 1, a, d, factorials))) {
                    ++unresolved;
                    std::cerr << "Unresolved (N,a,d): " << N << ' ' << a << ' ' << d << '\n';
                }
            }
        }
        total += row;
        std::cout << "Terminal " << N << ": " << row << '\n';
    }
    std::cout << "Total triples: " << total << '\n'
              << "Unresolved triples: " << unresolved << '\n';
    return total == 217952 && unresolved == 0 ? 0 : 1;
}
