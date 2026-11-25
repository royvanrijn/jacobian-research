// Exact characteristic-zero squarefreeness certificate for the finite constant case.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

using std::int64_t;
constexpr int first_gap = 6;
constexpr int last_gap = 29;
constexpr int last_order = 9999;
constexpr int large_prime = 1000000007;

struct PrimeField {
    int p;
    int product(int a, int b) const { return int64_t(a) * b % p; }

    int inverse(int a) const {
        int result = 1;
        for (int exponent = p - 2; exponent > 0; exponent /= 2) {
            if (exponent % 2 != 0) result = product(result, a);
            a = product(a, a);
        }
        return result;
    }

    std::vector<int> reciprocal(int order, int gap) const {
        std::vector<int> polynomial(order + 1);
        polynomial[order] = 1;
        int coefficient = 1;
        int argument = 0;
        for (int k = 1; k <= order; ++k) {
            coefficient = product(product(coefficient, order - k + 1), inverse(k));
            for (int j = 0; j < gap; ++j) {
                coefficient = product(coefficient, ++argument % p);
            }
            polynomial[order - k] = coefficient;
        }
        return polynomial;
    }

    bool gcd_degree_at_most_one(std::vector<int> a) const {
        if (a.size() <= 2) return true;
        std::vector<int> b(a.size() - 1);
        for (int i = 1; i < int(a.size()); ++i) b[i - 1] = product(i, a[i]);
        while (b.size() > 2) {
            const int degree = b.size() - 1;
            const int inverse_lead = inverse(b.back());
            for (int i = int(a.size()) - 1; i >= degree; --i) {
                if (!a[i]) continue;
                const int scale = product(a[i], inverse_lead);
                for (int j = 0; j <= degree; ++j) {
                    int& entry = a[i - degree + j];
                    entry -= product(scale, b[j]);
                    if (entry < 0) entry += p;
                }
            }
            while (!a.empty() && a.back() == 0) a.pop_back();
            if (a.empty()) return false;
            a.swap(b);
        }
        return true;
    }
};

int main() {
    for (int d = 2; int64_t(d) * d <= large_prime; ++d) {
        if (large_prime % d == 0) return 1;
    }
    if (large_prime <= last_order * last_gap) return 1;
    std::vector<bool> composite(last_order + 1);
    std::vector<int> primes;
    for (int p = 2; p <= last_order; ++p) {
        if (composite[p]) continue;
        primes.push_back(p);
        for (int j = p * p; j <= last_order; j += p) composite[j] = true;
    }

    int total = 0;
    int shifted = 0;
    int full = 0;
    int unresolved = 0;
    int max_small_order = 0;
    int max_full_order = 0;
    for (int gap = first_gap; gap <= last_gap; ++gap) {
        for (int order = 2; order <= last_order; ++order) {
            if ((int64_t(1) << gap) >= 4 * int64_t(order) * order) continue;
            ++total;
            bool certified = false;
            auto position = std::upper_bound(primes.begin(), primes.end(), order);
            while (position != primes.begin()) {
                const int p = *--position;
                const int small_order = order - p;
                if (p <= gap * std::max(1, small_order)) break;
                const PrimeField field{p};
                if (field.gcd_degree_at_most_one(field.reciprocal(small_order, gap))) {
                    max_small_order = std::max(max_small_order, small_order);
                    ++shifted;
                    certified = true;
                    break;
                }
            }
            if (certified) continue;
            const PrimeField field{large_prime};
            if (field.gcd_degree_at_most_one(field.reciprocal(order, gap))) {
                max_full_order = std::max(max_full_order, order);
                ++full;
            } else {
                ++unresolved;
                std::cerr << "Unresolved (order, gap): " << order << ' ' << gap << '\n';
            }
        }
    }
    std::cout << "Total pairs: " << total << '\n'
              << "Shifted-prime pairs: " << shifted << '\n'
              << "Largest shifted degree: " << max_small_order << '\n'
              << "Full-modular pairs: " << full << '\n'
              << "Largest full degree: " << max_full_order << '\n'
              << "Unresolved pairs: " << unresolved << '\n';
    return total == 202021 && total == full + shifted && unresolved == 0 ? 0 : 1;
}
