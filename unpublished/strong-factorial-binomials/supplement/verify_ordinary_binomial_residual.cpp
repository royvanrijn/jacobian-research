// Exact prime and Newton polygon certificate for 25 <= N < 3150.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

using std::int64_t;
constexpr int first_terminal = 25;
constexpr int last_terminal = 3149;

// Bounds (T33) and (T40) on d, independent of a and the factor degree.
int gap_cap(int terminal) {
    if (terminal <= 50) return 35;
    if (terminal <= 100) return 56;
    if (terminal <= 300) return 76;
    if (terminal <= 1000) return 110;
    return 140;
}

// Strict Gamma bound 6a < (5N-11)d.
int lower_exponent_cap(int terminal, int gap) {
    return ((5 * terminal - 11) * gap - 1) / 6;
}

int floor_sqrt(int value) {
    int root = static_cast<int>(std::sqrt(value));
    while (int64_t(root + 1) * (root + 1) <= value) ++root;
    while (int64_t(root) * root > value) --root;
    return root;
}

class PrimeFactors {
public:
    explicit PrimeFactors(int limit) : smallest_odd_factor(limit / 2 + 1) {
        const int root = floor_sqrt(limit);
        std::vector<int> primes;
        std::vector<bool> composite(root + 1);
        for (int p = 3; p <= root; p += 2) {
            if (composite[p]) continue;
            primes.push_back(p);
            for (int i = p * p; i <= root; i += 2 * p) composite[i] = true;
        }

        constexpr int segment = 1 << 20;
        #pragma omp parallel for schedule(dynamic)
        for (int start = 1; start <= limit; start += segment) {
            const int end = std::min(start + segment - 1, limit);
            for (int p : primes) {
                if (p * p > end) break;
                int multiple = std::max(p * p, (start + p - 1) / p * p);
                if ((multiple & 1) == 0) multiple += p;
                for (; multiple <= end; multiple += 2 * p) {
                    auto& entry = smallest_odd_factor[multiple / 2];
                    if (!entry) entry = p;
                }
            }
        }
    }

    bool is_odd_prime(int number) const {
        return number > 1 && (number & 1) && !smallest_odd_factor[number / 2];
    }

    bool has_factor_above(int number, int bound) const {
        number >>= __builtin_ctz(static_cast<unsigned>(number));
        while (number > bound) {
            const int p = smallest_odd_factor[number / 2];
            if (!p) return true;
            if (p > bound) return true;
            number /= p;
        }
        return false;
    }

    bool block_has_factor_above(int start, int width, int bound) const {
        for (int j = 1; j <= width; ++j) {
            if (has_factor_above(start + j, bound)) return true;
        }
        return false;
    }

private:
    std::vector<uint16_t> smallest_odd_factor;
};

int prime_cap(int terminal, int exponent, int gap) {
    const int unrestricted = gap * terminal;
    if (terminal < 275 && exponent < 2 * gap) return unrestricted;
    const int shape_bound = std::max(exponent + 2 * gap - 1,
                                     floor_sqrt((exponent + gap) * terminal));
    return std::min(unrestricted, shape_bound);
}


struct Edge {
    int rise, run, width;
};

std::vector<Edge> sparse_polygon(int m, int a, int d, int p) {
    std::vector<std::pair<int,int>> hull;
    auto add = [&](int index, int ordinate) {
        if (!hull.empty() && hull.back().first == index) return;
        while (hull.size() >= 2) {
            const auto [i, y] = hull[hull.size() - 2];
            const auto [j, z] = hull.back();
            if (int64_t(z - y) * (index - j) < int64_t(ordinate - z) * (j - i)) break;
            hull.pop_back();
        }
        hull.emplace_back(index, ordinate);
    };

    // p>aN, p^2>(a+d)N. The ordinates are floor((am+dk)/p).
    add(0, 0);
    const int jumps = (a + d) * m / p;
    for (int j = 1; j <= jumps; ++j) {
        add((j * p - a * m - 1) / d, j - 1);
    }
    add(m, jumps);

    std::vector<Edge> edges;
    for (size_t j = 1; j < hull.size(); ++j) {
        const int w = hull[j].first - hull[j - 1].first;
        const int r = hull[j].second - hull[j - 1].second;
        const int g = std::gcd(w, r);
        edges.push_back({r / g, w / g, w});
    }
    return edges;
}

bool possible_at_prime(std::vector<int>& degrees, int n, int a, int d, int p) {
    const auto lower = sparse_polygon(n - 1, a, d, p);
    const auto upper = sparse_polygon(n, a, d, p);
    const int limit = degrees.back();
    std::vector<char> allocations(limit + 1);
    allocations[0] = true;
    for (const Edge& left : lower) {
        for (const Edge& right : upper) {
            if (left.rise != right.rise || left.run != right.run) continue;
            const int width = std::min(left.width, right.width);
            for (int old = limit; old >= 0; --old) {
                if (!allocations[old]) continue;
                for (int added = left.run; added <= width && old + added <= limit; added += left.run) {
                    allocations[old + added] = true;
                }
            }
        }
    }
    std::erase_if(degrees, [&](int h) { return !allocations[h]; });
    return !degrees.empty();
}

enum class Exclusion { block, polygon, unresolved };

Exclusion exclude_remaining_support(const PrimeFactors& factors, int n, int a, int d) {
    const int x = a * n;
    const int m = (a + d) * n;
    const int shape_bound = a + 2 * d - 1;
    const int bound = std::max(shape_bound, floor_sqrt(m));
    std::vector<int> degrees;
    for (int h = 2; h < n - 1; h += 2) {
        const int width = d * h;
        const bool t10 = n >= 275 || (x + 1 > width + 13 && 262 * int64_t(x + 1) > 279 * width);
        if (t10 && 2 * width >= shape_bound && int64_t(2 * width) * (2 * width) >= m) continue;
        const int prime_bound = std::min(d * n, t10 ? bound : std::max(2 * width, bound));
        if (factors.block_has_factor_above(x, width, prime_bound)
            || factors.block_has_factor_above(a * (n - 1), width, prime_bound)) {
            continue;
        }
        degrees.push_back(h);
    }

    if (degrees.empty()) return Exclusion::block;

    for (int p = (x + 1) | 1; p <= m; p += 2) {
        if (factors.is_odd_prime(p) && !possible_at_prime(degrees, n, a, d, p)) {
            return Exclusion::polygon;
        }
    }
    return Exclusion::unresolved;
}


constexpr int modulus = 1000000007;

int product_mod(int a, int b) {
    return int64_t(a) * b % modulus;
}

int inverse_mod(int value) {
    int inverse = 1;
    for (int exponent = modulus - 2; exponent > 0; exponent /= 2) {
        if (exponent % 2 != 0) inverse = product_mod(inverse, value);
        value = product_mod(value, value);
    }
    return inverse;
}

std::vector<int> reciprocal_mod(int m, int a, int d, const std::vector<int>& inverses) {
    std::vector<int> polynomial(m + 1);
    polynomial[m] = 1;
    int coefficient = 1;
    int factorial_argument = a * m;
    for (int k = 1; k <= m; ++k) {
        coefficient = product_mod(product_mod(coefficient, m - k + 1), inverses[k]);
        for (int i = 0; i < d; ++i) {
            coefficient = product_mod(coefficient, ++factorial_argument % modulus);
        }
        polynomial[m - k] = coefficient;
    }
    return polynomial;
}

bool coprime_mod(std::vector<int> dividend, std::vector<int> divisor) {
    while (divisor.size() > 1) {
        const int inverse_lead = inverse_mod(divisor.back());
        const int degree = divisor.size() - 1;
        for (int i = int(dividend.size()) - 1; i >= degree; --i) {
            if (dividend[i] == 0) continue;
            const int scale = product_mod(dividend[i], inverse_lead);
            const int shift = i - degree;
            for (int j = 0; j < degree; ++j) {
                int value = dividend[shift + j] - product_mod(scale, divisor[j]);
                if (value < 0) value += modulus;
                dividend[shift + j] = value;
            }
            dividend[i] = 0;
        }
        while (!dividend.empty() && dividend.back() == 0) dividend.pop_back();
        if (dividend.empty()) return false;
        dividend.swap(divisor);
    }
    return divisor.size() == 1;
}


int main(int argc, char** argv) {
    const bool list_pairs = argc == 2 && std::string_view(argv[1]) == "--list-gcd-pairs";
    if (argc != 1 && !list_pairs) {
        std::cerr << "Usage: " << argv[0] << " [--list-gcd-pairs]\n";
        return 2;
    }
    for (int p = 2; int64_t(p) * p <= modulus; ++p) {
        if (modulus % p == 0) return 1;
    }
    if (last_terminal >= modulus) return 1;
    const int max_gap = gap_cap(last_terminal);
    const int limit = last_terminal * (lower_exponent_cap(last_terminal, max_gap) + max_gap);
    const PrimeFactors factors(limit);
    std::cout << "Factor table upper endpoint: " << limit << '\n';
    std::vector<std::vector<std::pair<int,int>>> survivors(last_terminal + 1);
    int64_t coarse_pairs = 0;
    int64_t refined_pairs = 0, block = 0, polygon = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:coarse_pairs,refined_pairs,block,polygon)
    for (int n = first_terminal; n <= last_terminal; ++n) {
        const int gap_upper = gap_cap(n);
        const int a_upper = lower_exponent_cap(n, gap_upper);
        for (int a = 5; a <= a_upper; ++a) {
            int gap_min = std::max(2, 6 * a / (5 * n - 11) + 1);
            int gap_max = std::min(gap_upper, 3 * a - 1);
            if (n >= 2200) {
                gap_min = std::max(gap_min, n - a + 1);
                gap_max = std::min(gap_max, (3 * a - 1) / 97);
            }
            if (gap_min > gap_max) continue;
            const int upper_start = a * n;
            const int lower_start = a * (n - 1);
            // The unrestricted branch can only increase the cap.
            const int uniform_cap = prime_cap(n, a, gap_max);
            if (factors.block_has_factor_above(upper_start, 2 * gap_min, uniform_cap)
                || factors.block_has_factor_above(lower_start, 2 * gap_min, uniform_cap)) {
                continue;
            }
            ++coarse_pairs;
            for (int d = gap_min; d <= gap_max; ++d) {
                const int b = a + d;
                if (a >= d * (n - 2) || d * d >= 19 * b) continue;
                const int bound = prime_cap(n, a, d);
                if (factors.block_has_factor_above(upper_start, 2 * d, bound)
                    || factors.block_has_factor_above(lower_start, 2 * d, bound)) {
                    continue;
                }
                ++refined_pairs;
                switch (exclude_remaining_support(factors, n, a, d)) {
                    case Exclusion::block: ++block; break;
                    case Exclusion::polygon: ++polygon; break;
                    case Exclusion::unresolved: survivors[n].emplace_back(a, d); break;
                }
            }
        }
    }
    std::cout << "Groups (N,a) passing the initial prime cap: " << coarse_pairs << '\n'
              << "Supports reaching degree tests: " << refined_pairs << '\n'
              << "Excluded by degree and prime blocks: " << block << '\n'
              << "Excluded by Newton allocations: " << polygon << '\n';
    if (list_pairs) {
        for (int n = first_terminal; n <= last_terminal; ++n) {
            for (const auto& [a, d] : survivors[n]) {
                std::cout << "Replay pair " << n << ' ' << a << ' ' << d << '\n';
            }
        }
    }
    std::vector<int> inverses(last_terminal + 1);
    for (int k = 1; k <= last_terminal; ++k) inverses[k] = inverse_mod(k);
    int64_t modular = 0, unresolved = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:modular,unresolved)
    for (int n = first_terminal; n <= last_terminal; ++n) {
        for (const auto& [a, d] : survivors[n]) {
            if (coprime_mod(reciprocal_mod(n, a, d, inverses),
                            reciprocal_mod(n - 1, a, d, inverses))) {
                ++modular;
            } else {
                #pragma omp critical
                std::cerr << "Unresolved (N,a,d): " << n << ' ' << a << ' ' << d << '\n';
                ++unresolved;
            }
        }
    }
    std::cout << "Excluded by modular coprimality: " << modular << '\n'
              << "Unresolved supports: " << unresolved << '\n';
    return unresolved == 0 && block + polygon + modular == refined_pairs ? 0 : 1;
}
