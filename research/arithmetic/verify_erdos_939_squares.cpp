#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

using i64 = int64_t;
using u64 = uint64_t;
using u128 = __uint128_t;

u128 power(u64 base, unsigned exponent) {
    u128 result = 1;
    for (unsigned i = 0; i < exponent; ++i) result *= base;
    return result;
}

// A floating approximation gives a starting point. The returned root is exact.
u64 floor_root(u128 n, unsigned exponent) {
    u64 root = std::pow((long double)n, 1.0L / exponent);
    while (power(root, exponent) > n) --root;
    while (power(root + 1, exponent) <= n) ++root;
    return root;
}

std::vector<u64> odd_squarefull(u64 limit) {
    const u64 max_conductor = floor_root(limit, 3);
    // Odd squarefull integers have a unique representation d^3 x^2
    // with d odd squarefree and x odd.
    std::vector<bool> squarefree(max_conductor + 1, true);
    for (u64 p = 3; p * p <= max_conductor; p += 2) {
        if (!squarefree[p]) continue;
        for (u64 multiple = p * p; multiple <= max_conductor; multiple += 2 * p * p) {
            squarefree[multiple] = false;
        }
    }
    std::vector<u64> values;
    for (u64 d = 1; d <= max_conductor; d += 2) {
        if (!squarefree[d]) continue;
        u64 d_cubed = d * d * d;
        for (u64 x = 1; (u128)d_cubed * x * x <= limit; x += 2) {
            values.push_back(d_cubed * x * x);
        }
    }
    std::sort(values.begin(), values.end());
    return values;
}

bool squarefull_by_trial_division(u128 n) {
    for (u64 p = 2; (u128)p * p <= n; p += (p == 2 ? 1 : 2)) {
        if (n % p != 0) continue;
        const u128 p_squared = (u128)p * p;
        if (n % p_squared != 0) return false;
        do {
            n /= p;
        } while (n % p == 0);
    }
    return n == 1;
}

// Euclid parameters for a primitive Pythagorean triple. The input is
// the fourth root of the upper bound for c.
void all_square_terms(u64 limit) {
    const u128 hypotenuse_limit = (u128)limit * limit;
    const u64 linear_limit = floor_root(2 * hypotenuse_limit, 2);
    const auto odds = odd_squarefull(linear_limit);
    const std::unordered_set<u64> odd_lookup(odds.begin(), odds.end());
    std::vector<std::vector<u64>> evens(225);
    for (const u64 odd : odds) {
        for (u64 even = 2 * odd; even < limit; even *= 2) {
            evens[even % 225].push_back(even);
        }
    }
    for (auto& bucket : evens) std::sort(bucket.begin(), bucket.end());

    // Modulo 9 and 25, no parameter, odd leg factor, or hypotenuse
    // may have p-adic valuation 1. Exclude nonprimitive pairs as well.
    std::vector<std::vector<int>> compatible(225);
    for (int m = 0; m < 225; ++m) {
        for (int n = 0; n < 225; ++n) {
            bool passes = true;
            for (const int p : {3, 5}) {
                if (m % p == 0 && n % p == 0) {
                    passes = false;
                    break;
                }
                for (const int value : {m, n, m+n, m-n, m*m+n*n}) {
                    if (value % p == 0 && value % (p*p) != 0) {
                        passes = false;
                        break;
                    }
                }
            }
            if (passes) compatible[m].push_back(n);
        }
    }

    const std::size_t odd_count = std::lower_bound(odds.begin(), odds.end(), limit) - odds.begin();
    std::size_t even_count = 0;
    for (const auto& bucket : evens) even_count += bucket.size();
    std::cerr << "ODDS=" << odd_count << " EVENS=" << even_count << " PAIRS=" << odd_count * even_count << "\n";
    std::size_t residue_passes = 0, linear_passes = 0, hits = 0;
    for (const u64 odd : odds) {
        if (odd >= limit) break;
        for (const int bucket : compatible[odd % 225]) {
            for (const u64 even : evens[bucket]) {
                const u128 hypotenuse = (u128)odd * odd + (u128)even * even;
                if (hypotenuse > hypotenuse_limit) break;
                ++residue_passes;
                if (!odd_lookup.contains(odd + even)) continue;
                if (!odd_lookup.contains(odd > even ? odd - even : even - odd)) continue;
                if (std::gcd(odd, even) != 1) continue;
                ++linear_passes;
                if (!squarefull_by_trial_division(hypotenuse)) continue;
                std::cout << "HIT " << odd << " " << even << "\n";
                ++hits;
            }
        }
    }
    std::cerr << "RESIDUE_AND_HEIGHT_PASSES=" << residue_passes << " LINEAR_PASSES=" << linear_passes << " HITS=" << hits << "\n";
}

std::vector<u64> odd_fourfull(u64 limit) {
    const u64 max_coeff_base = floor_root(limit, 5);
    std::vector<bool> squarefree(max_coeff_base + 1, true);
    for (u64 p = 3; p * p <= max_coeff_base; p += 2) {
        if (!squarefree[p]) continue;
        for (u64 multiple = p * p; multiple <= max_coeff_base; multiple += 2 * p * p) {
            squarefree[multiple] = false;
        }
    }
    // The exponents 5, 6, 7 represent the three nonzero classes mod 4.
    std::vector<u64> values;
    for (u64 e = 1; power(e, 5) <= limit; e += 2) {
        if (!squarefree[e]) continue;
        const u64 e5 = power(e, 5);
        for (u64 f = 1; power(f, 6) <= limit / e5; f += 2) {
            if (!squarefree[f] || std::gcd(e, f) != 1) continue;
            const u64 ef = e5 * power(f, 6);
            for (u64 g = 1; power(g, 7) <= limit / ef; g += 2) {
                if (!squarefree[g] || std::gcd(e*f, g) != 1) continue;
                const u64 coeff = ef * power(g, 7);
                const u64 max_d4 = limit / coeff;
                for (u64 d = 1; power(d, 4) <= max_d4; d += 2) {
                    values.push_back(coeff * power(d, 4));
                }
            }
        }
    }
    std::sort(values.begin(), values.end());
    return values;
}

bool local_squarefull(u64 n, u64 max_prime) {
    for (const u64 p : {2, 3, 5, 7, 11, 13, 17}) {
        if (p > max_prime) break;
        if (n % p == 0 && n % (p*p) != 0) return false;
    }
    return true;
}

// Search primitive solutions with a=x^2 and c=z^2. The input bounds z.
void square_c_and_summand(u64 limit) {
    const auto squarefull = odd_squarefull(limit);
    const auto fourfull = odd_fourfull(2 * limit);
    const std::size_t fourfull_even_max = std::lower_bound(fourfull.begin(), fourfull.end(), (limit + 3) / 4) - fourfull.begin();
    std::vector<u64> fourfull_times_4;
    for (std::size_t i = 0; i < fourfull_even_max; ++i) {
        for (u64 even = 4 * fourfull[i]; even < limit; even *= 2) {
            fourfull_times_4.push_back(even);
        }
    }
    std::sort(fourfull_times_4.begin(), fourfull_times_4.end());
    std::cerr << "ODD_SQUAREFULL=" << squarefull.size() << " ODD_FOURFULL=" << fourfull.size() << " EVEN_PARAMS=" << fourfull_times_4.size() << "\n";

    auto is_squarefull = [&](u64 n) {
        if (!local_squarefull(n, 11)) return false;
        n >>= __builtin_ctzll(n);
        return std::binary_search(squarefull.begin(), squarefull.end(), n);
    };
    std::size_t mixed_coprime = 0, odd_coprime = 0;
    // Odd x and z: use (z+x)/2 and (z-x)/2, with product b/4.
    for (const u64 odd : fourfull) {
        if (odd >= limit) break;
        for (const u64 even : fourfull_times_4) {
            if (odd + even > limit) break;
            const u64 smaller = odd > even ? odd - even : even - odd;
            if (!local_squarefull(smaller, 11)) continue;
            if (!is_squarefull(odd+even) || !is_squarefull(smaller)) continue;
            if (std::gcd(odd, even) != 1) continue;
            std::cout << "EVEN_GAP " << odd << " " << even << "\n";
            ++odd_coprime;
        }
    }
    // Opposite parity: use z+x and z-x, with product b.
    for (std::size_t i = 0; i < fourfull.size(); ++i) {
        const u64 lower = fourfull[i];
        if (lower >= limit) break;
        for (std::size_t j = i+1; j < fourfull.size(); ++j) {
            const u64 upper = fourfull[j];
            if (upper + lower > 2 * limit) break;
            const u64 smaller = (upper - lower) / 2;
            const u64 larger = (upper + lower) / 2;
            if (!local_squarefull(smaller, 11)) continue;
            if (!is_squarefull(larger) || !is_squarefull(smaller)) continue;
            if (std::gcd(upper, lower) != 1) continue;
            std::cout << "ODD_GAP " << upper << " " << lower << "\n";
            ++mixed_coprime;
        }
    }
    std::cerr << "ODD_SQUARE_HITS=" << odd_coprime << " MIXED_SQUARE_HITS=" << mixed_coprime << "\n";
}

// When a and b are squares, coprimality and 4-fullness of c imply
// opposite parity and a primitive Gaussian integer with odd norm.
struct Gaussian {
    i64 re;
    i64 im;
};

Gaussian operator*(Gaussian a, Gaussian b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

u64 mod_power(u64 base, u64 exponent, u64 prime) {
    u64 result = 1;
    while (exponent != 0) {
        if (exponent % 2 == 1) {
            result = (u128)result * base % prime;
        }
        base = (u128)base * base % prime;
        exponent /= 2;
    }
    return result;
}

Gaussian prime_above(u64 prime) {
    u64 nonresidue = 2;
    while (mod_power(nonresidue, (prime - 1) / 2, prime) == 1) {
        ++nonresidue;
    }
    u64 previous = prime;
    u64 current = mod_power(nonresidue, (prime - 1) / 4, prime);
    while (current * current > prime) {
        const u64 next = previous % current;
        previous = current;
        current = next;
    }
    const u64 other_sq = prime - current * current;
    const u64 other = floor_root(other_sq, 2);
    if (other * other != other_sq) throw std::runtime_error("Cornacchia failed");
    return {(i64)current, (i64)other};
}

class SquareSummands {
public:
    explicit SquareSummands(u64 root_limit)
        : norm_limit((u128)root_limit * root_limit),
          odd_values(odd_squarefull(root_limit)) {
        const u64 max_prime = floor_root(root_limit, 2);
        std::vector<bool> composite(max_prime + 1);
        for (u64 p = 2; p * p <= max_prime; ++p) {
            if (composite[p]) continue;
            for (u64 multiple = p * p; multiple <= max_prime; multiple += p) {
                composite[multiple] = true;
            }
        }
        for (u64 p = 5; p <= max_prime; p += 4) {
            if (!composite[p]) {
                const Gaussian g = prime_above(p);
                const Gaussian squared = g * g;
                split_primes.push_back({p, g, squared * squared, (u128)p*p*p*p});
            }
        }
    }

    void run() {
        std::cerr << "ODD_SQUAREFULL=" << odd_values.size() << " SPLIT_PRIMES=" << split_primes.size() << "\n";
        extend({1, 0}, 1, 0, true);
        std::cerr << "PRIMITIVE_REPRESENTATIONS=" << representations << " LOCAL_PASSES=" << local_passes
                  << " HITS=" << hits << "\n";
    }

private:
    struct SplitPrime {
        u64 prime;
        Gaussian choice;
        Gaussian fourth_power;
        u128 norm_fourth;
    };

    const u128 norm_limit;
    const std::vector<u64> odd_values;
    std::vector<SplitPrime> split_primes;
    std::size_t representations = 0, local_passes = 0, hits = 0;

    void inspect(Gaussian value) {
        u64 x = std::abs(value.re);
        u64 y = std::abs(value.im);
        if (x == 0 || y == 0) throw std::runtime_error("invalid primitive representation");
        ++representations;
        if (!local_squarefull(x, 17) || !local_squarefull(y, 17)) return;
        ++local_passes;
        const u64 odd_x = x >> __builtin_ctzll(x);
        const u64 odd_y = y >> __builtin_ctzll(y);
        if (!std::binary_search(odd_values.begin(), odd_values.end(), odd_x)) return;
        if (!std::binary_search(odd_values.begin(), odd_values.end(), odd_y)) return;
        std::cout << "HIT " << x << " " << y << "\n";
        ++hits;
    }

    // Choose one prime above each selected rational prime. Conjugate
    // globally so the first selection is fixed. Units only swap signs
    // and order of the coordinates.
    void extend(Gaussian current, u128 norm, std::size_t next_prime, bool first) {
        const u128 room = norm_limit / norm;
        for (std::size_t i = next_prime; i < split_primes.size(); ++i) {
            const SplitPrime& factor = split_primes[i];
            if (factor.norm_fourth > room) break;
            Gaussian value = factor.fourth_power;
            u128 factor_norm = factor.norm_fourth;
            while (factor_norm <= room) {
                const Gaussian product = current * value;
                const u128 next_norm = norm * factor_norm;
                inspect(product);
                extend(product, next_norm, i+1, false);
                if (!first) {
                    const Gaussian other = current * Gaussian{value.re, -value.im};
                    inspect(other);
                    extend(other, next_norm, i+1, false);
                }
                if (factor_norm > room / factor.prime) break;
                value = value * factor.choice;
                factor_norm *= factor.prime;
            }
        }
    }
};


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: 939_squares all-square|square-c|square-summands limit\n";
        return 2;
    }
    const std::string_view mode = argv[1];
    const std::string_view input = argv[2];
    if (input.empty() || input.find_first_not_of("0123456789") != std::string_view::npos) return 2;

    const u64 limit = std::stoull(argv[2]);
    if (limit == 0) return 2;
    if (mode == "all-square" && limit <= 10'000'000'000ULL) {
        all_square_terms(limit);
    } else if (mode == "square-c" && limit <= 1'000'000'000'000'000ULL) {
        square_c_and_summand(limit);
    } else if (mode == "square-summands" && limit <= 100'000'000'000'000'000ULL) {
        SquareSummands(limit).run();
    } else {
        std::cerr << "mode or limit out of range\n";
        return 2;
    }
}
