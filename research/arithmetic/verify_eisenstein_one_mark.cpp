#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

constexpr int modulus = 998244353;
constexpr int generator = 3;

int multiply(int a, int b) {
    return static_cast<std::int64_t>(a) * b % modulus;
}

int power(int base, int exponent) {
    int result = 1;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = multiply(result, base);
        }
        base = multiply(base, base);
        exponent /= 2;
    }
    return result;
}

void transform(std::vector<int>& values, bool inverse) {
    const int length = static_cast<int>(values.size());
    for (int i = 1, j = 0; i < length; ++i) {
        int bit = length / 2;
        while (j & bit) {
            j ^= bit;
            bit /= 2;
        }
        j ^= bit;
        if (i < j) {
            std::swap(values[i], values[j]);
        }
    }

    const int base = inverse ? power(generator, modulus - 2) : generator;
    for (int block = 2; block <= length; block *= 2) {
        const int step = power(base, (modulus - 1) / block);
        for (int offset = 0; offset < length; offset += block) {
            int root = 1;
            for (int i = offset; i < offset + block / 2; ++i) {
                const int j = i + block / 2;
                const int left = values[i];
                const int right = multiply(values[j], root);
                const int sum = left + right;
                const int difference = left - right;
                values[i] = sum < modulus ? sum : sum - modulus;
                values[j] = difference >= 0 ? difference : difference + modulus;
                root = multiply(root, step);
            }
        }
    }
    if (inverse) {
        const int inverse_length = power(length, modulus - 2);
        for (int& value : values) {
            value = multiply(value, inverse_length);
        }
    }
}

void square_convolution(std::vector<int>& values) {
    const auto length = std::bit_ceil(2 * values.size() - 1);
    assert(length <= (1U << 23));
    values.resize(length);
    transform(values, false);
    for (int& value : values) {
        value = multiply(value, value);
    }
    transform(values, true);
}

std::int64_t cube(int value) {
    return static_cast<std::int64_t>(value) * value * value;
}

int cube_floor(int value) {
    int root = 0;
    while (cube(root + 1) <= value) {
        ++root;
    }
    return root;
}

std::array<int, 2> direct_counts(int target, int cutoff, const std::vector<int>& least) {
    const auto kind = [&](int value) {
        const int prime = least[value];
        if (prime <= cutoff) {
            return 0;
        }
        const int cofactor = value / prime;
        if (cofactor == 1) {
            return 1;
        }
        return prime % 6 == 1 && cofactor % 6 == 1 && cofactor != prime &&
               least[cofactor] == cofactor ? 2 : 0;
    };

    std::array<int, 2> counts{};
    for (int a = 7; a < target / 2; a += 6) {
        const int b = target - a;
        const int type = kind(a);
        if (type != 0 && type == kind(b) && std::gcd(a, b) == 1) {
            counts[type - 1] += 2;
        }
    }
    return counts;
}

}  // namespace

int main(int argc, char** argv) {
    const int limit = argc == 2 ? std::stoi(argv[1]) : 20'000'000;
    assert(14 <= limit && limit <= 20'000'000);
    for (int divisor = 2; divisor * divisor <= modulus; ++divisor) {
        assert(modulus % divisor != 0);
    }
    assert(power(generator, modulus - 1) == 1);
    for (const int prime : {2, 7, 17}) {
        assert(power(generator, (modulus - 1) / prime) != 1);
    }

    std::vector<int> least(limit + 1);
    std::iota(least.begin(), least.end(), 0);
    for (int prime = 2; prime * prime <= limit; ++prime) {
        if (least[prime] == prime) {
            for (int multiple = prime * prime; multiple <= limit; multiple += prime) {
                if (least[multiple] == multiple) {
                    least[multiple] = prime;
                }
            }
        }
    }

    // The index i denotes the value 6i+1. A sum at index n is 6n+2.
    std::vector<int> primes;
    std::vector<int> prime_pairs((limit - 1) / 6 + 1);
    for (int value = 7; value <= limit; value += 6) {
        if (least[value] == value) {
            primes.push_back(value);
            prime_pairs[(value - 1) / 6] = 1;
        }
    }
    square_convolution(prime_pairs);

    int targets = 0;
    int fallbacks = 0;
    int crosschecks = 0;
    int nonpositive_weighted = 0;
    int nonpositive_unweighted = 0;
    int worst_target = 0;
    int worst_prime_pairs = 1;
    int worst_semiprime_upper = 0;
    for (int low = 8; low <= limit; low *= 2) {
        const int high = std::min(2 * low - 1, limit);
        const int cutoff = cube_floor(low);
        std::vector<int> semiprime_pairs((high - 1) / 6 + 1);
        for (auto p = std::upper_bound(primes.begin(), primes.end(), cutoff);
             p != primes.end() && static_cast<std::int64_t>(*p) * *p <= high; ++p) {
            for (auto q = p + 1; q != primes.end() && *q <= high / *p; ++q) {
                semiprime_pairs[(*p * *q - 1) / 6] = 1;
            }
        }
        square_convolution(semiprime_pairs);

        const int first = low + (8 - low % 6) % 6;
        for (int target = first; target <= high; target += 6) {
            ++targets;
            const int exact_cutoff = cube_floor(target);
            const int index = (target - 2) / 6;
            int prime_count = prime_pairs[index];
            for (int prime : primes) {
                if (prime > exact_cutoff) {
                    break;
                }
                const int partner = target - prime;
                if (least[partner] == partner) {
                    prime_count -= 2;
                }
            }
            const int midpoint = target / 2;
            if (midpoint % 6 == 1 && midpoint > exact_cutoff && least[midpoint] == midpoint) {
                --prime_count;
            }

            int semiprime_upper = semiprime_pairs[index];
            if (prime_count <= 4 * semiprime_upper || target <= 2000 || target % 60006 == 2) {
                fallbacks += prime_count <= 4 * semiprime_upper;
                ++crosschecks;
                const auto exact = direct_counts(target, exact_cutoff, least);
                assert(exact[0] == prime_count && exact[1] <= semiprime_upper);
                semiprime_upper = exact[1];
            }
            if (prime_count <= 4 * semiprime_upper) {
                assert((target == 8 || target == 14) && prime_count == 0 && semiprime_upper == 0);
                ++nonpositive_weighted;
                nonpositive_unweighted += prime_count <= semiprime_upper;
                std::cout << "nonpositive " << target << ' ' << prime_count << ' '
                          << semiprime_upper << '\n';
            }
            if (prime_count > 0 &&
                static_cast<std::int64_t>(semiprime_upper) * worst_prime_pairs >
                    static_cast<std::int64_t>(worst_semiprime_upper) * prime_count) {
                worst_target = target;
                worst_prime_pairs = prime_count;
                worst_semiprime_upper = semiprime_upper;
            }
        }
    }
    std::cout << "limit=" << limit << " targets=" << targets << " direct_checks=" << crosschecks
              << " fallbacks=" << fallbacks << " weighted_nonpositive=" << nonpositive_weighted
              << " unweighted_nonpositive=" << nonpositive_unweighted << '\n';
    std::cout << "worst_bound_target=" << worst_target << " pp=" << worst_prime_pairs
              << " ss_upper=" << worst_semiprime_upper << '\n';
    assert(targets == (limit - 2) / 6 && nonpositive_weighted == 2);
    assert(57 * 4 * static_cast<std::int64_t>(worst_semiprime_upper) <=
           44 * static_cast<std::int64_t>(worst_prime_pairs));
}
