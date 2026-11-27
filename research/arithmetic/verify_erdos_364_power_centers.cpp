// Exact least exponent for the prime tests on an even power center.
#include <charconv>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

using Word = std::uint64_t;

Word read_unsigned(std::string_view value) {
    Word parsed = 0;
    const char* end = value.data() + value.size();
    const auto result = std::from_chars(value.data(), end, parsed);
    if (result.ec != std::errc{} || result.ptr != end) {
        throw std::invalid_argument("expected an unsigned integer");
    }
    return parsed;
}

Word power_mod(Word base, Word exponent, Word modulus) {
    Word result = 1;
    base %= modulus;
    while (exponent != 0) {
        if (exponent % 2 != 0) {
            result = static_cast<__uint128_t>(result) * base % modulus;
        }
        base = static_cast<__uint128_t>(base) * base % modulus;
        exponent /= 2;
    }
    return result;
}

std::vector<std::uint32_t> odd_composite_factors(std::uint32_t bound) {
    std::vector<std::uint32_t> least(static_cast<Word>(bound) / 2 + 1);
    for (Word p = 3; p * p <= bound; p += 2) {
        if (least[p / 2] != 0) {
            continue;
        }
        for (Word n = p * p; n <= bound; n += 2 * p) {
            if (least[n / 2] == 0) {
                least[n / 2] = p;
            }
        }
    }
    return least;
}

std::vector<Word> distinct_factors(Word value, const std::vector<std::uint32_t>& least) {
    std::vector<Word> primes;
    if (value % 2 == 0) {
        primes.push_back(2);
        do {
            value /= 2;
        } while (value % 2 == 0);
    }
    while (value > 1) {
        const Word factor = least[value / 2] == 0 ? value : least[value / 2];
        primes.push_back(factor);
        do {
            value /= factor;
        } while (value % factor == 0);
    }
    return primes;
}

void scan(Word base, std::uint32_t bound) {
    const auto least = odd_composite_factors(bound);
    std::vector<bool> required(static_cast<Word>(bound) / 2 + 1);
    Word prime_count = 0;
    Word forced_count = 0;
    std::cout << "prime,order,half_order\n";
    for (Word p = 3; p <= bound; p += 2) {
        if (least[p / 2] != 0) {
            continue;
        }
        ++prime_count;
        if (base % p == 0) {
            continue;
        }

        // This is gcd(2*M, p-1) for the current squarefree exponent M.
        const auto factors = distinct_factors(p - 1, least);
        Word supported = 2;
        for (Word factor : factors) {
            if (factor != 2 && required[factor / 2]) {
                supported *= factor;
            }
        }
        if (power_mod(base, supported, p) != 1 ||
            power_mod(base, p - 1, p * p) == 1) {
            continue;
        }

        Word order = supported;
        for (Word factor : factors) {
            if (order % factor == 0 && power_mod(base, order / factor, p) == 1) {
                order /= factor;
            }
        }
        required[p / 2] = true;
        ++forced_count;
        std::cout << p << ',' << order << ',' << order / (order % 2 == 0 ? 2 : 1) << '\n';
    }
    std::cerr << "odd primes checked " << prime_count << ", forced " << forced_count << '\n';
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: verify_erdos_364_power_centers EVEN_BASE PRIME_BOUND\n";
        return 2;
    }
    try {
        const Word base = read_unsigned(argv[1]);
        const Word bound = read_unsigned(argv[2]);
        if (base < 2 || base % 2 != 0 || bound < 3 || bound > 1'000'000'000) {
            throw std::invalid_argument("expected an even base >=2 and 3 <= bound <= 10^9");
        }
        scan(base, bound);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
