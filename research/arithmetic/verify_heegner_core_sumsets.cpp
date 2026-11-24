#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Word = std::uint64_t;
using Bits = std::vector<Word>;

bool bit(const Bits& bits, Word index) {
    return (bits[index / 64] >> (index % 64)) & 1;
}

void clear(Bits& bits, Word index) {
    bits[index / 64] &= ~(Word{1} << (index % 64));
}

Bits split_core(Word bound, unsigned p) {
    // Bit i represents the odd integer 2i+1.
    const Word odd_count = (bound + 1) / 2;
    Bits primes((odd_count + 63) / 64, ~Word{0});
    clear(primes, 0);
    for (Word q = 3; q <= bound / q; q += 2)
        if (bit(primes, q / 2))
            for (Word j = q*q / 2; j < odd_count; j += q)
                clear(primes, j);

    std::vector<bool> square(p);
    for (unsigned j = 1; j < p; ++j)
        square[j * j % p] = true;
    Bits core(primes.size(), ~Word{0});
    for (Word q = 3; q <= bound; q += 2)
        if (bit(primes, q / 2) && !square[q % p])
            for (Word j = q / 2; j < odd_count; j += q)
                clear(core, j);
    if (odd_count % 64)
        core.back() &= (Word{1} << (odd_count % 64)) - 1;
    return core;
}

void check_coordinates(const Bits& core, Word bound, unsigned p) {
    const int limit = static_cast<int>(std::min(bound, Word{10000}));
    std::vector<bool> represented(limit + 1);
    for (int y = 0; static_cast<int>(p) * y * y <= 4 * limit; ++y) {
        for (int x = -y / 2; (2*x+y)*(2*x+y) + static_cast<int>(p)*y*y <= 4*limit; ++x) {
            const int n = x*x + x*y + static_cast<int>((p+1)/4)*y*y;
            if (n > 0 && n % p && std::gcd(x, y) == 1)
                represented[n] = true;
        }
    }
    for (int n = 1; n <= limit; ++n)
        if (represented[n] != (n % 2 && bit(core, n / 2)))
            throw std::runtime_error("split sieve and coordinate enumeration disagree");
}

Bits sum_witnesses(const Bits& core) {
    Bits covered(core.size());
    // If i and j encode a and b, bit i+j encodes the sum 2(i+j+1).
    for (Word a = 1; a <= 5000 && a / 2 < core.size() * 64; a += 2) {
        if (!bit(core, a / 2))
            continue;
        const Word offset = a / 2;
        const Word words = offset / 64;
        const unsigned shift = offset % 64;
        if (shift == 0) {
            for (Word j = 0; j + words < core.size(); ++j)
                covered[j + words] |= core[j];
        } else {
            covered[words] |= core[0] << shift;
            for (Word j = 1; j + words < core.size(); ++j)
                covered[j + words] |= (core[j] << shift) | (core[j - 1] >> (64 - shift));
        }
    }
    return covered;
}

void array(const std::vector<Word>& values) {
    std::cout << "[";
    for (size_t i = 0; i < values.size(); ++i)
        std::cout << (i ? "," : "") << values[i];
    std::cout << "]";
}

void scan(Word bound, unsigned p) {
    const Bits core = split_core(bound, p);
    check_coordinates(core, bound, p);
    const Bits covered = sum_witnesses(core);
    std::vector<Word> exceptions;
    std::vector<Word> unit_primitive_exceptions;
    Word residual_targets = 0;
    for (Word j = 0; j < covered.size(); ++j) {
        Word missing = ~covered[j];
        while (missing) {
            const Word n = 2 * (64*j + std::countr_zero(missing) + 1);
            missing &= missing - 1;
            if (n > bound || n % p == 0)
                continue;
            ++residual_targets;
            bool found = false;
            for (Word a = 1; a <= n / 2; a += 2) {
                if (bit(core, a / 2) && bit(core, (n-a) / 2)) {
                    found = true;
                    break;
                }
            }
            if (found)
                continue;
            exceptions.push_back(n);
            for (Word a = 1; p * a < n; a += 2) {
                if (bit(core, a / 2) && bit(core, (n-p*a) / 2)) {
                    found = true;
                    break;
                }
            }
            if (!found)
                unit_primitive_exceptions.push_back(n);
        }
    }

    std::cout << "{\"p\":" << p << ",\"bound\":" << bound
              << ",\"core_value_count\":" << std::accumulate(core.begin(), core.end(), Word{0},
                    [](Word count, Word word) { return count + std::popcount(word); })
              << ",\"exhaustive_residual_targets\":" << residual_targets
              << ",\"core_exceptions\":";
    array(exceptions);
    std::cout << ",\"unit_primitive_exceptions\":";
    array(unit_primitive_exceptions);
    std::cout << "}\n" << std::flush;
}

int main(int argc, char** argv) {
    if (argc < 3)
        throw std::invalid_argument("provide an upper bound and one or more ramified primes");
    const Word bound = std::stoull(argv[1]);
    if (bound < 2 || bound > (Word{1} << 36))
        throw std::invalid_argument("bound must be between 2 and 2^36");
    for (int i = 2; i < argc; ++i) {
        const auto p = std::stoul(argv[i]);
        if (p != 11 && p != 19 && p != 43 && p != 67 && p != 163)
            throw std::invalid_argument("unsupported discriminant");
        scan(bound, static_cast<unsigned>(p));
    }
}
