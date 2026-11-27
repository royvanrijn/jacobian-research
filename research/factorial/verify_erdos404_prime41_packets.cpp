#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using Word = std::uint64_t;
using Wide = unsigned __int128;
using Bits = std::vector<Word>;

constexpr Word prime = 41;
constexpr Word sixth_power = 4'750'104'241;
constexpr Word seventh_power = prime * sixth_power;
constexpr Word eighth_power = prime * seventh_power;
constexpr std::array<std::array<unsigned, 2>, 8> pairs{{
    {0, 1}, {5, 7}, {6, 19}, {8, 22}, {9, 26}, {11, 15}, {18, 32}, {23, 30},
}};

Word multiply_mod(Word left, Word right, Word modulus) {
    return static_cast<Wide>(left) * right % modulus;
}

Word subtract_mod(Word a, Word b, Word modulus) {
    return a >= b ? a - b : a + modulus - b;
}

Word last_mask(Word modulus) {
    const unsigned tail = modulus % 64;
    return tail == 0 ? ~Word{0} : (Word{1} << tail) - 1;
}

// Read one word from a cyclic bitset whose logical length may not be word-aligned.
Word cyclic_word(const Bits& bits, Word start, Word modulus) {
    const Word index = start / 64;
    const unsigned offset = start % 64;
    Word result = bits[index] >> offset;
    if (offset != 0 && index + 1 < bits.size()) {
        result |= bits[index + 1] << (64 - offset);
    }
    const Word tail = modulus - start;
    if (tail < 64) {
        result = (result & ((Word{1} << tail) - 1)) | (bits.front() << tail);
    }
    return result;
}

Bits last_block_support(Word modulus) {
    Bits support((modulus + 63) / 64);
    Bits next(support.size());
    Word weight = 1;
    support.front() = 1;
    for (Word j = 0; j < prime; ++j) {
        if (j != 0) {
            weight = multiply_mod(weight, modulus - prime + j, modulus);
        }
        #pragma omp parallel for schedule(static)
        for (Word index = 0; index < support.size(); ++index) {
            const Word start = subtract_mod(64 * index, weight, modulus);
            next[index] = support[index] | cyclic_word(support, start, modulus);
        }
        next.back() &= last_mask(modulus);
        support.swap(next);
    }
    return support;
}

std::vector<Word> missing_residues(const Bits& support, Word modulus) {
    std::vector<Word> missing;
    for (Word index = 0; index < support.size(); ++index) {
        Word absent = ~support[index];
        if (index + 1 == support.size()) {
            absent &= last_mask(modulus);
        }
        while (absent != 0) {
            missing.push_back(64 * index + std::countr_zero(absent));
            absent &= absent - 1;
        }
    }
    return missing;
}

Word inverse_unit(Word unit) {
    Word exponent = seventh_power - sixth_power - 1;
    Word inverse = 1;
    while (exponent != 0) {
        if (exponent % 2 != 0) {
            inverse = multiply_mod(inverse, unit, seventh_power);
        }
        unit = multiply_mod(unit, unit, seventh_power);
        exponent /= 2;
    }
    return inverse;
}

std::vector<Word> prefix_shifts(Word residue) {
    // This positive boundary has the specified high digit modulo 41^8.
    const Word boundary = (prime + residue) * seventh_power;
    const Word start = boundary - 2 * prime;
    const Word split = boundary - prime;
    std::array<Word, prime> weights{};
    weights[0] = 1;
    for (Word j = 1; j < weights.size(); ++j) {
        weights[j] = multiply_mod(weights[j - 1], start + j, eighth_power);
    }

    std::array<bool, prime> used{};
    for (auto [first, second] : pairs) {
        if (used[first] || used[second] || weights[first] % prime != weights[second] % prime) {
            return {};
        }
        used[first] = used[second] = true;
    }

    // The unused terms provide an anchor in every low residue, after the
    // common low residue of the pairs has been added.
    std::array<bool, prime> low_support{};
    low_support[0] = true;
    for (Word j = 0; j < weights.size(); ++j) {
        if (!used[j]) {
            auto previous = low_support;
            for (Word low = 0; low < prime; ++low) {
                if (previous[low]) {
                    low_support[(low + weights[j]) % prime] = true;
                }
            }
        }
    }
    if (!std::all_of(low_support.begin(), low_support.end(), [](bool bit) { return bit; })) {
        return {};
    }

    // split!/start! has valuation one. Its unit rescales prefix shifts.
    const Word unit = multiply_mod(weights.back() % seventh_power,
                                   split / prime, seventh_power);
    if (unit % prime == 0) {
        return {};
    }
    const Word inverse = inverse_unit(unit);
    std::vector<Word> shifts{0};
    for (auto [first, second] : pairs) {
        const Word difference = subtract_mod(weights[second], weights[first], eighth_power);
        const Word shift = multiply_mod(difference / prime, inverse, seventh_power);
        const std::size_t length = shifts.size();
        for (std::size_t j = 0; j < length; ++j) {
            shifts.push_back((shifts[j] + shift) % seventh_power);
        }
    }

    return shifts;
}

struct MissingColumn {
    Word low;
    Word high_bits;
};

std::vector<MissingColumn> missing_columns(const std::vector<Word>& missing) {
    std::vector<Word> transposed;
    transposed.reserve(missing.size());
    for (Word value : missing) {
        transposed.push_back(prime * (value % sixth_power) + value / sixth_power);
    }
    std::sort(transposed.begin(), transposed.end());
    std::vector<MissingColumn> columns;
    for (Word value : transposed) {
        const Word low = value / prime;
        if (columns.empty() || columns.back().low != low) {
            columns.push_back({low, 0});
        }
        columns.back().high_bits |= Word{1} << (value % prime);
    }
    return columns;
}

// All first shifts have the same reduction modulo 41^6. Find their
// intersections with the missing set in one pass, using the high digits.
std::vector<std::vector<Word>> first_intersections(
    const std::vector<MissingColumn>& columns,
    const std::array<std::vector<Word>, prime>& shifts) {
    const Word shift_low = shifts[1][1] % sixth_power;
    std::array<std::vector<Word>, prime> cases_by_high;
    for (Word residue = 1; residue < prime; ++residue) {
        const Word shift = shifts[residue][1];
        if (shift % sixth_power != shift_low) {
            return {};
        }
        cases_by_high[shift / sixth_power].push_back(residue);
    }

    std::vector<std::vector<Word>> hits(prime);
    std::size_t match = 0;
    Word previous_low = 0;
    for (const auto& column : columns) {
        const Word low = subtract_mod(column.low, shift_low, sixth_power);
        if (low < previous_low) {
            match = 0;
        }
        previous_low = low;
        while (match < columns.size() && columns[match].low < low) {
            ++match;
        }
        if (match == columns.size() || columns[match].low != low) {
            continue;
        }
        for (Word first = column.high_bits; first != 0; first &= first - 1) {
            const Word a = column.low + sixth_power * std::countr_zero(first);
            for (Word second = columns[match].high_bits; second != 0; second &= second - 1) {
                const Word b = low + sixth_power * std::countr_zero(second);
                const Word high = subtract_mod(a, b, seventh_power) / sixth_power;
                for (Word residue : cases_by_high[high]) {
                    hits[residue].push_back(a);
                }
            }
        }
    }
    return hits;
}

}  // namespace

int main() {
    if (!missing_residues(last_block_support(sixth_power), sixth_power).empty()) {
        std::cerr << "A single block fails to cover modulo 41^6\n";
        return 1;
    }
    std::cout << "One block covers all 41^6 = " << sixth_power << " residues\n" << std::flush;

    const auto suffix = last_block_support(seventh_power);
    const auto missing = missing_residues(suffix, seventh_power);
    if (missing.empty()) {
        std::cerr << "Expected a nonempty complement modulo 41^7\n";
        return 1;
    }
    std::cout << "One block misses " << missing.size() << " residues modulo 41^7\n" << std::flush;
    std::array<std::vector<Word>, prime> shifts;
    for (Word residue = 1; residue < prime; ++residue) {
        shifts[residue] = prefix_shifts(residue);
        if (shifts[residue].size() < 2) {
            std::cerr << "Prefix certificate failed for high digit " << residue << '\n';
            return 1;
        }
    }
    auto hits = first_intersections(missing_columns(missing), shifts);
    if (hits.size() != prime) {
        std::cerr << "First shifts have inconsistent low residues\n";
        return 1;
    }
    for (Word residue = 1; residue < prime; ++residue) {
        auto& uncovered = hits[residue];
        for (std::size_t j = 2; j < shifts[residue].size() && !uncovered.empty(); ++j) {
            const Word shift = shifts[residue][j];
            std::erase_if(uncovered, [&](Word value) {
                const Word source = subtract_mod(value, shift, seventh_power);
                return (suffix[source / 64] >> (source % 64)) & 1;
            });
        }
        if (!uncovered.empty()) {
            std::cerr << "Two-block certificate failed for high digit " << residue << '\n';
            return 1;
        }
    }
    std::cout << "Two blocks cover the size-seven jump in all 40 boundary classes\n";
}

