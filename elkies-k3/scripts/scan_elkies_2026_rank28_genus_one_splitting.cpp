// status: ACTIVE_SEARCH
// claim: complete residue-mask scan of the requested primitive projective box
// inputs/outputs: binary tables and candidates managed by the Sage orchestrator

#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t MAGIC = 0x47315331U;  // "G1S1"

template <typename T>
T read_value(std::ifstream& input) {
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!input) {
        throw std::runtime_error("truncated residue table");
    }
    return value;
}

struct Table {
    std::uint32_t prime;
    std::vector<std::uint16_t> masks;
};

int positive_mod(std::int64_t value, int modulus) {
    const int residue = static_cast<int>(value % modulus);
    return residue < 0 ? residue + modulus : residue;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " TABLE HEIGHT OUTPUT\n";
        return 2;
    }
    const std::int64_t height = std::stoll(argv[2]);
    if (height < 1) {
        std::cerr << "HEIGHT must be positive\n";
        return 2;
    }
    std::ifstream input(argv[1], std::ios::binary);
    if (!input || read_value<std::uint32_t>(input) != MAGIC) {
        std::cerr << "invalid residue table\n";
        return 2;
    }
    const auto curve_count = read_value<std::uint32_t>(input);
    const auto table_count = read_value<std::uint32_t>(input);
    if (curve_count == 0 || curve_count > 16 || table_count == 0) {
        std::cerr << "unsupported table dimensions\n";
        return 2;
    }
    std::vector<Table> tables;
    tables.reserve(table_count);
    for (std::uint32_t index = 0; index < table_count; ++index) {
        Table table;
        table.prime = read_value<std::uint32_t>(input);
        const std::size_t size = static_cast<std::size_t>(table.prime) * table.prime;
        table.masks.resize(size);
        input.read(reinterpret_cast<char*>(table.masks.data()),
                   static_cast<std::streamsize>(size * sizeof(std::uint16_t)));
        if (!input) {
            std::cerr << "truncated residue masks\n";
            return 2;
        }
        tables.push_back(std::move(table));
    }

    std::ofstream output(argv[3]);
    if (!output) {
        std::cerr << "cannot open output\n";
        return 2;
    }
    const std::uint16_t all = static_cast<std::uint16_t>((1U << curve_count) - 1U);
    std::uint64_t primitive_count = 0;
    std::uint64_t modular_candidate_count = 0;
    for (std::int64_t denominator = 1; denominator <= height; ++denominator) {
        for (std::int64_t numerator = -height; numerator <= height; ++numerator) {
            if (std::gcd(numerator < 0 ? -numerator : numerator, denominator) != 1) {
                continue;
            }
            ++primitive_count;
            std::uint16_t mask = all;
            for (const auto& table : tables) {
                const int p = static_cast<int>(table.prime);
                const int a = positive_mod(numerator, p);
                const int b = positive_mod(denominator, p);
                mask = static_cast<std::uint16_t>(
                    mask & table.masks[static_cast<std::size_t>(a) * p + b]);
                if (__builtin_popcount(static_cast<unsigned>(mask)) < 2) {
                    break;
                }
            }
            if (__builtin_popcount(static_cast<unsigned>(mask)) >= 2) {
                ++modular_candidate_count;
                output << numerator << ' ' << denominator << ' ' << mask << '\n';
            }
        }
    }
    std::cout << "GENUS1SPLITSCAN|height=" << height
              << "|primitive=" << primitive_count
              << "|modular_candidates=" << modular_candidate_count
              << "|status=PASS\n";
    return 0;
}
