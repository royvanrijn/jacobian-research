// status: ACTIVE_SEARCH
// claim: complete mixed-trace residue scan of a primitive projective box
// inputs/outputs: binary tables and candidates managed by the Sage orchestrator

#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t MAGIC = 0x4d475331U;  // "MGS1"

template <typename T>
T read_value(std::ifstream& input) {
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!input) {
        throw std::runtime_error("truncated residue table");
    }
    return value;
}

struct Mask {
    std::uint64_t low;
    std::uint64_t high;
};

struct Table {
    std::uint32_t prime;
    std::vector<Mask> masks;
};

int positive_mod(std::int64_t value, int modulus) {
    const int residue = static_cast<int>(value % modulus);
    return residue < 0 ? residue + modulus : residue;
}

bool mixed_trace(const Mask& mask, std::uint32_t targets_per_trace) {
    if (mask.low == 0 && mask.high == 0) {
        return false;
    }
    const unsigned first = mask.low != 0
        ? static_cast<unsigned>(__builtin_ctzll(mask.low))
        : 64U + static_cast<unsigned>(__builtin_ctzll(mask.high));
    const unsigned trace = first / targets_per_trace;
    std::uint64_t trace_low = 0;
    std::uint64_t trace_high = 0;
    const unsigned begin = trace * targets_per_trace;
    const unsigned end = begin + targets_per_trace;
    for (unsigned index = begin; index < end; ++index) {
        if (index < 64) {
            trace_low |= std::uint64_t{1} << index;
        } else {
            trace_high |= std::uint64_t{1} << (index - 64);
        }
    }
    return (mask.low & ~trace_low) != 0 || (mask.high & ~trace_high) != 0;
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
    const auto trace_count = read_value<std::uint32_t>(input);
    const auto targets_per_trace = read_value<std::uint32_t>(input);
    const auto table_count = read_value<std::uint32_t>(input);
    if (curve_count == 0 || curve_count > 128 || trace_count < 2 ||
        curve_count != trace_count * targets_per_trace || table_count == 0) {
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
                   static_cast<std::streamsize>(size * sizeof(Mask)));
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
    Mask all{
        curve_count >= 64 ? ~std::uint64_t{0}
                          : ((std::uint64_t{1} << curve_count) - 1),
        curve_count <= 64 ? 0
                          : ((std::uint64_t{1} << (curve_count - 64)) - 1),
    };
    std::uint64_t primitive_count = 0;
    std::uint64_t modular_candidate_count = 0;
    for (std::int64_t denominator = 1; denominator <= height; ++denominator) {
        for (std::int64_t numerator = -height; numerator <= height; ++numerator) {
            if (std::gcd(numerator < 0 ? -numerator : numerator, denominator) != 1) {
                continue;
            }
            ++primitive_count;
            Mask mask = all;
            for (const auto& table : tables) {
                const int p = static_cast<int>(table.prime);
                const int a = positive_mod(numerator, p);
                const int b = positive_mod(denominator, p);
                const auto& local = table.masks[static_cast<std::size_t>(a) * p + b];
                mask.low &= local.low;
                mask.high &= local.high;
                if (!mixed_trace(mask, targets_per_trace)) {
                    break;
                }
            }
            if (mixed_trace(mask, targets_per_trace)) {
                ++modular_candidate_count;
                output << numerator << ' ' << denominator << ' '
                       << mask.low << ' ' << mask.high << '\n';
            }
        }
    }
    std::cout << "MIXEDTRACEG1SCAN|height=" << height
              << "|primitive=" << primitive_count
              << "|modular_candidates=" << modular_candidate_count
              << "|status=PASS\n";
    return 0;
}
