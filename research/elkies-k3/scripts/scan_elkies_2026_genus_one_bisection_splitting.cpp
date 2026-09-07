// status: ACTIVE_SEARCH
// Dynamic simultaneous Legendre-mask scanner for the R17 genus-one covers.

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr std::uint32_t MAGIC = 0x47315332U;  // "G1S2"

template <typename T> T read_value(std::ifstream& input) {
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!input) throw std::runtime_error("truncated residue table");
    return value;
}

struct PrimeTable {
    std::uint32_t prime;
    std::uint32_t block;
    std::vector<std::uint64_t> masks;
};

struct Score {
    std::uint32_t weakest_rank{};
    std::uint32_t weakest_count{};
    std::uint32_t rank_sum{};
    std::uint32_t count_sum{};
    std::uint32_t full_count{};
    auto tuple() const {
        return std::tie(weakest_rank, weakest_count, rank_sum, count_sum, full_count);
    }
    bool operator<(const Score& other) const { return tuple() < other.tuple(); }
};

struct RankedRow {
    Score score;
    std::int64_t numerator{};
    std::int64_t denominator{};
};

struct WorseFirst {
    bool operator()(const RankedRow& left, const RankedRow& right) const {
        return right.score < left.score;
    }
};

std::int64_t positive_mod(std::int64_t value, std::int64_t modulus) {
    const auto residue = value % modulus;
    return residue < 0 ? residue + modulus : residue;
}

std::uint32_t mask_count(const std::vector<std::uint64_t>& words) {
    std::uint32_t count = 0;
    for (auto word : words) count += static_cast<std::uint32_t>(__builtin_popcountll(word));
    return count;
}

std::uint32_t lattice_rank(
    const std::vector<std::uint64_t>& words,
    const std::vector<std::uint32_t>& lattice_masks,
    std::uint32_t curve_count
) {
    std::uint32_t pivots[17]{};
    std::uint32_t rank = 0;
    for (std::uint32_t index = 0; index < curve_count; ++index) {
        if (((words[index / 64] >> (index % 64)) & 1U) == 0) continue;
        std::uint32_t value = lattice_masks[index];
        while (value) {
            const unsigned bit = 31U - static_cast<unsigned>(__builtin_clz(value));
            if (pivots[bit]) {
                value ^= pivots[bit];
            } else {
                pivots[bit] = value;
                ++rank;
                break;
            }
        }
    }
    return rank;
}

std::uint64_t xorshift64star(std::uint64_t& state) {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * UINT64_C(2685821657736338717);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "usage: " << argv[0]
                  << " TABLE RANDOM_COUNT COORD_BOUND BOX_HEIGHT TOP_K SEED OUTPUT\n";
        return 2;
    }
    const std::uint64_t random_count = std::stoull(argv[2]);
    const std::int64_t coordinate_bound = std::stoll(argv[3]);
    const std::int64_t box_height = std::stoll(argv[4]);
    const std::size_t top_k = static_cast<std::size_t>(std::stoull(argv[5]));
    std::uint64_t random_state = std::stoull(argv[6]);
    if (coordinate_bound < 1 || box_height < 0 || top_k < 1 || random_state == 0) {
        std::cerr << "invalid scan bounds\n";
        return 2;
    }

    std::ifstream input(argv[1], std::ios::binary);
    if (!input || read_value<std::uint32_t>(input) != MAGIC) {
        std::cerr << "invalid residue table\n";
        return 2;
    }
    const auto curve_count = read_value<std::uint32_t>(input);
    const auto word_count = read_value<std::uint32_t>(input);
    const auto prime_count = read_value<std::uint32_t>(input);
    const auto block_count = read_value<std::uint32_t>(input);
    if (!curve_count || word_count != (curve_count + 63) / 64 || !prime_count || block_count < 2) {
        std::cerr << "unsupported residue-table dimensions\n";
        return 2;
    }
    std::vector<std::uint32_t> lattice_masks(curve_count);
    input.read(reinterpret_cast<char*>(lattice_masks.data()),
               static_cast<std::streamsize>(curve_count * sizeof(std::uint32_t)));
    if (!input) throw std::runtime_error("truncated lattice-mask header");

    std::vector<PrimeTable> tables;
    for (std::uint32_t index = 0; index < prime_count; ++index) {
        PrimeTable table;
        table.prime = read_value<std::uint32_t>(input);
        table.block = read_value<std::uint32_t>(input);
        if (table.block >= block_count) throw std::runtime_error("invalid prime block");
        const std::size_t entries = static_cast<std::size_t>(table.prime) * table.prime * word_count;
        table.masks.resize(entries);
        input.read(reinterpret_cast<char*>(table.masks.data()),
                   static_cast<std::streamsize>(entries * sizeof(std::uint64_t)));
        if (!input) throw std::runtime_error("truncated prime table");
        tables.push_back(std::move(table));
    }

    std::ofstream output(argv[7]);
    if (!output) {
        std::cerr << "cannot open output\n";
        return 2;
    }
    std::priority_queue<RankedRow, std::vector<RankedRow>, WorseFirst> top;
    std::uint64_t evaluated = 0;
    std::uint64_t extremes = 0;
    const std::uint64_t final_word_mask = curve_count % 64
        ? ((UINT64_C(1) << (curve_count % 64)) - 1) : ~UINT64_C(0);

    auto evaluate = [&](std::int64_t numerator, std::int64_t denominator) {
        if (denominator < 1 || std::gcd(numerator < 0 ? -numerator : numerator, denominator) != 1) return;
        ++evaluated;
        std::vector<std::vector<std::uint64_t>> blocks(
            block_count, std::vector<std::uint64_t>(word_count, ~UINT64_C(0)));
        for (auto& block : blocks) block.back() &= final_word_mask;
        for (const auto& table : tables) {
            const auto p = static_cast<std::int64_t>(table.prime);
            const auto a = positive_mod(numerator, p);
            const auto b = positive_mod(denominator, p);
            const std::size_t offset = (static_cast<std::size_t>(a) * p + b) * word_count;
            for (std::uint32_t word = 0; word < word_count; ++word) {
                blocks[table.block][word] &= table.masks[offset + word];
            }
        }
        Score score{17U, curve_count, 0U, 0U, 0U};
        std::vector<std::uint64_t> full(word_count, ~UINT64_C(0));
        full.back() &= final_word_mask;
        for (const auto& block : blocks) {
            const auto count = mask_count(block);
            const auto rank = lattice_rank(block, lattice_masks, curve_count);
            score.weakest_count = std::min(score.weakest_count, count);
            score.weakest_rank = std::min(score.weakest_rank, rank);
            score.count_sum += count;
            score.rank_sum += rank;
            for (std::uint32_t word = 0; word < word_count; ++word) full[word] &= block[word];
        }
        score.full_count = mask_count(full);
        if (score.full_count >= 2) {
            ++extremes;
            output << "E " << numerator << ' ' << denominator << ' '
                   << score.weakest_rank << ' ' << score.weakest_count << ' '
                   << score.rank_sum << ' ' << score.count_sum << ' ' << score.full_count;
            output << std::hex;
            for (auto word : full) output << ' ' << word;
            output << std::dec << '\n';
        }
        RankedRow row{score, numerator, denominator};
        if (top.size() < top_k) top.push(row);
        else if (top.top().score < score) { top.pop(); top.push(row); }
    };

    if (box_height > 0) {
        for (std::int64_t denominator = 1; denominator <= box_height; ++denominator) {
            for (std::int64_t numerator = -box_height; numerator <= box_height; ++numerator) {
                evaluate(numerator, denominator);
            }
        }
    }
    for (std::uint64_t index = 0; index < random_count; ++index) {
        const auto numerator = static_cast<std::int64_t>(
            xorshift64star(random_state) % (2 * static_cast<std::uint64_t>(coordinate_bound) + 1))
            - coordinate_bound;
        const auto denominator = static_cast<std::int64_t>(
            xorshift64star(random_state) % static_cast<std::uint64_t>(coordinate_bound)) + 1;
        evaluate(numerator, denominator);
    }

    std::vector<RankedRow> ranked;
    while (!top.empty()) { ranked.push_back(top.top()); top.pop(); }
    std::sort(ranked.begin(), ranked.end(), [](const auto& left, const auto& right) {
        return right.score < left.score;
    });
    for (const auto& row : ranked) {
        output << "R " << row.numerator << ' ' << row.denominator << ' '
               << row.score.weakest_rank << ' ' << row.score.weakest_count << ' '
               << row.score.rank_sum << ' ' << row.score.count_sum << ' '
               << row.score.full_count << '\n';
    }
    std::cout << "G1SPLITV2|curves=" << curve_count << "|words=" << word_count
              << "|primes=" << prime_count << "|blocks=" << block_count
              << "|evaluated=" << evaluated << "|extremes=" << extremes
              << "|status=PASS\n";
    return 0;
}
