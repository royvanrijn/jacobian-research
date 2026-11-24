
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <omp.h>

namespace {

using Integer = __uint128_t;
using Word = std::uint64_t;
using Clock = std::chrono::steady_clock;

constexpr Integer max_bound = static_cast<Integer>(1'000'000'000'000ULL) * 1'000'000'000'000ULL;
constexpr unsigned modulus = 16 * 81;

std::string decimal(Integer value) {
    std::string digits;
    do {
        digits.push_back('0' + value % 10);
        value /= 10;
    } while (value != 0);
    std::reverse(digits.begin(), digits.end());
    return digits;
}

Integer parse_bound(const std::string& text) {
    Integer bound = 0;
    for (char digit : text) {
        if (digit < '0' || digit > '9' || bound > (max_bound - (digit - '0')) / 10) {
            throw std::invalid_argument("bound must be an integer from 1 through 10^24");
        }
        bound = 10 * bound + (digit - '0');
    }
    if (bound == 0) {
        throw std::invalid_argument("bound must be an integer from 1 through 10^24");
    }
    return bound;
}

struct IntegerHash {
    std::size_t operator()(Integer value) const {
        return static_cast<Word>(value)
            ^ (static_cast<Word>(value >> 64) * 0x9e3779b97f4a7c15ULL);
    }
};

Integer gcd(Integer a, Integer b) {
    while (b != 0) {
        Integer remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

std::vector<Integer> enumerate_4full(Integer bound) {
    unsigned limit = 1;
    while (static_cast<Integer>(limit + 1) * (limit + 1) * (limit + 1) * (limit + 1) <= bound) {
        ++limit;
    }

    std::vector<bool> composite(limit + 1);
    std::vector<unsigned> primes;
    std::vector<Integer> fourth_powers;
    for (unsigned prime = 2; prime <= limit; ++prime) {
        if (composite[prime]) {
            continue;
        }
        primes.push_back(prime);
        fourth_powers.push_back(static_cast<Integer>(prime) * prime * prime * prime);
        if (static_cast<Word>(prime) * prime <= limit) {
            for (unsigned multiple = prime * prime; multiple <= limit; multiple += prime) {
                composite[multiple] = true;
            }
        }
    }

    std::vector<Integer> values;
    // Increasing primes and exponents >= 4 visit each factorization once.
    const auto visit = [&](auto&& self, Integer value, std::size_t first_prime) -> void {
        values.push_back(value);
        Integer remaining = bound / value;
        for (std::size_t i = first_prime;
             i < primes.size() && fourth_powers[i] <= remaining;
             ++i) {
            Integer power = fourth_powers[i];
            while (true) {
                self(self, value * power, i + 1);
                if (power > remaining / primes[i]) {
                    break;
                }
                power *= primes[i];
            }
        }
    };
    visit(visit, 1, 0);
    std::sort(values.begin(), values.end());
    return values;
}

bool locally_possible(unsigned a, unsigned b) {
    unsigned c = a + b;
    // A primitive pair cannot share 2 or 3. If either prime divides c,
    // fourth-fullness requires its fourth power to divide c as well.
    return !((a % 2 == 0 && b % 2 == 0) || (a % 3 == 0 && b % 3 == 0)
          || (c % 2 == 0 && c % 16 != 0) || (c % 3 == 0 && c % 81 != 0));
}

unsigned bit_index(Integer value) {
    Word key = static_cast<Word>(value >> 4);
    return key ^ (key >> 23) ^ (key >> 45);
}

void write_words(const std::string& path, const std::vector<Integer>& values, std::size_t count) {
    std::ofstream output(path, std::ios::binary);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    for (std::size_t i = 0; i < count; ++i) {
        std::array<unsigned char, 16> bytes;
        for (unsigned byte = 0; byte < bytes.size(); ++byte) {
            bytes[byte] = static_cast<unsigned char>(values[i] >> (8 * byte));
        }
        output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    }
    output.close();
}

struct Task {
    unsigned first;
    unsigned second;
    std::size_t start;
    std::size_t stop;
};

void search(Integer bound, const char* dump_prefix) {
    auto begin = Clock::now();
    const auto values = enumerate_4full(2 * bound);
    const auto last_summand = std::upper_bound(values.begin(), values.end(), bound);
    const Word summand_count = last_summand - values.begin();

    if (dump_prefix != nullptr) {
        write_words(std::string(dump_prefix) + "-summands.bin", values, summand_count);
        write_words(std::string(dump_prefix) + "-targets.bin", values, values.size());
    }
    std::cerr << "B " << decimal(bound) << ", summands " << summand_count
              << ", targets " << values.size() << ", threads " << omp_get_max_threads()
              << std::endl;

    std::unordered_set<Integer, IntegerHash> targets;
    targets.max_load_factor(0.65);
    targets.reserve(values.size());

    std::array<Word, modulus> target_counts{};
    for (Integer value : values) {
        ++target_counts[value % modulus];
    }

    // Each residue has its own bit filter. An unset bit proves absence.
    // A set bit may be a collision and always requires exact lookup.
    std::array<std::vector<Word>, modulus> filters;
    for (unsigned residue = 0; residue < modulus; ++residue) {
        Word bits = 64;
        while (bits < 128 * target_counts[residue]) {
            bits *= 2;
        }
        filters[residue].resize(bits / 64);
    }
    for (Integer value : values) {
        targets.insert(value);
        auto& filter = filters[value % modulus];
        unsigned index = bit_index(value) & (64 * filter.size() - 1);
        filter[index >> 6] |= Word{1} << (index & 63);
    }

    std::array<std::vector<Integer>, modulus> buckets;
    for (auto it = values.begin(); it != last_summand; ++it) {
        buckets[*it % modulus].push_back(*it);
    }

    std::vector<Task> tasks;
    Word surviving_pairs = 0;
    unsigned active_buckets = 0;
    for (unsigned i = 0; i < modulus; ++i) {
        if (buckets[i].empty()) {
            continue;
        }
        ++active_buckets;
        for (unsigned j = i; j < modulus; ++j) {
            if (buckets[j].empty() || target_counts[(i + j) % modulus] == 0 || !locally_possible(i, j)) {
                continue;
            }
            const auto m = buckets[i].size();
            const auto k = buckets[j].size();
            surviving_pairs += i == j ? m * (m + 1) / 2 : m * k;

            // Visit each unordered pair once. Keep the larger bucket in the
            // inner loop and split the outer bucket to balance thread work.
            unsigned left = i;
            unsigned right = j;
            if (m > k) {
                std::swap(left, right);
            }
            const auto stride = std::max<std::size_t>(1, 500'000 / buckets[right].size());
            for (std::size_t start = 0; start < buckets[left].size(); start += stride) {
                tasks.push_back({left, right, start, std::min(start + stride, buckets[left].size())});
            }
        }
    }
    std::cerr << "unordered_pairs " << summand_count * (summand_count + 1) / 2
              << ", surviving_pairs_after_target_residues " << surviving_pairs
              << ", active_buckets " << active_buckets << ", work_units " << tasks.size() << std::endl;

    const auto& exact_targets = targets;
    Word powerful_sum_candidates = 0;
    Word coprime_hits = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:powerful_sum_candidates,coprime_hits)
    for (std::size_t task_index = 0; task_index < tasks.size(); ++task_index) {
        const auto& task = tasks[task_index];
        const auto& left = buckets[task.first];
        const auto& right = buckets[task.second];
        const bool diagonal = task.first == task.second;
        const auto& filter = filters[(task.first + task.second) % modulus];
        const unsigned mask = 64 * filter.size() - 1;
        for (std::size_t i = task.start; i < task.stop; ++i) {
            const Integer a = left[i];
            for (std::size_t j = diagonal ? i : 0; j < right.size(); ++j) {
                const Integer b = right[j];
                const Integer c = a + b;
                const unsigned index = bit_index(c) & mask;
                if (!(filter[index >> 6] & (Word{1} << (index & 63)))) {
                    continue;
                }
                if (exact_targets.find(c) == exact_targets.end()) {
                    continue;
                }
                ++powerful_sum_candidates;
                if (gcd(a, b) != 1) {
                    continue;
                }
                ++coprime_hits;
                #pragma omp critical
                std::cout << "hit " << decimal(std::min(a, b)) << " " << decimal(std::max(a, b))
                          << " " << decimal(c) << std::endl;
            }
        }
    }

    std::cout << "powerful_sum_candidates " << powerful_sum_candidates
              << ", coprime_hits " << coprime_hits
              << ", seconds " << std::chrono::duration<double>(Clock::now() - begin).count()
              << std::endl;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " max-summand [binary-dump-prefix]\n";
        return 2;
    }

    try {
        search(parse_bound(argv[1]), argc == 3 ? argv[2] : nullptr);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}

