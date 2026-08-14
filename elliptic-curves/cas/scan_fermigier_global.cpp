// Exhaustive, leakage-free global scanner for the even Fermigier--Mestre family.
//
// The Python driver supplies two disjoint prime bands.  Each table contains a
// good-reduction rank score modulo p and an independently useful repeated-
// discriminant-prime score modulo p^2.  Only the discovery band participates
// in retention.  Held-band values are computed for retained rows but never
// enter a discovery heap.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::int64_t kCompositePowerDivisor = 16;

struct Table {
  int prime = 0;
  int power_modulus = 0;
  std::vector<std::int64_t> rank_weights;
  std::vector<std::int64_t> power_weights;
};

struct Candidate {
  int numerator = 0;
  int denominator = 0;
  std::int64_t discovery_rank = 0;
  std::int64_t discovery_power = 0;
  std::int64_t held_rank = 0;
  std::int64_t held_power = 0;
};

enum class Metric { Rank, Power, Composite };

std::int64_t metric_value(const Candidate& candidate, Metric metric) {
  if (metric == Metric::Rank) return candidate.discovery_rank;
  if (metric == Metric::Power) return candidate.discovery_power;
  return candidate.discovery_rank +
         candidate.discovery_power / kCompositePowerDivisor;
}

auto quality(const Candidate& candidate, Metric metric) {
  return std::tuple(metric_value(candidate, metric),
                    -std::max(candidate.numerator, candidate.denominator),
                    -candidate.denominator, -candidate.numerator);
}

struct Better {
  Metric metric = Metric::Composite;
  bool operator()(const Candidate& left, const Candidate& right) const {
    return quality(left, metric) > quality(right, metric);
  }
};

using Heap = std::priority_queue<Candidate, std::vector<Candidate>, Better>;

long parse_long(const char* text, const char* label) {
  char* end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (end == text || *end != '\0') {
    throw std::invalid_argument(std::string("invalid ") + label);
  }
  return value;
}

std::vector<Table> read_band() {
  int count = 0;
  if (!(std::cin >> count) || count < 1) {
    throw std::runtime_error("invalid or empty prime-band header");
  }
  std::vector<Table> result;
  result.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index) {
    Table table;
    if (!(std::cin >> table.prime) || table.prime < 5) {
      throw std::runtime_error("invalid table prime");
    }
    table.power_modulus = table.prime * table.prime;
    table.rank_weights.resize(static_cast<std::size_t>(table.prime + 1));
    table.power_weights.resize(
        static_cast<std::size_t>(table.power_modulus + 1));
    for (auto& value : table.rank_weights) {
      if (!(std::cin >> value)) throw std::runtime_error("truncated rank table");
    }
    for (auto& value : table.power_weights) {
      if (!(std::cin >> value)) throw std::runtime_error("truncated power table");
    }
    result.push_back(std::move(table));
  }
  return result;
}

int inverse_mod(int value, int modulus) {
  int old_r = value;
  int r = modulus;
  int old_s = 1;
  int s = 0;
  while (r != 0) {
    const int quotient = old_r / r;
    const int next_r = old_r - quotient * r;
    old_r = r;
    r = next_r;
    const int next_s = old_s - quotient * s;
    old_s = s;
    s = next_s;
  }
  if (old_r != 1) throw std::runtime_error("noninvertible denominator");
  old_s %= modulus;
  if (old_s < 0) old_s += modulus;
  return old_s;
}

void retain(Heap& heap, const Candidate& candidate, std::size_t limit,
            Metric metric) {
  if (limit == 0) return;
  if (heap.size() < limit) {
    heap.push(candidate);
  } else if (quality(candidate, metric) > quality(heap.top(), metric)) {
    heap.pop();
    heap.push(candidate);
  }
}

std::vector<Candidate> drain(Heap& heap, Metric metric) {
  std::vector<Candidate> answer;
  answer.reserve(heap.size());
  while (!heap.empty()) {
    answer.push_back(heap.top());
    heap.pop();
  }
  std::sort(answer.begin(), answer.end(), [metric](const Candidate& left,
                                                   const Candidate& right) {
    if (quality(left, metric) != quality(right, metric))
      return quality(left, metric) > quality(right, metric);
    return std::tie(left.denominator, left.numerator) <
           std::tie(right.denominator, right.numerator);
  });
  return answer;
}

struct BandState {
  std::vector<int> rank_indices;
  std::vector<int> rank_steps;
  std::vector<int> power_indices;
  std::vector<int> power_steps;
};

BandState initial_state(const std::vector<Table>& tables, int denominator) {
  BandState state;
  state.rank_indices.reserve(tables.size());
  state.rank_steps.reserve(tables.size());
  state.power_indices.reserve(tables.size());
  state.power_steps.reserve(tables.size());
  for (const auto& table : tables) {
    if (denominator % table.prime == 0) {
      state.rank_indices.push_back(table.prime);  // projective infinity
      state.rank_steps.push_back(0);
      state.power_indices.push_back(table.power_modulus);
      state.power_steps.push_back(0);
    } else {
      const int rank_step = inverse_mod(denominator % table.prime, table.prime);
      const int power_step =
          inverse_mod(denominator % table.power_modulus, table.power_modulus);
      state.rank_indices.push_back(0);  // numerator starts at zero
      state.rank_steps.push_back(rank_step);
      state.power_indices.push_back(0);
      state.power_steps.push_back(power_step);
    }
  }
  return state;
}

std::pair<std::int64_t, std::int64_t> score_band(
    const std::vector<Table>& tables, const BandState& state) {
  std::int64_t rank = 0;
  std::int64_t power = 0;
  for (std::size_t index = 0; index < tables.size(); ++index) {
    rank += tables[index].rank_weights[
        static_cast<std::size_t>(state.rank_indices[index])];
    power += tables[index].power_weights[
        static_cast<std::size_t>(state.power_indices[index])];
  }
  return {rank, power};
}

void advance(const std::vector<Table>& tables, BandState& state) {
  for (std::size_t index = 0; index < tables.size(); ++index) {
    if (state.rank_steps[index] != 0) {
      state.rank_indices[index] += state.rank_steps[index];
      if (state.rank_indices[index] >= tables[index].prime)
        state.rank_indices[index] -= tables[index].prime;
      state.power_indices[index] += state.power_steps[index];
      if (state.power_indices[index] >= tables[index].power_modulus)
        state.power_indices[index] -= tables[index].power_modulus;
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 7) {
      std::cerr << "usage: scan_fermigier_global A_MAX B_MAX RANK_KEEP "
                   "POWER_KEEP COMPOSITE_KEEP PER_DENOMINATOR_KEEP\n";
      return 2;
    }
    const long a_max = parse_long(argv[1], "A_MAX");
    const long b_max = parse_long(argv[2], "B_MAX");
    const long rank_keep = parse_long(argv[3], "RANK_KEEP");
    const long power_keep = parse_long(argv[4], "POWER_KEEP");
    const long composite_keep = parse_long(argv[5], "COMPOSITE_KEEP");
    const long per_denominator_keep = parse_long(argv[6], "PER_DENOMINATOR_KEEP");
    if (a_max < 1 || b_max < 1 || rank_keep < 1 || power_keep < 1 ||
        composite_keep < 1 || per_denominator_keep < 1 ||
        a_max > 2000000000L || b_max > 2000000000L) {
      throw std::invalid_argument("all bounds must be positive 32-bit integers");
    }

    const auto discovery = read_band();
    const auto held = read_band();

    Heap by_rank(Better{Metric::Rank});
    Heap by_power(Better{Metric::Power});
    Heap by_composite(Better{Metric::Composite});
    std::vector<Candidate> denominator_frontier;
    std::uint64_t primitive_count = 0;

    for (int denominator = 1; denominator <= b_max; ++denominator) {
      auto discovery_state = initial_state(discovery, denominator);
      auto held_state = initial_state(held, denominator);
      Heap local(Better{Metric::Composite});
      for (int numerator = 0; numerator <= a_max; ++numerator) {
        if (std::gcd(numerator, denominator) == 1) {
          ++primitive_count;
          const auto discovery_scores = score_band(discovery, discovery_state);
          const auto held_scores = score_band(held, held_state);
          const Candidate candidate{numerator, denominator,
                                    discovery_scores.first,
                                    discovery_scores.second,
                                    held_scores.first, held_scores.second};
          retain(by_rank, candidate, static_cast<std::size_t>(rank_keep),
                 Metric::Rank);
          retain(by_power, candidate, static_cast<std::size_t>(power_keep),
                 Metric::Power);
          retain(by_composite, candidate,
                 static_cast<std::size_t>(composite_keep), Metric::Composite);
          retain(local, candidate,
                 static_cast<std::size_t>(per_denominator_keep),
                 Metric::Composite);
        }
        advance(discovery, discovery_state);
        advance(held, held_state);
      }
      auto rows = drain(local, Metric::Composite);
      denominator_frontier.insert(denominator_frontier.end(), rows.begin(), rows.end());
    }

    auto rank_rows = drain(by_rank, Metric::Rank);
    auto power_rows = drain(by_power, Metric::Power);
    auto composite_rows = drain(by_composite, Metric::Composite);
    std::unordered_map<std::uint64_t, Candidate> unique;
    auto insert = [&unique](const Candidate& candidate) {
      const auto key = (static_cast<std::uint64_t>(candidate.denominator) << 32) |
                       static_cast<std::uint32_t>(candidate.numerator);
      unique.emplace(key, candidate);
    };
    for (const auto& row : rank_rows) insert(row);
    for (const auto& row : power_rows) insert(row);
    for (const auto& row : composite_rows) insert(row);
    for (const auto& row : denominator_frontier) insert(row);
    std::vector<Candidate> output;
    output.reserve(unique.size());
    for (const auto& item : unique) output.push_back(item.second);
    std::sort(output.begin(), output.end(), [](const Candidate& left,
                                               const Candidate& right) {
      if (quality(left, Metric::Composite) != quality(right, Metric::Composite))
        return quality(left, Metric::Composite) > quality(right, Metric::Composite);
      return std::tie(left.denominator, left.numerator) <
             std::tie(right.denominator, right.numerator);
    });

    std::cout << "SUMMARY\t" << primitive_count << '\t' << rank_rows.size()
              << '\t' << power_rows.size() << '\t' << composite_rows.size()
              << '\t' << denominator_frontier.size() << '\t' << output.size()
              << '\n';
    for (const auto& row : output) {
      std::cout << "ROW\t" << row.numerator << '\t' << row.denominator << '\t'
                << row.discovery_rank << '\t' << row.discovery_power << '\t'
                << row.held_rank << '\t' << row.held_power << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
