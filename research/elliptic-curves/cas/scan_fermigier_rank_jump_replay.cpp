// Exhaustive score-rank replay for two certified Fermigier escape anchors.
//
// Python supplies exact integer lookup tables for two disjoint prime bands.
// This program enumerates the declared primitive projective box and counts
// how many parameters sort ahead of each anchor under four already-frozen
// scores.  It retains no candidates and uses no Mordell--Weil labels.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr std::int64_t kCompositePowerDivisor = 16;

struct Table {
  int prime = 0;
  int power_modulus = 0;
  std::vector<std::int64_t> rank_weights;
  std::vector<std::int64_t> power_weights;
};

struct Scores {
  std::int64_t discovery_rank = 0;
  std::int64_t discovery_power = 0;
  std::int64_t held_rank = 0;
  std::int64_t held_power = 0;
};

struct Target {
  int numerator = 0;
  int denominator = 0;
  Scores scores;
  std::uint64_t ahead_discovery_rank = 0;
  std::uint64_t ahead_discovery_composite = 0;
  std::uint64_t ahead_held_rank = 0;
  std::uint64_t ahead_held_composite = 0;
  std::uint64_t equal_discovery_rank = 0;
  std::uint64_t equal_discovery_composite = 0;
  std::uint64_t equal_held_rank = 0;
  std::uint64_t equal_held_composite = 0;
};

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

struct BandState {
  std::vector<int> rank_indices;
  std::vector<int> rank_steps;
  std::vector<int> power_indices;
  std::vector<int> power_steps;
};

BandState initial_state(const std::vector<Table>& tables, int denominator) {
  BandState state;
  for (const auto& table : tables) {
    if (denominator % table.prime == 0) {
      state.rank_indices.push_back(table.prime);
      state.rank_steps.push_back(0);
      state.power_indices.push_back(table.power_modulus);
      state.power_steps.push_back(0);
    } else {
      state.rank_indices.push_back(0);
      state.rank_steps.push_back(
          inverse_mod(denominator % table.prime, table.prime));
      state.power_indices.push_back(0);
      state.power_steps.push_back(
          inverse_mod(denominator % table.power_modulus,
                      table.power_modulus));
    }
  }
  return state;
}

std::pair<std::int64_t, std::int64_t> score_band(
    const std::vector<Table>& tables, const BandState& state) {
  std::int64_t rank = 0;
  std::int64_t power = 0;
  for (std::size_t index = 0; index < tables.size(); ++index) {
    rank += tables[index].rank_weights[state.rank_indices[index]];
    power += tables[index].power_weights[state.power_indices[index]];
  }
  return {rank, power};
}

void advance(const std::vector<Table>& tables, BandState& state) {
  for (std::size_t index = 0; index < tables.size(); ++index) {
    if (state.rank_steps[index] == 0) continue;
    state.rank_indices[index] += state.rank_steps[index];
    if (state.rank_indices[index] >= tables[index].prime)
      state.rank_indices[index] -= tables[index].prime;
    state.power_indices[index] += state.power_steps[index];
    if (state.power_indices[index] >= tables[index].power_modulus)
      state.power_indices[index] -= tables[index].power_modulus;
  }
}

Scores score_parameter(int numerator, int denominator,
                       const std::vector<Table>& discovery,
                       const std::vector<Table>& held) {
  const auto score_one = [numerator, denominator](
                             const std::vector<Table>& tables) {
    std::int64_t rank = 0;
    std::int64_t power = 0;
    for (const auto& table : tables) {
      if (denominator % table.prime == 0) continue;
      const int rank_index = static_cast<int>(
          static_cast<std::int64_t>(numerator % table.prime) *
          inverse_mod(denominator % table.prime, table.prime) % table.prime);
      const int power_index = static_cast<int>(
          static_cast<std::int64_t>(numerator % table.power_modulus) *
          inverse_mod(denominator % table.power_modulus,
                      table.power_modulus) %
          table.power_modulus);
      rank += table.rank_weights[rank_index];
      power += table.power_weights[power_index];
    }
    return std::pair(rank, power);
  };
  const auto d = score_one(discovery);
  const auto h = score_one(held);
  return Scores{d.first, d.second, h.first, h.second};
}

std::int64_t composite(std::int64_t rank, std::int64_t power) {
  return rank + power / kCompositePowerDivisor;
}

auto quality(std::int64_t score, int numerator, int denominator) {
  return std::tuple(score, -std::max(numerator, denominator), -denominator,
                    -numerator);
}

void compare_one(std::int64_t score, std::int64_t target_score, int numerator,
                 int denominator, const Target& target, std::uint64_t& ahead,
                 std::uint64_t& equal) {
  if (score == target_score) ++equal;
  if (quality(score, numerator, denominator) >
      quality(target_score, target.numerator, target.denominator))
    ++ahead;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 8) {
      std::cerr << "usage: scan A_MAX B_MAX A1 B1 A2 B2 EXPECTED_COUNT\n";
      return 2;
    }
    const int a_max = static_cast<int>(parse_long(argv[1], "A_MAX"));
    const int b_max = static_cast<int>(parse_long(argv[2], "B_MAX"));
    Target targets[2];
    targets[0].numerator = static_cast<int>(parse_long(argv[3], "A1"));
    targets[0].denominator = static_cast<int>(parse_long(argv[4], "B1"));
    targets[1].numerator = static_cast<int>(parse_long(argv[5], "A2"));
    targets[1].denominator = static_cast<int>(parse_long(argv[6], "B2"));
    const std::uint64_t expected =
        static_cast<std::uint64_t>(parse_long(argv[7], "EXPECTED_COUNT"));
    if (a_max < 1 || b_max < 1) throw std::invalid_argument("invalid box");
    const auto discovery = read_band();
    const auto held = read_band();
    for (auto& target : targets) {
      if (std::gcd(target.numerator, target.denominator) != 1 ||
          target.numerator < 0 || target.numerator > a_max ||
          target.denominator < 1 || target.denominator > b_max)
        throw std::invalid_argument("target is outside the primitive box");
      target.scores = score_parameter(target.numerator, target.denominator,
                                      discovery, held);
    }

    std::uint64_t enumerated = 0;
    for (int denominator = 1; denominator <= b_max; ++denominator) {
      auto d_state = initial_state(discovery, denominator);
      auto h_state = initial_state(held, denominator);
      for (int numerator = 0; numerator <= a_max; ++numerator) {
        if (std::gcd(numerator, denominator) == 1) {
          ++enumerated;
          const auto d = score_band(discovery, d_state);
          const auto h = score_band(held, h_state);
          for (auto& target : targets) {
            compare_one(d.first, target.scores.discovery_rank, numerator,
                        denominator, target, target.ahead_discovery_rank,
                        target.equal_discovery_rank);
            compare_one(composite(d.first, d.second),
                        composite(target.scores.discovery_rank,
                                  target.scores.discovery_power),
                        numerator, denominator, target,
                        target.ahead_discovery_composite,
                        target.equal_discovery_composite);
            compare_one(h.first, target.scores.held_rank, numerator,
                        denominator, target, target.ahead_held_rank,
                        target.equal_held_rank);
            compare_one(composite(h.first, h.second),
                        composite(target.scores.held_rank,
                                  target.scores.held_power),
                        numerator, denominator, target,
                        target.ahead_held_composite,
                        target.equal_held_composite);
          }
        }
        advance(discovery, d_state);
        advance(held, h_state);
      }
    }
    if (enumerated != expected)
      throw std::runtime_error("enumerated count differs from frozen population");
    std::cout << "SUMMARY\t" << enumerated << '\n';
    for (const auto& target : targets) {
      std::cout << "TARGET\t" << target.numerator << '\t' << target.denominator
                << '\t' << target.scores.discovery_rank << '\t'
                << target.scores.discovery_power << '\t'
                << target.scores.held_rank << '\t'
                << target.scores.held_power << '\t'
                << target.ahead_discovery_rank + 1 << '\t'
                << target.equal_discovery_rank << '\t'
                << target.ahead_discovery_composite + 1 << '\t'
                << target.equal_discovery_composite << '\t'
                << target.ahead_held_rank + 1 << '\t'
                << target.equal_held_rank << '\t'
                << target.ahead_held_composite + 1 << '\t'
                << target.equal_held_composite << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
