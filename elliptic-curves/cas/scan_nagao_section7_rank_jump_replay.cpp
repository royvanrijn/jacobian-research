// Complete score-rank replay for the certified Nagao section-7 rank-20 fibre.
//
// Python supplies exact integer lookup tables for the two frozen prime bands.
// This helper enumerates every positive primitive T=a/b in the archived box
// and counts the target's position.  It retains no candidates and sees no
// Mordell--Weil labels.

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

struct Table {
  int prime = 0;
  std::vector<std::int64_t> weights;
};

struct BandState {
  std::vector<int> indices;
  std::vector<int> steps;
};

long parse_long(const char* text, const char* label) {
  char* end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (end == text || *end != '\0')
    throw std::invalid_argument(std::string("invalid ") + label);
  return value;
}

std::vector<Table> read_band() {
  int count = 0;
  if (!(std::cin >> count) || count < 1)
    throw std::runtime_error("invalid or empty prime-band header");
  std::vector<Table> result;
  result.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index) {
    Table table;
    if (!(std::cin >> table.prime) || table.prime < 2)
      throw std::runtime_error("invalid table prime");
    table.weights.resize(static_cast<std::size_t>(table.prime + 1));
    for (auto& value : table.weights)
      if (!(std::cin >> value))
        throw std::runtime_error("truncated score table");
    result.push_back(std::move(table));
  }
  return result;
}

int inverse_mod(int value, int modulus) {
  int old_r = value, r = modulus, old_s = 1, s = 0;
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

BandState initial_state(const std::vector<Table>& tables, int denominator) {
  BandState state;
  for (const auto& table : tables) {
    if (denominator % table.prime == 0) {
      state.indices.push_back(table.prime);
      state.steps.push_back(0);
    } else {
      const int step = inverse_mod(denominator % table.prime, table.prime);
      state.indices.push_back(step);  // numerator starts at one
      state.steps.push_back(step);
    }
  }
  return state;
}

std::int64_t score(const std::vector<Table>& tables, const BandState& state) {
  std::int64_t total = 0;
  for (std::size_t index = 0; index < tables.size(); ++index)
    total += tables[index].weights[state.indices[index]];
  return total;
}

void advance(const std::vector<Table>& tables, BandState& state) {
  for (std::size_t index = 0; index < tables.size(); ++index) {
    if (state.steps[index] == 0) continue;
    state.indices[index] += state.steps[index];
    if (state.indices[index] >= tables[index].prime)
      state.indices[index] -= tables[index].prime;
  }
}

std::int64_t score_parameter(int numerator, int denominator,
                             const std::vector<Table>& tables) {
  std::int64_t total = 0;
  for (const auto& table : tables) {
    if (denominator % table.prime == 0) continue;
    const int residue = static_cast<int>(
        static_cast<std::int64_t>(numerator % table.prime) *
        inverse_mod(denominator % table.prime, table.prime) % table.prime);
    total += table.weights[residue];
  }
  return total;
}

auto quality(std::int64_t value, int numerator, int denominator) {
  return std::tuple(value, -std::max(numerator, denominator), -numerator,
                    -denominator);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 6) {
      std::cerr << "usage: scan A_MAX B_MAX TARGET_A TARGET_B EXPECTED_COUNT\n";
      return 2;
    }
    const int a_max = static_cast<int>(parse_long(argv[1], "A_MAX"));
    const int b_max = static_cast<int>(parse_long(argv[2], "B_MAX"));
    const int target_a = static_cast<int>(parse_long(argv[3], "TARGET_A"));
    const int target_b = static_cast<int>(parse_long(argv[4], "TARGET_B"));
    const std::uint64_t expected =
        static_cast<std::uint64_t>(parse_long(argv[5], "EXPECTED_COUNT"));
    if (a_max < 1 || b_max < 1 || target_a < 1 || target_a > a_max ||
        target_b < 1 || target_b > b_max || std::gcd(target_a, target_b) != 1)
      throw std::invalid_argument("invalid box or target");

    const auto training = read_band();
    const auto validation = read_band();
    const auto target_training = score_parameter(target_a, target_b, training);
    const auto target_validation = score_parameter(target_a, target_b, validation);
    std::uint64_t enumerated = 0, ahead_training = 0, ahead_validation = 0;
    std::uint64_t equal_training = 0, equal_validation = 0;
    for (int denominator = 1; denominator <= b_max; ++denominator) {
      auto training_state = initial_state(training, denominator);
      auto validation_state = initial_state(validation, denominator);
      for (int numerator = 1; numerator <= a_max; ++numerator) {
        if (std::gcd(numerator, denominator) == 1) {
          ++enumerated;
          const auto training_value = score(training, training_state);
          const auto validation_value = score(validation, validation_state);
          if (training_value == target_training) ++equal_training;
          if (validation_value == target_validation) ++equal_validation;
          if (quality(training_value, numerator, denominator) >
              quality(target_training, target_a, target_b))
            ++ahead_training;
          if (quality(validation_value, numerator, denominator) >
              quality(target_validation, target_a, target_b))
            ++ahead_validation;
        }
        advance(training, training_state);
        advance(validation, validation_state);
      }
    }
    if (enumerated != expected)
      throw std::runtime_error("enumerated count differs from frozen population");
    std::cout << "SUMMARY\t" << enumerated << '\n';
    std::cout << "TARGET\t" << target_a << '\t' << target_b << '\t'
              << target_training << '\t' << target_validation << '\t'
              << ahead_training + 1 << '\t' << equal_training << '\t'
              << ahead_validation + 1 << '\t' << equal_validation << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
