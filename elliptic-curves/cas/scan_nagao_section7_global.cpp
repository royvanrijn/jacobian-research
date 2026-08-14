// Exhaustive, residue-score-only scanner for Nagao's 1994 section-7 family.
//
// The Python driver writes two disjoint prime-band lookup tables to stdin.
// This helper enumerates every positive primitive (a,b) in the declared box.
// Candidate retention uses the first (training) band only.  The second band
// is merely carried forward so the driver can perform a genuinely held-out
// validation after the complete scan.  No curve, conductor, or known-record
// information is compiled into this helper.

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

struct Table {
  int prime = 0;
  std::vector<std::int64_t> weights;
};

struct Candidate {
  int numerator = 0;
  int denominator = 0;
  std::int64_t training_score = 0;
  std::int64_t validation_score = 0;
};

// A larger tuple is a better training-only candidate.  The deterministic
// height/numerator tie breaks contain no held-out information.
auto quality(const Candidate& candidate) {
  return std::tuple(candidate.training_score,
                    -std::max(candidate.numerator, candidate.denominator),
                    -candidate.numerator, -candidate.denominator);
}

struct Better {
  bool operator()(const Candidate& left, const Candidate& right) const {
    return quality(left) > quality(right);
  }
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
  if (!(std::cin >> count) || count < 0) {
    throw std::runtime_error("invalid prime-band header");
  }
  std::vector<Table> result;
  result.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index) {
    Table table;
    if (!(std::cin >> table.prime) || table.prime < 2) {
      throw std::runtime_error("invalid score-table prime");
    }
    table.weights.resize(static_cast<std::size_t>(table.prime + 1));
    for (auto& value : table.weights) {
      if (!(std::cin >> value)) {
        throw std::runtime_error("truncated score table");
      }
    }
    result.push_back(std::move(table));
  }
  return result;
}

int inverse_mod(int value, int prime) {
  int old_r = value;
  int r = prime;
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
  old_s %= prime;
  if (old_s < 0) old_s += prime;
  return old_s;
}

void retain(std::priority_queue<Candidate, std::vector<Candidate>, Better>& heap,
            const Candidate& candidate, std::size_t limit) {
  if (limit == 0) return;
  if (heap.size() < limit) {
    heap.push(candidate);
  } else if (quality(candidate) > quality(heap.top())) {
    heap.pop();
    heap.push(candidate);
  }
}

std::vector<Candidate> drain(
    std::priority_queue<Candidate, std::vector<Candidate>, Better>& heap) {
  std::vector<Candidate> result;
  result.reserve(heap.size());
  while (!heap.empty()) {
    result.push_back(heap.top());
    heap.pop();
  }
  std::sort(result.begin(), result.end(), [](const Candidate& left,
                                              const Candidate& right) {
    return quality(left) > quality(right);
  });
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 5) {
      std::cerr << "usage: scan_nagao_section7_global "
                   "A_MAX B_MAX GLOBAL_KEEP PER_DENOMINATOR_KEEP\n";
      return 2;
    }
    const long a_max = parse_long(argv[1], "A_MAX");
    const long b_max = parse_long(argv[2], "B_MAX");
    const long global_keep = parse_long(argv[3], "GLOBAL_KEEP");
    const long per_denominator_keep = parse_long(argv[4], "PER_DENOMINATOR_KEEP");
    if (a_max < 1 || b_max < 1 || global_keep < 1 ||
        per_denominator_keep < 1 || a_max > 2000000000L ||
        b_max > 2000000000L) {
      throw std::invalid_argument("all bounds must be positive 32-bit integers");
    }

    const auto training = read_band();
    const auto validation = read_band();
    if (training.empty() || validation.empty()) {
      throw std::runtime_error("both disjoint score bands must be nonempty");
    }

    std::priority_queue<Candidate, std::vector<Candidate>, Better> global;
    std::vector<Candidate> denominator_frontier;
    std::uint64_t primitive_count = 0;

    for (int denominator = 1; denominator <= b_max; ++denominator) {
      std::vector<int> train_indices;
      std::vector<int> train_steps;
      std::vector<int> validation_indices;
      std::vector<int> validation_steps;
      train_indices.reserve(training.size());
      train_steps.reserve(training.size());
      validation_indices.reserve(validation.size());
      validation_steps.reserve(validation.size());
      for (const auto& table : training) {
        if (denominator % table.prime == 0) {
          train_indices.push_back(table.prime);
          train_steps.push_back(0);
        } else {
          const int step = inverse_mod(denominator % table.prime, table.prime);
          train_indices.push_back(step);  // numerator starts at one
          train_steps.push_back(step);
        }
      }
      for (const auto& table : validation) {
        if (denominator % table.prime == 0) {
          validation_indices.push_back(table.prime);
          validation_steps.push_back(0);
        } else {
          const int step = inverse_mod(denominator % table.prime, table.prime);
          validation_indices.push_back(step);
          validation_steps.push_back(step);
        }
      }

      std::priority_queue<Candidate, std::vector<Candidate>, Better> local;
      for (int numerator = 1; numerator <= a_max; ++numerator) {
        if (std::gcd(numerator, denominator) == 1) {
          ++primitive_count;
          std::int64_t train_score = 0;
          std::int64_t validation_score = 0;
          for (std::size_t index = 0; index < training.size(); ++index) {
            train_score += training[index].weights[
                static_cast<std::size_t>(train_indices[index])];
          }
          for (std::size_t index = 0; index < validation.size(); ++index) {
            validation_score += validation[index].weights[
                static_cast<std::size_t>(validation_indices[index])];
          }
          const Candidate candidate{numerator, denominator, train_score,
                                    validation_score};
          retain(global, candidate, static_cast<std::size_t>(global_keep));
          retain(local, candidate,
                 static_cast<std::size_t>(per_denominator_keep));
        }

        for (std::size_t index = 0; index < training.size(); ++index) {
          if (train_steps[index] != 0) {
            train_indices[index] += train_steps[index];
            if (train_indices[index] >= training[index].prime) {
              train_indices[index] -= training[index].prime;
            }
          }
        }
        for (std::size_t index = 0; index < validation.size(); ++index) {
          if (validation_steps[index] != 0) {
            validation_indices[index] += validation_steps[index];
            if (validation_indices[index] >= validation[index].prime) {
              validation_indices[index] -= validation[index].prime;
            }
          }
        }
      }
      auto local_values = drain(local);
      denominator_frontier.insert(denominator_frontier.end(),
                                  local_values.begin(), local_values.end());
    }

    auto global_values = drain(global);
    std::unordered_map<std::uint64_t, Candidate> unique;
    auto insert = [&unique](const Candidate& candidate) {
      const auto key = (static_cast<std::uint64_t>(candidate.denominator) << 32) |
                       static_cast<std::uint32_t>(candidate.numerator);
      unique.emplace(key, candidate);
    };
    for (const auto& candidate : global_values) insert(candidate);
    for (const auto& candidate : denominator_frontier) insert(candidate);
    std::vector<Candidate> output;
    output.reserve(unique.size());
    for (const auto& item : unique) output.push_back(item.second);
    std::sort(output.begin(), output.end(), [](const Candidate& left,
                                               const Candidate& right) {
      if (quality(left) != quality(right)) return quality(left) > quality(right);
      return std::tie(left.denominator, left.numerator) <
             std::tie(right.denominator, right.numerator);
    });

    std::cout << "SUMMARY\t" << primitive_count << '\t' << global_values.size()
              << '\t' << denominator_frontier.size() << '\t' << output.size()
              << '\n';
    for (const auto& candidate : output) {
      std::cout << "ROW\t" << candidate.numerator << '\t'
                << candidate.denominator << '\t' << candidate.training_score
                << '\t' << candidate.validation_score << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
