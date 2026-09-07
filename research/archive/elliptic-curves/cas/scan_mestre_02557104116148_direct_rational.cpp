// Exact-local direct rational-T scanner for the Nagao rank-13 root tuple.
//
// The six-root Mestre family has roots (0,25,57,104,116,148) and primitive
// Jacobian y^2=x^3+A(T)x+B(T), where A and B are even of degrees 8 and 12.
// For primitive T=a/b>0 this program evaluates the weighted homogeneous
// coefficients at fresh discovery primes 601..701.  Only the fixed discovery
// tail is then read at the disjoint held band 709..809.  A caller-supplied
// exact exclusion file contains one reduced positive numerator/denominator
// pair per line; it keeps this direct-T population disjoint from every
// discoverable earlier Nagao rank-13 specialization.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

struct LocalSymbol {
  bool good = false;
  int trace = 0;
  std::int64_t score_units = 0;
};

struct PrimeTable {
  int prime = 0;
  std::vector<LocalSymbol> symbols;
};

struct Candidate {
  int numerator = 0;
  int denominator = 1;
  std::int64_t discovery_score = 0;
  std::int64_t held_score = 0;
  int discovery_good = 0;
  int held_good = 0;
};

static constexpr std::array<int, 17> DISCOVERY_PRIMES{
    601, 607, 613, 617, 619, 631, 641, 643, 647,
    653, 659, 661, 673, 677, 683, 691, 701};
static constexpr std::array<int, 14> HELD_PRIMES{
    709, 719, 727, 733, 739, 743, 751,
    757, 761, 769, 773, 787, 797, 809};

static const std::array<std::string, 9> A_COEFFICIENTS{{
    "-27546462334108146267", "0", "-1530834958134000", "0",
    "3106634557536", "0", "-359251200", "0", "-34992"}};

static const std::array<std::string, 13> B_COEFFICIENTS{{
    "55624396621360883431446459126", "0",
    "5508337118494541464141200", "0",
    "-13510101631102695979884", "0",
    "4296246530145998400", "0",
    "-634269830133888", "0", "38799129600", "0", "2519424"}};

static int decimal_mod(const std::string& text, int prime) {
  const bool negative = !text.empty() && text[0] == '-';
  int answer = 0;
  for (std::size_t index = negative ? 1 : 0; index < text.size(); ++index) {
    answer = (answer * 10 + (text[index] - '0')) % prime;
  }
  return negative && answer ? prime - answer : answer;
}

static int multiply_mod(std::int64_t left, std::int64_t right, int prime) {
  return static_cast<int>((left * right) % prime);
}

static int power_mod(int base, int exponent, int prime) {
  int answer = 1;
  while (exponent) {
    if (exponent & 1) answer = multiply_mod(answer, base, prime);
    base = multiply_mod(base, base, prime);
    exponent >>= 1;
  }
  return answer;
}

template <std::size_t N>
static int evaluate_finite(
    const std::array<std::string, N>& coefficients, int value, int prime) {
  int answer = 0;
  for (std::size_t offset = 0; offset < N; ++offset) {
    const std::size_t index = N - 1 - offset;
    answer = (multiply_mod(answer, value, prime)
              + decimal_mod(coefficients[index], prime)) % prime;
  }
  return answer;
}

static int trace_of_frobenius(int coefficient_a, int coefficient_b, int prime) {
  int character_sum = 0;
  for (int x = 0; x < prime; ++x) {
    const int rhs = (multiply_mod(multiply_mod(x, x, prime), x, prime)
                     + multiply_mod(coefficient_a, x, prime) + coefficient_b)
                    % prime;
    if (rhs == 0) continue;
    const int symbol = power_mod(rhs, (prime - 1) / 2, prime);
    character_sum += symbol == 1 ? 1 : -1;
  }
  return -character_sum;
}

static LocalSymbol local_symbol(int residue, bool infinity, int prime) {
  const int coefficient_a = infinity
                                ? decimal_mod(A_COEFFICIENTS.back(), prime)
                                : evaluate_finite(A_COEFFICIENTS, residue, prime);
  const int coefficient_b = infinity
                                ? decimal_mod(B_COEFFICIENTS.back(), prime)
                                : evaluate_finite(B_COEFFICIENTS, residue, prime);
  const int discriminant_core =
      (4 * multiply_mod(multiply_mod(coefficient_a, coefficient_a, prime),
                        coefficient_a, prime)
       + 27 * multiply_mod(coefficient_b, coefficient_b, prime)) % prime;
  if (discriminant_core == 0) return {};
  const int trace = trace_of_frobenius(coefficient_a, coefficient_b, prime);
  const int group_order = prime + 1 - trace;
  const double score =
      (static_cast<double>(2 - trace) / static_cast<double>(group_order))
      * std::log(static_cast<double>(prime));
  return {true, trace,
          static_cast<std::int64_t>(std::llround(score * 1.0e12))};
}

template <std::size_t N>
static std::vector<PrimeTable> build_tables(const std::array<int, N>& primes) {
  std::vector<PrimeTable> answer;
  answer.reserve(N);
  for (const int prime : primes) {
    PrimeTable table;
    table.prime = prime;
    table.symbols.reserve(prime + 1);
    for (int residue = 0; residue < prime; ++residue) {
      table.symbols.push_back(local_symbol(residue, false, prime));
    }
    table.symbols.push_back(local_symbol(0, true, prime));
    answer.push_back(std::move(table));
  }
  return answer;
}

static bool better(const Candidate& left, const Candidate& right) {
  if (left.discovery_score != right.discovery_score)
    return left.discovery_score > right.discovery_score;
  if (left.discovery_good != right.discovery_good)
    return left.discovery_good > right.discovery_good;
  if (left.denominator != right.denominator)
    return left.denominator < right.denominator;
  return left.numerator < right.numerator;
}

struct BetterComparator {
  bool operator()(const Candidate& left, const Candidate& right) const {
    return better(left, right);
  }
};

static void add_band_score(
    Candidate& candidate, const std::vector<PrimeTable>& tables, bool held) {
  std::int64_t score = 0;
  int good = 0;
  for (const PrimeTable& table : tables) {
    const int prime = table.prime;
    int symbol_index = prime;
    if (candidate.denominator % prime != 0) {
      const int inverse = power_mod(candidate.denominator % prime, prime - 2, prime);
      symbol_index = multiply_mod(candidate.numerator % prime, inverse, prime);
    }
    const LocalSymbol& symbol = table.symbols[symbol_index];
    if (symbol.good) {
      score += symbol.score_units;
      ++good;
    }
  }
  if (held) {
    candidate.held_score = score;
    candidate.held_good = good;
  } else {
    candidate.discovery_score = score;
    candidate.discovery_good = good;
  }
}

static std::uint64_t table_digest(const std::vector<PrimeTable>& tables) {
  std::uint64_t digest = 1469598103934665603ULL;
  auto mix = [&digest](std::uint64_t value) {
    for (int offset = 0; offset < 8; ++offset) {
      digest ^= (value >> (8 * offset)) & 255ULL;
      digest *= 1099511628211ULL;
    }
  };
  for (const PrimeTable& table : tables) {
    mix(static_cast<std::uint64_t>(table.prime));
    for (const LocalSymbol& symbol : table.symbols) {
      mix(static_cast<std::uint64_t>(symbol.good));
      mix(static_cast<std::uint64_t>(static_cast<std::int64_t>(symbol.trace)));
    }
  }
  return digest;
}

static std::string score_text(std::int64_t units) {
  const bool negative = units < 0;
  const std::uint64_t absolute = negative ? -units : units;
  std::string answer = negative ? "-" : "";
  answer += std::to_string(absolute / 1000000000000ULL) + ".";
  std::string tail = std::to_string(absolute % 1000000000000ULL);
  answer += std::string(12 - tail.size(), '0') + tail;
  return answer;
}

static std::uint64_t key(int numerator, int denominator) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(numerator)) << 32)
         | static_cast<std::uint32_t>(denominator);
}

static std::unordered_set<std::uint64_t> load_exclusions(
    const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("could not open exclusion file");
  std::unordered_set<std::uint64_t> answer;
  int numerator = 0;
  int denominator = 0;
  while (input >> numerator >> denominator) {
    if (numerator <= 0 || denominator <= 0
        || std::gcd(numerator, denominator) != 1) {
      throw std::runtime_error("malformed non-reduced exclusion");
    }
    answer.insert(key(numerator, denominator));
  }
  if (!input.eof()) throw std::runtime_error("malformed exclusion file");
  return answer;
}

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: scan_mestre_02557104116148_direct_rational "
                 "NUM_BOUND DEN_BOUND KEEP EXCLUSIONS\n";
    return 2;
  }
  const int numerator_bound = std::atoi(argv[1]);
  const int denominator_bound = std::atoi(argv[2]);
  const int keep = std::atoi(argv[3]);
  if (numerator_bound < 1 || numerator_bound > 100000
      || denominator_bound < 1 || denominator_bound > 5000 || keep < 1
      || keep > 20000) {
    std::cerr << "require NUM<=100000 DEN<=5000 KEEP<=20000\n";
    return 2;
  }
  const auto exclusions = load_exclusions(argv[4]);
  const std::vector<PrimeTable> discovery = build_tables(DISCOVERY_PRIMES);
  const std::vector<PrimeTable> held = build_tables(HELD_PRIMES);
  std::priority_queue<Candidate, std::vector<Candidate>, BetterComparator> heap;
  std::uint64_t primitive_population = 0;
  std::uint64_t prior_excluded = 0;
  std::uint64_t evaluated_population = 0;

  for (int denominator = 1; denominator <= denominator_bound; ++denominator) {
    std::vector<int> finite_multiplier;
    finite_multiplier.reserve(discovery.size());
    for (const PrimeTable& table : discovery) {
      finite_multiplier.push_back(
          denominator % table.prime == 0
              ? -1
              : power_mod(denominator % table.prime, table.prime - 2, table.prime));
    }
    std::vector<int> residues(discovery.size(), 0);
    for (int numerator = 1; numerator <= numerator_bound; ++numerator) {
      for (std::size_t index = 0; index < discovery.size(); ++index) {
        if (finite_multiplier[index] >= 0) {
          residues[index] += finite_multiplier[index];
          if (residues[index] >= discovery[index].prime)
            residues[index] -= discovery[index].prime;
        }
      }
      if (std::gcd(numerator, denominator) != 1) continue;
      ++primitive_population;
      if (exclusions.count(key(numerator, denominator))) {
        ++prior_excluded;
        continue;
      }
      ++evaluated_population;
      Candidate candidate;
      candidate.numerator = numerator;
      candidate.denominator = denominator;
      for (std::size_t index = 0; index < discovery.size(); ++index) {
        const PrimeTable& table = discovery[index];
        const int symbol_index = finite_multiplier[index] < 0
                                     ? table.prime
                                     : residues[index];
        const LocalSymbol& symbol = table.symbols[symbol_index];
        if (symbol.good) {
          candidate.discovery_score += symbol.score_units;
          ++candidate.discovery_good;
        }
      }
      if (static_cast<int>(heap.size()) < keep) {
        heap.push(candidate);
      } else if (better(candidate, heap.top())) {
        heap.pop();
        heap.push(candidate);
      }
    }
  }

  std::vector<Candidate> retained;
  retained.reserve(heap.size());
  while (!heap.empty()) {
    retained.push_back(heap.top());
    heap.pop();
  }
  std::sort(retained.begin(), retained.end(), better);
  for (Candidate& candidate : retained) add_band_score(candidate, held, true);

  Candidate calibration;
  calibration.numerator = 1;
  calibration.denominator = 1;
  add_band_score(calibration, discovery, false);
  add_band_score(calibration, held, true);

  std::cout << "MESTRE_02557104116148_DIRECT_SCAN_V1\nD";
  for (const int prime : DISCOVERY_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nH";
  for (const int prime : HELD_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nL " << table_digest(discovery) << ' ' << table_digest(held)
            << "\nK 1 1 " << score_text(calibration.discovery_score) << ' '
            << score_text(calibration.held_score) << ' '
            << calibration.discovery_good << ' ' << calibration.held_good << '\n';
  for (const Candidate& candidate : retained) {
    std::cout << "C " << candidate.numerator << ' ' << candidate.denominator << ' '
              << score_text(candidate.discovery_score) << ' '
              << score_text(candidate.held_score) << ' '
              << candidate.discovery_good << ' ' << candidate.held_good << '\n';
  }
  std::cout << "S " << numerator_bound << ' ' << denominator_bound << ' ' << keep
            << ' ' << primitive_population << ' ' << prior_excluded << ' '
            << evaluated_population << ' ' << retained.size() << ' '
            << exclusions.size() << '\n';
  return 0;
}
