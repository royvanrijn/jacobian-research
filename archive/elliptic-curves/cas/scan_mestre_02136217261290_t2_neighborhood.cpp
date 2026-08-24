// Exact-local rational scanner around T=2 for roots (0,2,136,217,261,290).
//
// The short Jacobian is y^2=x^3+A(T)x+B(T), with even A,B of degrees 8,12.
// For T=a/b the local model uses the weighted homogeneous values modulo p,
// including the exact b=0 projective symbol.  The discovery and held-out
// prime bands are disjoint and were not used to select the max-root-300 T=2
// calibration fiber.  Every trace and good-reduction decision is exact; only
// the logged score is quantized to 10^-12.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
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

static constexpr std::array<int, 12> DISCOVERY_PRIMES{
    587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647};
static constexpr std::array<int, 12> HELD_PRIMES{
    653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733};

static const std::array<std::string, 9> A_COEFFICIENTS{{
    "-542564766112960201029552", "0", "169665434230038056352", "0",
    "-17668537074117936", "0", "742147660800", "0", "-21676032"}};
static const std::array<std::string, 13> B_COEFFICIENTS{{
    "153797356987326659323597489852634496", "0",
    "-71778413690052019619903289482112", "0",
    "12434077250691567747253858944", "0",
    "-867706293306445539767424", "0", "3733083805979308032", "0",
    "1994892912230400", "0", "-38843449344"}};

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
    for (int residue = 0; residue < prime; ++residue)
      table.symbols.push_back(local_symbol(residue, false, prime));
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

static std::vector<int> denominator_multipliers(
    int denominator, const std::vector<PrimeTable>& tables) {
  std::vector<int> answer;
  answer.reserve(tables.size());
  for (const PrimeTable& table : tables) {
    answer.push_back(
        denominator % table.prime == 0
            ? -1
            : power_mod(denominator % table.prime, table.prime - 2, table.prime));
  }
  return answer;
}

static void add_discovery_score(
    Candidate& candidate,
    const std::vector<PrimeTable>& tables,
    const std::vector<int>& multipliers) {
  for (std::size_t index = 0; index < tables.size(); ++index) {
    const PrimeTable& table = tables[index];
    const int symbol_index = multipliers[index] < 0
                                 ? table.prime
                                 : multiply_mod(
                                       candidate.numerator % table.prime,
                                       multipliers[index], table.prime);
    const LocalSymbol& symbol = table.symbols[symbol_index];
    if (symbol.good) {
      candidate.discovery_score += symbol.score_units;
      ++candidate.discovery_good;
    }
  }
}

static void add_held_score(
    Candidate& candidate, const std::vector<PrimeTable>& tables) {
  const std::vector<int> multipliers =
      denominator_multipliers(candidate.denominator, tables);
  for (std::size_t index = 0; index < tables.size(); ++index) {
    const PrimeTable& table = tables[index];
    const int symbol_index = multipliers[index] < 0
                                 ? table.prime
                                 : multiply_mod(
                                       candidate.numerator % table.prime,
                                       multipliers[index], table.prime);
    const LocalSymbol& symbol = table.symbols[symbol_index];
    if (symbol.good) {
      candidate.held_score += symbol.score_units;
      ++candidate.held_good;
    }
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

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: scan_mestre_02136217261290_t2_neighborhood "
                 "MODE DEN_BOUND KEEP\n";
    return 2;
  }
  const std::string mode = argv[1];
  const int denominator_bound = std::atoi(argv[2]);
  const int keep = std::atoi(argv[3]);
  if ((mode != "near" && mode != "ordinary") || denominator_bound < 2
      || denominator_bound > (mode == "near" ? 20000 : 2000)
      || keep < 1 || keep > 20000) {
    std::cerr << "mode near requires DEN<=20000; ordinary requires DEN<=2000; "
                 "KEEP<=20000\n";
    return 2;
  }

  const std::vector<PrimeTable> discovery = build_tables(DISCOVERY_PRIMES);
  const std::vector<PrimeTable> held = build_tables(HELD_PRIMES);
  std::priority_queue<Candidate, std::vector<Candidate>, BetterComparator> heap;
  std::uint64_t proposed_population = 0;
  std::uint64_t primitive_population = 0;

  auto consider = [&](int numerator, int denominator,
                      const std::vector<int>& multipliers) {
    ++proposed_population;
    if (numerator <= 0 || std::gcd(numerator, denominator) != 1) return;
    ++primitive_population;
    Candidate candidate;
    candidate.numerator = numerator;
    candidate.denominator = denominator;
    add_discovery_score(candidate, discovery, multipliers);
    if (static_cast<int>(heap.size()) < keep) {
      heap.push(candidate);
    } else if (better(candidate, heap.top())) {
      heap.pop();
      heap.push(candidate);
    }
  };

  for (int denominator = 2; denominator <= denominator_bound; ++denominator) {
    const std::vector<int> multipliers =
        denominator_multipliers(denominator, discovery);
    if (mode == "near") {
      for (int offset = -32; offset <= 32; ++offset) {
        if (offset != 0) consider(2 * denominator + offset, denominator, multipliers);
      }
    } else {
      const int first = (3 * denominator + 1) / 2;
      const int last = (5 * denominator) / 2;
      for (int numerator = first; numerator <= last; ++numerator) {
        if (std::abs(numerator - 2 * denominator) >= 33)
          consider(numerator, denominator, multipliers);
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
  for (Candidate& candidate : retained) add_held_score(candidate, held);

  Candidate calibration;
  calibration.numerator = 2;
  calibration.denominator = 1;
  const std::vector<int> calibration_discovery =
      denominator_multipliers(1, discovery);
  add_discovery_score(calibration, discovery, calibration_discovery);
  add_held_score(calibration, held);

  std::cout << "MESTRE_02136217261290_T2_NEIGHBORHOOD_SCAN_V1\nD";
  for (const int prime : DISCOVERY_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nH";
  for (const int prime : HELD_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nL " << table_digest(discovery) << ' ' << table_digest(held)
            << "\nA " << score_text(calibration.discovery_score) << ' '
            << score_text(calibration.held_score) << ' '
            << calibration.discovery_good << ' ' << calibration.held_good << '\n';
  for (const Candidate& candidate : retained) {
    std::cout << "C " << candidate.numerator << ' ' << candidate.denominator << ' '
              << score_text(candidate.discovery_score) << ' '
              << score_text(candidate.held_score) << ' '
              << candidate.discovery_good << ' ' << candidate.held_good << '\n';
  }
  std::cout << "S " << mode << ' ' << denominator_bound << ' ' << keep << ' '
            << proposed_population << ' ' << primitive_population << ' '
            << retained.size() << '\n';
  return 0;
}
