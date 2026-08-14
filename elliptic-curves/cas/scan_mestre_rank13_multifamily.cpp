// Exact local-trace scanner for a supplied even Mestre Jacobian family.
//
// The coefficient manifest contains the nine ascending coefficients of A(T)
// followed by the thirteen ascending coefficients of B(T), where
// y^2=x^3+A(T)x+B(T).  At T=a/b the scanner evaluates the corresponding
// degree-8/degree-12 homogeneous forms modulo p.  Thus denominators divisible
// by p are handled by the exact projective symbol T=infinity.  All traces and
// good-reduction decisions are exact; only the Nagao-style score is quantized
// to 10^-12 for deterministic ranking.

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

static constexpr std::array<int, 8> DISCOVERY_PRIMES{
    811, 821, 823, 827, 829, 839, 853, 857};
static constexpr std::array<int, 8> HELD_PRIMES{
    859, 863, 877, 881, 883, 887, 907, 911};

struct Coefficients {
  std::array<std::string, 9> a;
  std::array<std::string, 13> b;
};

static Coefficients read_coefficients(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("could not open coefficient manifest");
  std::string header;
  input >> header;
  if (header != "MESTRE_EVEN_COEFFICIENTS_V1")
    throw std::runtime_error("bad coefficient-manifest header");
  Coefficients answer;
  for (std::string& value : answer.a) input >> value;
  for (std::string& value : answer.b) input >> value;
  std::string trailing;
  if (!input || (input >> trailing))
    throw std::runtime_error("malformed coefficient manifest");
  return answer;
}

static int decimal_mod(const std::string& text, int prime) {
  const bool negative = !text.empty() && text[0] == '-';
  int answer = 0;
  for (std::size_t index = negative ? 1 : 0; index < text.size(); ++index)
    answer = (answer * 10 + (text[index] - '0')) % prime;
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

static std::vector<signed char> quadratic_characters(int prime) {
  std::vector<signed char> answer(prime, -1);
  answer[0] = 0;
  for (int value = 1; value < prime; ++value)
    answer[multiply_mod(value, value, prime)] = 1;
  return answer;
}

static int trace_of_frobenius(
    int coefficient_a, int coefficient_b, int prime,
    const std::vector<signed char>& characters) {
  int character_sum = 0;
  for (int x = 0; x < prime; ++x) {
    const int rhs = (multiply_mod(multiply_mod(x, x, prime), x, prime)
                     + multiply_mod(coefficient_a, x, prime) + coefficient_b)
                    % prime;
    character_sum += characters[rhs];
  }
  return -character_sum;
}

static LocalSymbol local_symbol(
    const Coefficients& coefficients, int residue, bool infinity, int prime,
    const std::vector<signed char>& characters) {
  const int coefficient_a = infinity
      ? decimal_mod(coefficients.a.back(), prime)
      : evaluate_finite(coefficients.a, residue, prime);
  const int coefficient_b = infinity
      ? decimal_mod(coefficients.b.back(), prime)
      : evaluate_finite(coefficients.b, residue, prime);
  const int discriminant_core =
      (4 * multiply_mod(multiply_mod(coefficient_a, coefficient_a, prime),
                        coefficient_a, prime)
       + 27 * multiply_mod(coefficient_b, coefficient_b, prime)) % prime;
  if (discriminant_core == 0) return {};
  const int trace = trace_of_frobenius(
      coefficient_a, coefficient_b, prime, characters);
  const int group_order = prime + 1 - trace;
  const double score =
      (static_cast<double>(2 - trace) / static_cast<double>(group_order))
      * std::log(static_cast<double>(prime));
  return {true, trace,
          static_cast<std::int64_t>(std::llround(score * 1.0e12))};
}

template <std::size_t N>
static std::vector<PrimeTable> build_tables(
    const Coefficients& coefficients, const std::array<int, N>& primes) {
  std::vector<PrimeTable> answer;
  answer.reserve(N);
  for (const int prime : primes) {
    const std::vector<signed char> characters = quadratic_characters(prime);
    PrimeTable table;
    table.prime = prime;
    table.symbols.reserve(prime + 1);
    for (int residue = 0; residue < prime; ++residue)
      table.symbols.push_back(
          local_symbol(coefficients, residue, false, prime, characters));
    table.symbols.push_back(
        local_symbol(coefficients, 0, true, prime, characters));
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
      const int inverse = power_mod(
          candidate.denominator % prime, prime - 2, prime);
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

int main(int argc, char** argv) {
  if (argc != 8) {
    std::cerr << "usage: scan_mestre_rank13_multifamily COEFF_FILE FAMILY "
                 "CAL_NUM CAL_DEN NUM_BOUND DEN_BOUND KEEP\n";
    return 2;
  }
  const std::string coefficient_path = argv[1];
  const int family = std::atoi(argv[2]);
  const int calibration_numerator = std::atoi(argv[3]);
  const int calibration_denominator = std::atoi(argv[4]);
  const int numerator_bound = std::atoi(argv[5]);
  const int denominator_bound = std::atoi(argv[6]);
  const int keep = std::atoi(argv[7]);
  if (family < 0 || family > 100 || calibration_denominator < 1
      || numerator_bound < 1 || numerator_bound > 100000
      || denominator_bound < 1 || denominator_bound > 5000
      || keep < 1 || keep > 20000) {
    std::cerr << "invalid family/calibration/bounds\n";
    return 2;
  }
  const Coefficients coefficients = read_coefficients(coefficient_path);
  const std::vector<PrimeTable> discovery =
      build_tables(coefficients, DISCOVERY_PRIMES);
  const std::vector<PrimeTable> held = build_tables(coefficients, HELD_PRIMES);
  std::priority_queue<Candidate, std::vector<Candidate>, BetterComparator> heap;
  std::uint64_t primitive_population = 0;
  std::uint64_t panel_excluded = 0;
  std::uint64_t evaluated_population = 0;

  for (int denominator = 1; denominator <= denominator_bound; ++denominator) {
    std::vector<int> multiplier;
    multiplier.reserve(discovery.size());
    for (const PrimeTable& table : discovery) {
      multiplier.push_back(
          denominator % table.prime == 0
              ? -1
              : power_mod(
                    denominator % table.prime, table.prime - 2, table.prime));
    }
    std::vector<int> residues(discovery.size(), 0);
    for (int numerator = 1; numerator <= numerator_bound; ++numerator) {
      for (std::size_t index = 0; index < discovery.size(); ++index) {
        if (multiplier[index] >= 0) {
          residues[index] += multiplier[index];
          if (residues[index] >= discovery[index].prime)
            residues[index] -= discovery[index].prime;
        }
      }
      if (std::gcd(numerator, denominator) != 1) continue;
      ++primitive_population;
      // The frozen max-root-200 panel already screened T=1,...,8.
      if (denominator == 1 && numerator <= 8) {
        ++panel_excluded;
        continue;
      }
      ++evaluated_population;
      Candidate candidate;
      candidate.numerator = numerator;
      candidate.denominator = denominator;
      for (std::size_t index = 0; index < discovery.size(); ++index) {
        const PrimeTable& table = discovery[index];
        const int symbol_index = multiplier[index] < 0
            ? table.prime : residues[index];
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
  for (Candidate& candidate : retained)
    add_band_score(candidate, held, true);

  Candidate calibration;
  calibration.numerator = calibration_numerator;
  calibration.denominator = calibration_denominator;
  add_band_score(calibration, discovery, false);
  add_band_score(calibration, held, true);

  std::cout << "MESTRE_RANK13_MULTIFAMILY_SCAN_V1\nF " << family << '\n' << "D";
  for (const int prime : DISCOVERY_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nH";
  for (const int prime : HELD_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nL " << table_digest(discovery) << ' '
            << table_digest(held)
            << "\nK " << calibration.numerator << ' '
            << calibration.denominator << ' '
            << score_text(calibration.discovery_score) << ' '
            << score_text(calibration.held_score) << ' '
            << calibration.discovery_good << ' ' << calibration.held_good << '\n';
  for (const Candidate& candidate : retained) {
    std::cout << "C " << candidate.numerator << ' ' << candidate.denominator << ' '
              << score_text(candidate.discovery_score) << ' '
              << score_text(candidate.held_score) << ' '
              << candidate.discovery_good << ' ' << candidate.held_good << '\n';
  }
  std::cout << "S " << numerator_bound << ' ' << denominator_bound << ' '
            << keep << ' ' << primitive_population << ' ' << panel_excluded << ' '
            << evaluated_population << ' ' << retained.size() << '\n';
  return 0;
}
