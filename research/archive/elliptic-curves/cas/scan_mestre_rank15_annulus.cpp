// Exact-local scanner for a disjoint annulus around the Mestre T=490/9 lead.
//
// Population: primitive (a,b), 4097 <= b <= 16000,
//   a = nearest(490*b/9) + delta, 5 <= |delta| <= 16.
// The two peer Farey rays are excluded exactly after rational reduction.  To
// prevent overlap with the peer CRT grids before their numerator manifests
// are known, every denominator dividing any raw grid denominator 552+37*j or
// 553+41*j (0 <= j <= 255) is conservatively excluded.

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
  int offset = 0;
  std::int64_t discovery_score = 0;
  std::int64_t held_score = 0;
  int discovery_good = 0;
  int held_good = 0;
};

static constexpr int DENOMINATOR_MIN = 4097;
static constexpr int DENOMINATOR_MAX = 16000;
static constexpr std::array<int, 14> DISCOVERY_PRIMES{
    587, 593, 599, 601, 607, 613, 617,
    619, 631, 641, 643, 647, 653, 659};
static constexpr std::array<int, 14> HELD_PRIMES{
    661, 673, 677, 683, 691, 701, 709,
    719, 727, 733, 739, 743, 751, 757};

// Primitive Jacobian of roots (0,7,121,128,183,194).
static const std::array<std::string, 9> A_COEFFICIENTS{{
    "-23176885126717023443712", "0", "-30716437929772416480", "0",
    "5349668724315261", "0", "-44248113000", "0", "-29428272"}};
static const std::array<std::string, 13> B_COEFFICIENTS{{
    "867772568274887992065902700355584", "0",
    "4115037444689925520425570654720", "0",
    "-945618543766059469223433696", "0", "142515300626499263347710", "0",
    "-16911444014164334952", "0", "138585089916000", "0", "61446231936"}};

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

static int trace_of_frobenius(int coefficient_a, int coefficient_b, int prime) {
  int character_sum = 0;
  for (int x = 0; x < prime; ++x) {
    const int rhs = (multiply_mod(multiply_mod(x, x, prime), x, prime)
                     + multiply_mod(coefficient_a, x, prime) + coefficient_b)
                    % prime;
    if (rhs == 0) continue;
    character_sum += power_mod(rhs, (prime - 1) / 2, prime) == 1 ? 1 : -1;
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

static std::uint64_t pair_key(int numerator, int denominator) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(numerator)) << 32)
         | static_cast<std::uint32_t>(denominator);
}

static std::unordered_set<std::uint64_t> peer_farey_manifest() {
  std::unordered_set<std::uint64_t> answer;
  for (int m = 61; m <= 2048; ++m) {
    for (const auto& offsets : std::array<std::array<int, 2>, 2>{{
             {{381, 7}}, {{109, 2}}}}) {
      int numerator = 490 * m + offsets[0];
      int denominator = 9 * m + offsets[1];
      const int divisor = std::gcd(numerator, denominator);
      numerator /= divisor;
      denominator /= divisor;
      answer.insert(pair_key(numerator, denominator));
    }
  }
  return answer;
}

static std::vector<bool> excluded_grid_divisors() {
  std::vector<bool> answer(DENOMINATOR_MAX + 1, false);
  for (int j = 0; j <= 255; ++j) {
    for (const int raw : {552 + 37 * j, 553 + 41 * j}) {
      for (int divisor = 1; divisor * divisor <= raw; ++divisor) {
        if (raw % divisor) continue;
        if (divisor <= DENOMINATOR_MAX) answer[divisor] = true;
        const int other = raw / divisor;
        if (other <= DENOMINATOR_MAX) answer[other] = true;
      }
    }
  }
  return answer;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: scan_mestre_rank15_annulus KEEP\n";
    return 2;
  }
  const int keep = std::atoi(argv[1]);
  if (keep < 1 || keep > 20000) {
    std::cerr << "KEEP must lie in [1,20000]\n";
    return 2;
  }
  const auto discovery = build_tables(DISCOVERY_PRIMES);
  const auto held = build_tables(HELD_PRIMES);
  const auto farey = peer_farey_manifest();
  const auto grid_divisors = excluded_grid_divisors();
  std::priority_queue<Candidate, std::vector<Candidate>, BetterComparator> heap;
  std::uint64_t raw_population = 0;
  std::uint64_t nonprimitive = 0;
  std::uint64_t grid_excluded = 0;
  std::uint64_t farey_excluded = 0;
  std::uint64_t evaluated = 0;

  for (int denominator = DENOMINATOR_MIN; denominator <= DENOMINATOR_MAX;
       ++denominator) {
    const int nearest = (490 * denominator + 4) / 9;
    for (int offset = -16; offset <= 16; ++offset) {
      if (offset > -5 && offset < 5) continue;
      ++raw_population;
      const int numerator = nearest + offset;
      if (std::gcd(numerator, denominator) != 1) {
        ++nonprimitive;
        continue;
      }
      if (grid_divisors[denominator]) {
        ++grid_excluded;
        continue;
      }
      if (farey.count(pair_key(numerator, denominator))) {
        ++farey_excluded;
        continue;
      }
      ++evaluated;
      Candidate candidate;
      candidate.numerator = numerator;
      candidate.denominator = denominator;
      candidate.offset = offset;
      add_band_score(candidate, discovery, false);
      if (static_cast<int>(heap.size()) < keep) {
        heap.push(candidate);
      } else if (better(candidate, heap.top())) {
        heap.pop();
        heap.push(candidate);
      }
    }
  }

  std::vector<Candidate> retained;
  while (!heap.empty()) {
    retained.push_back(heap.top());
    heap.pop();
  }
  std::sort(retained.begin(), retained.end(), better);
  for (Candidate& candidate : retained) add_band_score(candidate, held, true);

  std::cout << "MESTRE_RANK15_ANNULUS_SCAN_V1\nD";
  for (const int prime : DISCOVERY_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nH";
  for (const int prime : HELD_PRIMES) std::cout << ' ' << prime;
  std::cout << "\nL " << table_digest(discovery) << ' ' << table_digest(held)
            << '\n';
  for (const Candidate& candidate : retained) {
    std::cout << "C " << candidate.numerator << ' ' << candidate.denominator << ' '
              << candidate.offset << ' ' << score_text(candidate.discovery_score)
              << ' ' << score_text(candidate.held_score) << ' '
              << candidate.discovery_good << ' ' << candidate.held_good << '\n';
  }
  std::cout << "S " << DENOMINATOR_MIN << ' ' << DENOMINATOR_MAX << ' ' << keep
            << ' ' << raw_population << ' ' << nonprimitive << ' '
            << grid_excluded << ' ' << farey_excluded << ' ' << evaluated << ' '
            << retained.size() << '\n';
  return 0;
}
