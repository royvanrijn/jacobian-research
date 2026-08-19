// Family-independent rational Nagao scanner.
//
// Input table format is produced by make_rational_nagao_tables.py.  For each
// prime p it stores p+1 projective local symbols: residues 0..p-1 and infinity.
// The hot loop therefore does no point counting and no large-integer arithmetic.
//
// Usage:
//   scan_rational_nagao_tables TABLE NUM DEN KEEP [SHARD_ID SHARDS]
//
// The scan visits positive reduced fractions a/b with
//   1 <= a <= NUM, 1 <= b <= DEN, gcd(a,b)=1.
// Optional sharding assigns denominator b to shard (b-1) mod SHARDS.
// Discovery scores choose the heap; held scores are evaluated only for survivors.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
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
  std::vector<LocalSymbol> symbols;  // p finite residues, then infinity
};

struct Candidate {
  int numerator = 0;
  int denominator = 1;
  std::int64_t discovery_score = 0;
  std::int64_t held_score = 0;
  int discovery_good = 0;
  int held_good = 0;
};

static int multiply_mod(std::int64_t a, std::int64_t b, int p) {
  return static_cast<int>((a * b) % p);
}

static int power_mod(int base, int exponent, int p) {
  int result = 1;
  while (exponent) {
    if (exponent & 1) result = multiply_mod(result, base, p);
    base = multiply_mod(base, base, p);
    exponent >>= 1;
  }
  return result;
}

static std::string score_text(std::int64_t units) {
  const bool negative = units < 0;
  const std::uint64_t absolute = negative
      ? static_cast<std::uint64_t>(-units)
      : static_cast<std::uint64_t>(units);
  const std::uint64_t whole = absolute / 1000000000000ULL;
  const std::uint64_t fraction = absolute % 1000000000000ULL;
  std::ostringstream out;
  if (negative) out << '-';
  out << whole << '.' << std::setw(12) << std::setfill('0') << fraction;
  return out.str();
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
    // priority_queue::top is the worst retained candidate.
    return better(left, right);
  }
};

struct Tables {
  std::string family;
  std::vector<PrimeTable> discovery;
  std::vector<PrimeTable> held;
};

static PrimeTable read_prime(std::istream& in) {
  std::string marker;
  int p = 0;
  if (!(in >> marker >> p) || marker != "P" || p < 2)
    throw std::runtime_error("malformed P record");
  PrimeTable table;
  table.prime = p;
  table.symbols.reserve(static_cast<std::size_t>(p + 1));
  for (int i = 0; i <= p; ++i) {
    int good = 0;
    int trace = 0;
    std::int64_t score = 0;
    if (!(in >> good >> trace >> score))
      throw std::runtime_error("truncated local-symbol table");
    table.symbols.push_back({good != 0, trace, score});
  }
  return table;
}

static Tables read_tables(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open table file: " + path);
  std::string header;
  std::getline(in, header);
  if (header != "RATIONAL_NAGAO_LOCAL_TABLE_V1")
    throw std::runtime_error("unexpected table header");

  Tables result;
  std::string fmark;
  int adeg = 0, bdeg = 0;
  if (!(in >> fmark >> result.family >> adeg >> bdeg) || fmark != "F")
    throw std::runtime_error("malformed family record");

  for (int band_index = 0; band_index < 2; ++band_index) {
    std::string bmark, label;
    int count = 0;
    if (!(in >> bmark >> label >> count) || bmark != "B" || count < 1)
      throw std::runtime_error("malformed band record");
    std::vector<PrimeTable>* target = nullptr;
    if (label == "D") target = &result.discovery;
    if (label == "H") target = &result.held;
    if (!target) throw std::runtime_error("unknown band label");
    target->reserve(count);
    for (int i = 0; i < count; ++i) target->push_back(read_prime(in));
  }

  std::string end;
  if (!(in >> end) || end != "END")
    throw std::runtime_error("missing END record");
  return result;
}

static void add_band_score(
    Candidate& candidate, const std::vector<PrimeTable>& tables, bool held) {
  std::int64_t score = 0;
  int good = 0;
  for (const PrimeTable& table : tables) {
    const int p = table.prime;
    int index = p;  // infinity
    if (candidate.denominator % p != 0) {
      const int inverse = power_mod(candidate.denominator % p, p - 2, p);
      index = multiply_mod(candidate.numerator % p, inverse, p);
    }
    const LocalSymbol& symbol = table.symbols[index];
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

int main(int argc, char** argv) {
  if (argc != 5 && argc != 7) {
    std::cerr << "usage: scan_rational_nagao_tables TABLE NUM DEN KEEP [SHARD_ID SHARDS]\n";
    return 2;
  }

  try {
    const std::string table_path = argv[1];
    const int numerator_bound = std::atoi(argv[2]);
    const int denominator_bound = std::atoi(argv[3]);
    const int keep = std::atoi(argv[4]);
    const int shard_id = argc == 7 ? std::atoi(argv[5]) : 0;
    const int shards = argc == 7 ? std::atoi(argv[6]) : 1;
    if (numerator_bound < 1 || denominator_bound < 1 || keep < 1) {
      std::cerr << "NUM, DEN and KEEP must be positive\n";
      return 2;
    }
    if (shards < 1 || shard_id < 0 || shard_id >= shards) {
      std::cerr << "require SHARDS>=1 and 0<=SHARD_ID<SHARDS\n";
      return 2;
    }

    const Tables tables = read_tables(table_path);
    std::priority_queue<Candidate, std::vector<Candidate>, BetterComparator> heap;
    std::uint64_t primitive_population = 0;
    std::uint64_t evaluated_population = 0;

    // Denominator-major order lets us update a/b mod p by repeated addition of
    // b^{-1}. If p|b every candidate on that denominator maps to infinity.
    // Sharding by denominator keeps shards disjoint and deterministic.
    for (int denominator = 1; denominator <= denominator_bound; ++denominator) {
      if ((denominator - 1) % shards != shard_id) continue;

      const std::size_t nprimes = tables.discovery.size();
      std::vector<int> step(nprimes, -1);
      std::vector<int> residue(nprimes, 0);
      for (std::size_t i = 0; i < nprimes; ++i) {
        const int p = tables.discovery[i].prime;
        if (denominator % p != 0)
          step[i] = power_mod(denominator % p, p - 2, p);
      }

      for (int numerator = 1; numerator <= numerator_bound; ++numerator) {
        for (std::size_t i = 0; i < nprimes; ++i) {
          if (step[i] >= 0) {
            residue[i] += step[i];
            if (residue[i] >= tables.discovery[i].prime)
              residue[i] -= tables.discovery[i].prime;
          }
        }

        if (std::gcd(numerator, denominator) != 1) continue;
        ++primitive_population;
        ++evaluated_population;

        Candidate candidate;
        candidate.numerator = numerator;
        candidate.denominator = denominator;
        for (std::size_t i = 0; i < nprimes; ++i) {
          const PrimeTable& table = tables.discovery[i];
          const int index = step[i] < 0 ? table.prime : residue[i];
          const LocalSymbol& symbol = table.symbols[index];
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
      add_band_score(candidate, tables.held, true);

    std::cout << "RATIONAL_NAGAO_SCAN_V1\n";
    std::cout << "F " << tables.family << '\n';
    std::cout << "D";
    for (const auto& table : tables.discovery) std::cout << ' ' << table.prime;
    std::cout << "\nH";
    for (const auto& table : tables.held) std::cout << ' ' << table.prime;
    std::cout << "\nR " << shard_id << ' ' << shards << '\n';
    for (const Candidate& candidate : retained) {
      std::cout << "C " << candidate.numerator << ' ' << candidate.denominator << ' '
                << score_text(candidate.discovery_score) << ' '
                << score_text(candidate.held_score) << ' '
                << candidate.discovery_good << ' ' << candidate.held_good << '\n';
    }
    std::cout << "S " << numerator_bound << ' ' << denominator_bound << ' '
              << keep << ' ' << primitive_population << ' '
              << evaluated_population << ' ' << retained.size() << '\n';
    return 0;
  } catch (const std::exception& exc) {
    std::cerr << "error: " << exc.what() << '\n';
    return 1;
  }
}
