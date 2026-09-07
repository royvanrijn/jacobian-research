// Deterministic staged rank-jump sweep for the Fermigier adapter family.
//
// Build with:
//   g++ -O3 -march=native -fopenmp -std=c++20 \
//     -o /tmp/fermigier-score-sweep \
//     elliptic-curves/ecsearch/fermigier_score_sweep.cpp
//
// At every scoring prime this computes the exact mean trace on the good
// projective fibres of the family, and scores a specialization by its trace
// residual from that mean.  Thus the cheap signal is aimed at the exceptional
// quotient over the generic rank-12 subgroup, rather than at total rank.
//
// The output is a heuristic ranking, never a rank or conductor certificate.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <omp.h>

namespace {

struct TraceTable {
  int prime;
  std::vector<std::int16_t> trace;
  std::vector<std::uint8_t> good;
  std::vector<std::uint16_t> inverse;
  double logarithm;
  int good_trace_sum;
  int good_fibre_count;
};

struct Candidate {
  int numerator;
  int denominator;
  float score;
};

struct ScoreFeatures {
  double s0 = 0;
  double s5 = 0;
  int primes_used = 0;
  int bad_specialization_primes = 0;
};

long long power_mod(long long base, long long exponent, int modulus) {
  long long result = 1;
  while (exponent != 0) {
    if ((exponent & 1) != 0) {
      result = result * base % modulus;
    }
    base = base * base % modulus;
    exponent >>= 1;
  }
  return result;
}

int decimal_mod(const std::string& value, int modulus) {
  long long residue = 0;
  int sign = 1;
  for (char character : value) {
    if (character == '-') {
      sign = -1;
    } else {
      residue = (10 * residue + character - '0') % modulus;
    }
  }
  return static_cast<int>((sign * residue % modulus + modulus) % modulus);
}

std::vector<int> primes_through(int limit) {
  std::vector<bool> is_prime(limit + 1, true);
  is_prime[0] = false;
  is_prime[1] = false;
  for (int prime = 2; prime * prime <= limit; ++prime) {
    if (!is_prime[prime]) {
      continue;
    }
    for (int multiple = prime * prime; multiple <= limit; multiple += prime) {
      is_prime[multiple] = false;
    }
  }
  std::vector<int> result;
  for (int value = 2; value <= limit; ++value) {
    if (is_prime[value]) {
      result.push_back(value);
    }
  }
  return result;
}

std::vector<TraceTable> build_trace_tables() {
  const std::vector<std::string> a2 = {
      "298803565660", "1718550", "-8"};
  const std::vector<std::string> a4 = {
      "27856983036916830925012", "317946481466562025",
      "1028009011008", "-18151200", "-64"};
  const std::vector<std::string> a6 = {
      "829998138277457737118423455411406",
      "14321756340366264921294086000",
      "109667821527431621677482", "573566223236700400",
      "-31029122010320", "35222400", "512"};

  std::vector<TraceTable> tables;
  for (int prime : primes_through(2000)) {
    // The family model has no good fibre at its excluded primes 2 and 3.
    // Ordinary good-reduction trace residuals are undefined there.
    if (prime == 2 || prime == 3) {
      continue;
    }
    std::vector<int> character(prime);
    for (int value = 1; value < prime; ++value) {
      character[value] =
          power_mod(value, (prime - 1) / 2, prime) == 1 ? 1 : -1;
    }
    std::vector<int> c2, c4, c6;
    for (const auto& value : a2) c2.push_back(decimal_mod(value, prime));
    for (const auto& value : a4) c4.push_back(decimal_mod(value, prime));
    for (const auto& value : a6) c6.push_back(decimal_mod(value, prime));

    TraceTable table{prime,
                     std::vector<std::int16_t>(prime + 1),
                     std::vector<std::uint8_t>(prime + 1),
                     std::vector<std::uint16_t>(prime),
                     std::log(static_cast<double>(prime)),
                     0,
                     0};
    for (int value = 1; value < prime; ++value) {
      table.inverse[value] =
          static_cast<std::uint16_t>(power_mod(value, prime - 2, prime));
    }

#pragma omp parallel for schedule(dynamic)
    for (int parameter = 0; parameter <= prime; ++parameter) {
      long long a1_coefficient;
      long long a2_coefficient;
      long long a3_coefficient;
      long long a4_coefficient;
      long long a6_coefficient;
      if (parameter < prime) {
        const long long square =
            static_cast<long long>(parameter) * parameter % prime;
        long long power = 1;
        a2_coefficient = 0;
        for (int coefficient : c2) {
          a2_coefficient = (a2_coefficient + coefficient * power) % prime;
          power = power * square % prime;
        }
        power = 1;
        a4_coefficient = 0;
        for (int coefficient : c4) {
          a4_coefficient = (a4_coefficient + coefficient * power) % prime;
          power = power * square % prime;
        }
        power = 1;
        a6_coefficient = 0;
        for (int coefficient : c6) {
          a6_coefficient = (a6_coefficient + coefficient * power) % prime;
          power = power * square % prime;
        }
        a1_coefficient = 1;
        a3_coefficient = 1;
      } else {
        // Projective parameter [a:b]=[1:0] in the integral homogeneous model.
        a1_coefficient = 0;
        a2_coefficient = (prime - 8 % prime) % prime;
        a3_coefficient = 0;
        a4_coefficient = (prime - 64 % prime) % prime;
        a6_coefficient = 512 % prime;
      }

      int character_sum = 0;
      for (int x_coordinate = 0; x_coordinate < prime; ++x_coordinate) {
        const long long linear_y =
            (a1_coefficient * x_coordinate + a3_coefficient) % prime;
        const long long rhs =
            (((static_cast<long long>(x_coordinate) * x_coordinate % prime *
                   x_coordinate +
               a2_coefficient * static_cast<long long>(x_coordinate) % prime *
                   x_coordinate) %
                  prime +
              a4_coefficient * x_coordinate) %
                 prime +
             a6_coefficient) %
            prime;
        character_sum += character[(linear_y * linear_y + 4 * rhs) % prime];
      }
      table.trace[parameter] = static_cast<std::int16_t>(-character_sum);

      // A specialization can be scored by an ordinary a_p only at good
      // reduction.  The family baseline U(F_p) is exactly this good-fibre
      // locus, including the projective infinity fibre when it is smooth.
      const long long b2 =
          (a1_coefficient * a1_coefficient + 4 * a2_coefficient) % prime;
      const long long b4 =
          (a1_coefficient * a3_coefficient + 2 * a4_coefficient) % prime;
      const long long b6 =
          (a3_coefficient * a3_coefficient + 4 * a6_coefficient) % prime;
      const long long b8 =
          (a1_coefficient * a1_coefficient * a6_coefficient +
           4 * a2_coefficient * a6_coefficient -
           a1_coefficient * a3_coefficient * a4_coefficient +
           a2_coefficient * a3_coefficient * a3_coefficient -
           a4_coefficient * a4_coefficient) %
          prime;
      const long long discriminant =
          (-b2 * b2 * b8 - 8 * b4 * b4 * b4 - 27 * b6 * b6 +
           9 * b2 * b4 * b6) %
          prime;
      table.good[parameter] = discriminant != 0;
    }
    for (int parameter = 0; parameter <= prime; ++parameter) {
      if (table.good[parameter] == 0) {
        continue;
      }
      table.good_trace_sum += table.trace[parameter];
      ++table.good_fibre_count;
    }
    if (table.good_fibre_count == 0) {
      // Fixed bad primes of this integral family model carry no ordinary
      // good-fibre trace signal.  Exclude them exactly as a specialization-
      // bad prime is excluded below; aborting would make the declared sweep
      // impossible before reaching its first candidate.
      std::cerr << "skipping family-bad prime=" << prime << '\n';
      continue;
    }
    tables.push_back(std::move(table));
  }
  return tables;
}

int table_count_through(const std::vector<TraceTable>& tables, int bound) {
  int count = 0;
  while (count < static_cast<int>(tables.size()) &&
         tables[count].prime <= bound) {
    ++count;
  }
  return count;
}

ScoreFeatures score_features(const std::vector<TraceTable>& tables,
                             int begin,
                             int end,
                             int numerator,
                             int denominator) {
  ScoreFeatures result;
  for (int index = begin; index < end; ++index) {
    const auto& table = tables[index];
    const int denominator_residue = denominator % table.prime;
    const int parameter = denominator_residue == 0
                              ? table.prime
                              : static_cast<int>(
                                    static_cast<long long>(numerator % table.prime) *
                                    table.inverse[denominator_residue] % table.prime);
    if (table.good[parameter] == 0) {
      ++result.bad_specialization_primes;
      continue;
    }
    const int trace = table.trace[parameter];
    const double family_mean = static_cast<double>(table.good_trace_sum) /
                               static_cast<double>(table.good_fibre_count);
    // S0 is the linear trace residual used in Mestre--Nagao scoring.  S5 is
    // the corresponding local Euler-factor residual.  Both vanish after
    // averaging their first-order trace signal over the family.
    result.s0 += (family_mean - trace) * table.logarithm / table.prime;
    result.s5 += std::log((table.prime + 1.0 - trace) /
                          (table.prime + 1.0 - family_mean));
    ++result.primes_used;
  }
  return result;
}

double rank_jump_score(const ScoreFeatures& cumulative,
                       const ScoreFeatures& window,
                       int bound,
                       int previous_bound) {
  // This fixed, declared combination is only a staging rule.  It deliberately
  // keeps a disjoint window term: residual signals need not improve
  // monotonically with a larger prime cutoff.
  const double cumulative_s0 = cumulative.s0 / std::log(bound);
  if (previous_bound == 0) {
    return cumulative_s0 + cumulative.s5;
  }
  const double window_scale =
      std::log(static_cast<double>(bound) / previous_bound);
  const double window_s0 = window.s0 / window_scale;
  return cumulative_s0 + cumulative.s5 + 0.5 * (window_s0 + window.s5);
}

void retain_top(std::vector<Candidate>& candidates, std::size_t count) {
  if (candidates.size() <= count) {
    return;
  }
  auto comparator = [](const Candidate& left, const Candidate& right) {
    return left.score > right.score;
  };
  std::nth_element(candidates.begin(), candidates.begin() + count,
                   candidates.end(), comparator);
  candidates.resize(count);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "usage: " << argv[0]
              << " MAX_NUMERATOR MAX_DENOMINATOR [OUTPUT_COUNT]"
                 " [MIN_NUMERATOR]\n";
    return 2;
  }
  const int maximum_numerator = std::atoi(argv[1]);
  const int maximum_denominator = std::atoi(argv[2]);
  const int output_count = argc >= 4 ? std::atoi(argv[3]) : 20000;
  const int minimum_numerator = argc == 5 ? std::atoi(argv[4]) : 1;
  if (maximum_numerator <= 0 || maximum_denominator <= 0 || output_count <= 0 ||
      minimum_numerator <= 0 || minimum_numerator > maximum_numerator) {
    std::cerr << "all bounds must be positive\n";
    return 2;
  }

  const auto start = std::chrono::steady_clock::now();
  const auto tables = build_trace_tables();
  const std::vector<int> bounds = {100, 200, 400, 1000, 2000};
  std::vector<int> table_counts;
  std::vector<double> e22_scores;
  std::vector<double> rank20_scores;
  for (std::size_t stage = 0; stage < bounds.size(); ++stage) {
    const int bound = bounds[stage];
    const int count = table_count_through(tables, bound);
    table_counts.push_back(count);
    const int previous_count = stage == 0 ? 0 : table_counts[stage - 1];
    const ScoreFeatures cumulative =
        score_features(tables, 0, count, 19754, 39);
    const ScoreFeatures window = stage == 0
                                     ? cumulative
                                     : score_features(tables, previous_count,
                                                      count, 19754, 39);
    e22_scores.push_back(
        rank_jump_score(cumulative, window, bound,
                        stage == 0 ? 0 : bounds[stage - 1]));
    const ScoreFeatures rank20_cumulative =
        score_features(tables, 0, count, 28917, 20);
    const ScoreFeatures rank20_window =
        stage == 0
            ? rank20_cumulative
            : score_features(tables, previous_count, count, 28917, 20);
    rank20_scores.push_back(
        rank_jump_score(rank20_cumulative, rank20_window, bound,
                        stage == 0 ? 0 : bounds[stage - 1]));
  }

  const int thread_count = omp_get_max_threads();
  std::vector<std::vector<Candidate>> thread_candidates(thread_count);
#pragma omp parallel
  {
    auto& local = thread_candidates[omp_get_thread_num()];
    local.reserve(static_cast<std::size_t>(maximum_numerator - minimum_numerator + 1) *
                  maximum_denominator / thread_count * 7 / 10);
#pragma omp for schedule(dynamic, 1)
    for (int denominator = 1; denominator <= maximum_denominator;
         ++denominator) {
      for (int numerator = minimum_numerator; numerator <= maximum_numerator;
           ++numerator) {
        if (std::gcd(numerator, denominator) != 1) {
          continue;
        }
        const ScoreFeatures initial = score_features(
            tables, 0, table_counts[0], numerator, denominator);
        local.push_back(Candidate{numerator, denominator,
                                  static_cast<float>(rank_jump_score(
                                      initial, initial, bounds[0], 0))});
      }
    }
  }

  std::size_t total = 0;
  for (const auto& local : thread_candidates) total += local.size();
  std::vector<Candidate> candidates;
  candidates.reserve(total);
  for (auto& local : thread_candidates) {
    candidates.insert(candidates.end(), local.begin(), local.end());
    std::vector<Candidate>().swap(local);
  }

  // Wide first-stage retention is intentional: the known E22 has only a
  // moderately exceptional score at the smallest prime bounds.  Scale the
  // later caps with an explicitly larger requested output so held-out-prime
  // reranking can inspect a deeper p<=2000 tail without silently receiving
  // only the historical 100,000 survivors.
  const std::size_t requested = static_cast<std::size_t>(output_count);
  const std::vector<std::size_t> caps = {
      std::max<std::size_t>(8000000, requested),
      std::max<std::size_t>(2000000, requested),
      std::max<std::size_t>(500000, requested),
      std::max<std::size_t>(100000, requested),
      std::max<std::size_t>(20000, requested)};
  for (std::size_t stage = 0; stage < bounds.size(); ++stage) {
    if (stage != 0) {
#pragma omp parallel for schedule(static)
      for (std::size_t index = 0; index < candidates.size(); ++index) {
        const ScoreFeatures cumulative = score_features(
            tables, 0, table_counts[stage], candidates[index].numerator,
            candidates[index].denominator);
        const ScoreFeatures window = score_features(
            tables, table_counts[stage - 1], table_counts[stage],
            candidates[index].numerator, candidates[index].denominator);
        candidates[index].score = static_cast<float>(rank_jump_score(
            cumulative, window, bounds[stage], bounds[stage - 1]));
      }
    }
    const std::size_t above_e22 = static_cast<std::size_t>(std::count_if(
        candidates.begin(), candidates.end(), [&](const Candidate& candidate) {
          return candidate.score >= e22_scores[stage] - 1e-6;
        }));
    const std::size_t above_rank20 = static_cast<std::size_t>(std::count_if(
        candidates.begin(), candidates.end(), [&](const Candidate& candidate) {
          return candidate.score >= rank20_scores[stage] - 1e-6;
        }));
    std::cerr << "stage=" << bounds[stage]
              << " candidates=" << candidates.size()
              << " e22_rank_jump_score=" << e22_scores[stage]
              << " at_or_above_e22=" << above_e22
              << " rank20_rank_jump_score=" << rank20_scores[stage]
              << " at_or_above_rank20=" << above_rank20 << '\n';
    retain_top(candidates, caps[stage]);
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& left, const Candidate& right) {
              if (left.score != right.score) return left.score > right.score;
              if (left.numerator != right.numerator) {
                return left.numerator < right.numerator;
              }
              return left.denominator < right.denominator;
            });
  if (candidates.size() > static_cast<std::size_t>(output_count)) {
    candidates.resize(output_count);
  }
  const int final_count = table_counts.back();
  const int final_window_begin = table_counts[table_counts.size() - 2];
  std::cout << "# rank\tnumerator\tdenominator\trank_jump_score"
               "\tcumulative_S0\tcumulative_S5\twindow_S0\twindow_S5"
               "\tgood_primes\tbad_specialization_primes\n";
  std::cout.precision(10);
  for (std::size_t index = 0; index < candidates.size(); ++index) {
    const auto& candidate = candidates[index];
    const ScoreFeatures cumulative = score_features(
        tables, 0, final_count, candidate.numerator, candidate.denominator);
    const ScoreFeatures window = score_features(
        tables, final_window_begin, final_count, candidate.numerator,
        candidate.denominator);
    std::cout << index + 1 << '\t' << candidate.numerator << '\t'
              << candidate.denominator << '\t' << candidate.score << '\t'
              << cumulative.s0 << '\t' << cumulative.s5 << '\t'
              << window.s0 << '\t' << window.s5 << '\t'
              << cumulative.primes_used << '\t'
              << cumulative.bad_specialization_primes << '\n';
  }
  const double seconds = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count();
  std::cerr << "elapsed_seconds=" << seconds << '\n';
}
