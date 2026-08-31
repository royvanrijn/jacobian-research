// Fast table-lookup scan for the exact q12/orbit5867 rootless R17 family.
//
// The Python companion exports every point of P^1(F_p).  This hot loop does no
// coefficient or section evaluation: it enumerates signed primitive (a:b),
// including zero and infinity, and exactly reproduces the staged bucketed
// Pareto retention used by search_h92_q12o5867_rootless_nagao.py.
//
// Usage:
//   scan TABLE NUM DEN BUCKET_WIDTH KEEP1,KEEP2,... FINALISTS OUTPUT.json
//        [PARAMETER_SCALE [CONTROL_A/B,...]]
//
// The optional scale searches the rational chart u=PARAMETER_SCALE*v while
// retaining v in balanced height buckets.  It must be nonzero and invertible
// modulo every table prime.  This is useful when the exact binary A_8,B_12
// model is badly conditioned in its inherited projective coordinate.  When a
// control list is supplied, the scanner scores the *complete* box using every
// block, ranks by the minimum centered/standardized block signal, and reports
// exact population ranks for the controls.  No staged pruning is used in that
// calibration mode.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kMaximumBlocks = 16;

struct LocalSymbol {
  bool good = false;
  bool singular = true;
  int point_count = -1;
  int trace = 0;
  std::int64_t score_units = 0;
};

struct PrimeTable {
  int prime = 0;
  std::vector<LocalSymbol> symbols;  // affine 0..p-1, then infinity at p
  double good_mean = 0.0;
  double good_standard_deviation = 0.0;
};

struct Tables {
  std::string model_sha256;
  int a_degree = 0;
  int b_degree = 0;
  std::int64_t score_scale = 0;
  std::vector<std::vector<PrimeTable>> blocks;
};

struct Candidate {
  std::int64_t numerator = 0;
  std::int64_t denominator = 1;
  std::int64_t height = 1;
  std::array<std::int64_t, kMaximumBlocks> block_scores{};
  std::array<double, kMaximumBlocks> standardized_block_scores{};
  int completed_blocks = 0;
  int good_primes = 0;
  int bad_primes = 0;

  std::int64_t total_score() const {
    std::int64_t total = 0;
    for (int i = 0; i < completed_blocks; ++i) total += block_scores[i];
    return total;
  }

  std::int64_t minimum_block_score() const {
    if (completed_blocks == 0) throw std::logic_error("unscored candidate");
    std::int64_t answer = block_scores[0];
    for (int i = 1; i < completed_blocks; ++i)
      answer = std::min(answer, block_scores[i]);
    return answer;
  }

  double minimum_standardized_block_score() const {
    if (completed_blocks == 0) throw std::logic_error("unscored candidate");
    double answer = standardized_block_scores[0];
    for (int i = 1; i < completed_blocks; ++i)
      answer = std::min(answer, standardized_block_scores[i]);
    return answer;
  }

  double mean_standardized_block_score() const {
    if (completed_blocks == 0) throw std::logic_error("unscored candidate");
    double answer = 0.0;
    for (int i = 0; i < completed_blocks; ++i)
      answer += standardized_block_scores[i];
    return answer / completed_blocks;
  }

  std::string parameter() const {
    if (denominator == 0) return "infinity";
    return std::to_string(numerator) + "/" + std::to_string(denominator);
  }
};

struct StageSummary {
  int stage = 0;
  std::vector<int> primes;
  int cap_per_bucket = 0;
  std::uint64_t population_scored = 0;
  std::uint64_t bucket_count = 0;
  std::uint64_t pareto_insertions = 0;
  std::uint64_t retained_count = 0;
  double seconds = 0.0;
};

int multiply_mod(std::int64_t left, std::int64_t right, int prime) {
  return static_cast<int>((left * right) % prime);
}

int power_mod(int base, int exponent, int prime) {
  int result = 1;
  while (exponent) {
    if (exponent & 1) result = multiply_mod(result, base, prime);
    base = multiply_mod(base, base, prime);
    exponent >>= 1;
  }
  return result;
}

int normalized_mod(std::int64_t value, int prime) {
  int answer = static_cast<int>(value % prime);
  return answer < 0 ? answer + prime : answer;
}

Tables read_tables(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open table file: " + path);
  std::string header;
  std::getline(input, header);
  if (header != "H92_Q12O5867_PROJECTIVE_NAGAO_TABLE_V1")
    throw std::runtime_error("unexpected table header");

  Tables tables;
  std::string marker;
  if (!(input >> marker >> tables.model_sha256 >> tables.a_degree >>
        tables.b_degree >> tables.score_scale) || marker != "M")
    throw std::runtime_error("malformed model record");
  int block_count = 0;
  if (!(input >> marker >> block_count) || marker != "C" || block_count < 1 ||
      block_count > kMaximumBlocks)
    throw std::runtime_error("invalid block count");
  tables.blocks.reserve(static_cast<std::size_t>(block_count));
  for (int block_index = 1; block_index <= block_count; ++block_index) {
    int stored_index = 0;
    int prime_count = 0;
    if (!(input >> marker >> stored_index >> prime_count) || marker != "B" ||
        stored_index != block_index || prime_count < 1)
      throw std::runtime_error("malformed block record");
    std::vector<PrimeTable> block;
    block.reserve(static_cast<std::size_t>(prime_count));
    for (int prime_index = 0; prime_index < prime_count; ++prime_index) {
      PrimeTable table;
      int symbol_count = 0;
      if (!(input >> marker >> table.prime >> symbol_count) || marker != "P" ||
          symbol_count != table.prime + 1)
        throw std::runtime_error("incomplete projective prime record");
      table.symbols.reserve(static_cast<std::size_t>(symbol_count));
      for (int symbol_index = 0; symbol_index < symbol_count; ++symbol_index) {
        int good = 0;
        int singular = 0;
        LocalSymbol symbol;
        if (!(input >> good >> singular >> symbol.point_count >> symbol.trace >>
              symbol.score_units))
          throw std::runtime_error("truncated local-symbol table");
        symbol.good = good != 0;
        symbol.singular = singular != 0;
        if (symbol.good == symbol.singular)
          throw std::runtime_error("inconsistent good/singular flags");
        table.symbols.push_back(symbol);
      }
      double sum = 0.0;
      int good_count = 0;
      for (const LocalSymbol& symbol : table.symbols) {
        if (!symbol.good) continue;
        sum += static_cast<double>(symbol.score_units);
        ++good_count;
      }
      if (good_count < 2)
        throw std::runtime_error("a prime table has fewer than two good fibres");
      table.good_mean = sum / good_count;
      double squared_deviation_sum = 0.0;
      for (const LocalSymbol& symbol : table.symbols) {
        if (!symbol.good) continue;
        const double deviation =
            static_cast<double>(symbol.score_units) - table.good_mean;
        squared_deviation_sum += deviation * deviation;
      }
      table.good_standard_deviation =
          std::sqrt(squared_deviation_sum / good_count);
      if (!(table.good_standard_deviation > 0.0))
        throw std::runtime_error("a prime table has zero Nagao variance");
      block.push_back(std::move(table));
    }
    tables.blocks.push_back(std::move(block));
  }
  if (!(input >> marker) || marker != "END")
    throw std::runtime_error("missing END record");
  return tables;
}

std::vector<int> parse_keeps(const std::string& text) {
  std::vector<int> result;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    const int value = std::stoi(token);
    if (value < 1) throw std::runtime_error("keep counts must be positive");
    result.push_back(value);
  }
  if (result.empty()) throw std::runtime_error("no keep counts supplied");
  return result;
}

bool dominates(const Candidate& left, const Candidate& right) {
  if (left.completed_blocks != right.completed_blocks)
    throw std::logic_error("incomparable sieve stages");
  const bool weak =
      left.total_score() >= right.total_score() &&
      left.good_primes >= right.good_primes &&
      left.bad_primes <= right.bad_primes && left.height <= right.height;
  const bool strict =
      left.total_score() != right.total_score() ||
      left.good_primes != right.good_primes ||
      left.bad_primes != right.bad_primes || left.height != right.height;
  return weak && strict;
}

bool better(const Candidate& left, const Candidate& right) {
  if (left.total_score() != right.total_score())
    return left.total_score() > right.total_score();
  if (left.minimum_block_score() != right.minimum_block_score())
    return left.minimum_block_score() > right.minimum_block_score();
  if (left.good_primes != right.good_primes)
    return left.good_primes > right.good_primes;
  if (left.bad_primes != right.bad_primes)
    return left.bad_primes < right.bad_primes;
  if (left.height != right.height) return left.height < right.height;
  if (left.denominator != right.denominator)
    return left.denominator < right.denominator;
  return left.numerator < right.numerator;
}

bool calibration_better(const Candidate& left, const Candidate& right) {
  if (left.minimum_standardized_block_score() !=
      right.minimum_standardized_block_score())
    return left.minimum_standardized_block_score() >
           right.minimum_standardized_block_score();
  if (left.mean_standardized_block_score() !=
      right.mean_standardized_block_score())
    return left.mean_standardized_block_score() >
           right.mean_standardized_block_score();
  if (left.good_primes != right.good_primes)
    return left.good_primes > right.good_primes;
  if (left.bad_primes != right.bad_primes)
    return left.bad_primes < right.bad_primes;
  if (left.height != right.height) return left.height < right.height;
  if (left.denominator != right.denominator)
    return left.denominator < right.denominator;
  return left.numerator < right.numerator;
}

std::vector<std::pair<std::int64_t, std::int64_t>> parse_controls(
    const std::string& text) {
  std::vector<std::pair<std::int64_t, std::int64_t>> result;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    const std::size_t slash = token.find('/');
    if (slash == std::string::npos)
      throw std::runtime_error("a control parameter is not a/b");
    std::int64_t numerator = std::stoll(token.substr(0, slash));
    std::int64_t denominator = std::stoll(token.substr(slash + 1));
    if (denominator <= 0 || std::gcd(numerator < 0 ? -numerator : numerator,
                                     denominator) != 1)
      throw std::runtime_error("control parameters must be primitive with b>0");
    result.emplace_back(numerator, denominator);
  }
  if (result.empty()) throw std::runtime_error("no positive controls supplied");
  return result;
}

void add_calibration_symbol(Candidate& candidate, const PrimeTable& table,
                            const LocalSymbol& symbol, int block_index) {
  if (symbol.good) {
    candidate.block_scores[block_index] += symbol.score_units;
    candidate.standardized_block_scores[block_index] +=
        (static_cast<double>(symbol.score_units) - table.good_mean) /
        table.good_standard_deviation;
    ++candidate.good_primes;
  } else {
    // Mean imputation gives standardized contribution zero.  The bad count is
    // retained as a visible heuristic choice, never as a rank conclusion.
    ++candidate.bad_primes;
  }
}

Candidate score_calibration_pair(std::int64_t numerator,
                                 std::int64_t denominator,
                                 const Tables& tables) {
  Candidate candidate;
  candidate.numerator = numerator;
  candidate.denominator = denominator;
  candidate.height = denominator == 0
                         ? 1
                         : std::max(numerator < 0 ? -numerator : numerator,
                                    denominator);
  candidate.completed_blocks = static_cast<int>(tables.blocks.size());
  for (std::size_t block_index = 0; block_index < tables.blocks.size();
       ++block_index) {
    const auto& block = tables.blocks[block_index];
    for (const PrimeTable& table : block) {
      const int prime = table.prime;
      int index = prime;
      if (denominator != 0 && denominator % prime != 0) {
        const int inverse =
            power_mod(normalized_mod(denominator, prime), prime - 2, prime);
        index = multiply_mod(normalized_mod(numerator, prime), inverse, prime);
      }
      add_calibration_symbol(
          candidate, table, table.symbols[static_cast<std::size_t>(index)],
          static_cast<int>(block_index));
    }
    candidate.standardized_block_scores[block_index] /=
        std::sqrt(static_cast<double>(block.size()));
  }
  return candidate;
}

struct FullCalibrationSummary {
  std::uint64_t population = 0;
  std::vector<Candidate> controls;
  std::vector<std::uint64_t> control_ranks;
  std::vector<Candidate> finalists;
  double seconds = 0.0;
};

FullCalibrationSummary scan_full_worst_block(
    std::int64_t numerator_bound, std::int64_t denominator_bound,
    const Tables& tables, std::int64_t parameter_scale,
    const std::vector<std::pair<std::int64_t, std::int64_t>>& control_pairs,
    int finalist_count) {
  const auto started = std::chrono::steady_clock::now();
  FullCalibrationSummary summary;
  for (const auto& pair : control_pairs)
    summary.controls.push_back(
        score_calibration_pair(pair.first, pair.second, tables));
  summary.control_ranks.assign(summary.controls.size(), 1);
  std::vector<Candidate> heap;
  heap.reserve(static_cast<std::size_t>(finalist_count + 1));

  auto observe = [&](Candidate candidate) {
    ++summary.population;
    for (std::size_t i = 0; i < summary.controls.size(); ++i)
      if (calibration_better(candidate, summary.controls[i]))
        ++summary.control_ranks[i];
    if (static_cast<int>(heap.size()) < finalist_count) {
      heap.push_back(std::move(candidate));
      std::push_heap(heap.begin(), heap.end(), calibration_better);
    } else if (calibration_better(candidate, heap.front())) {
      std::pop_heap(heap.begin(), heap.end(), calibration_better);
      heap.back() = std::move(candidate);
      std::push_heap(heap.begin(), heap.end(), calibration_better);
    }
  };

  Candidate infinity = score_calibration_pair(1, 0, tables);
  observe(infinity);

  struct IncrementalTable {
    const PrimeTable* table = nullptr;
    int block_index = 0;
    int residue = 0;
    int step = 0;
  };
  std::vector<IncrementalTable> incremental;
  for (std::size_t block_index = 0; block_index < tables.blocks.size();
       ++block_index)
    for (const PrimeTable& table : tables.blocks[block_index])
      incremental.push_back({&table, static_cast<int>(block_index), 0, 0});

  for (std::int64_t denominator = 1; denominator <= denominator_bound;
       ++denominator) {
    for (IncrementalTable& entry : incremental) {
      const int prime = entry.table->prime;
      if (denominator % prime == 0) {
        entry.step = -1;
        entry.residue = prime;
      } else {
        const int inverse = power_mod(normalized_mod(denominator, prime),
                                      prime - 2, prime);
        entry.step = multiply_mod(normalized_mod(parameter_scale, prime),
                                  inverse, prime);
        entry.residue = multiply_mod(
            multiply_mod(normalized_mod(parameter_scale, prime),
                         normalized_mod(-numerator_bound, prime), prime),
            inverse, prime);
      }
    }
    for (std::int64_t numerator = -numerator_bound;
         numerator <= numerator_bound; ++numerator) {
      if (std::gcd(numerator < 0 ? -numerator : numerator, denominator) == 1) {
        Candidate candidate;
        candidate.numerator = numerator * parameter_scale;
        candidate.denominator = denominator;
        const std::int64_t common = std::gcd(
            candidate.numerator < 0 ? -candidate.numerator : candidate.numerator,
            candidate.denominator);
        candidate.numerator /= common;
        candidate.denominator /= common;
        candidate.height = std::max(
            candidate.numerator < 0 ? -candidate.numerator : candidate.numerator,
            candidate.denominator);
        candidate.completed_blocks = static_cast<int>(tables.blocks.size());
        for (const IncrementalTable& entry : incremental) {
          const LocalSymbol& symbol = entry.table->symbols[
              static_cast<std::size_t>(entry.step < 0 ? entry.table->prime
                                                     : entry.residue)];
          add_calibration_symbol(candidate, *entry.table, symbol,
                                 entry.block_index);
        }
        for (std::size_t block_index = 0; block_index < tables.blocks.size();
             ++block_index)
          candidate.standardized_block_scores[block_index] /=
              std::sqrt(static_cast<double>(tables.blocks[block_index].size()));
        observe(std::move(candidate));
      }
      for (IncrementalTable& entry : incremental) {
        if (entry.step < 0) continue;
        entry.residue += entry.step;
        if (entry.residue >= entry.table->prime)
          entry.residue -= entry.table->prime;
      }
    }
  }
  std::sort(heap.begin(), heap.end(), calibration_better);
  summary.finalists = std::move(heap);
  summary.seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started)
                        .count();
  return summary;
}

bool insert_pareto(std::vector<Candidate>& frontier, const Candidate& candidate,
                   int cap) {
  for (const Candidate& retained : frontier)
    if (dominates(retained, candidate)) return false;
  frontier.erase(
      std::remove_if(frontier.begin(), frontier.end(),
                     [&](const Candidate& retained) {
                       return dominates(candidate, retained);
                     }),
      frontier.end());
  frontier.push_back(candidate);
  if (static_cast<int>(frontier.size()) > cap) {
    std::sort(frontier.begin(), frontier.end(), better);
    frontier.resize(static_cast<std::size_t>(cap));
  }
  return true;
}

void add_symbol(Candidate& candidate, const LocalSymbol& symbol,
                std::int64_t& block_score) {
  if (symbol.good) {
    block_score += symbol.score_units;
    ++candidate.good_primes;
  } else {
    ++candidate.bad_primes;
  }
}

Candidate score_later_block(Candidate candidate,
                            const std::vector<PrimeTable>& tables,
                            std::int64_t parameter_scale) {
  std::int64_t block_score = 0;
  for (const PrimeTable& table : tables) {
    const int prime = table.prime;
    int index = prime;
    if (candidate.denominator != 0 && candidate.denominator % prime != 0) {
      const int denominator = normalized_mod(candidate.denominator, prime);
      const int inverse = power_mod(denominator, prime - 2, prime);
      index = multiply_mod(
          multiply_mod(normalized_mod(parameter_scale, prime),
                       normalized_mod(candidate.numerator, prime), prime),
          inverse, prime);
    }
    add_symbol(candidate, table.symbols[static_cast<std::size_t>(index)],
               block_score);
  }
  candidate.block_scores[candidate.completed_blocks++] = block_score;
  return candidate;
}

std::vector<Candidate> flatten_and_sort(
    const std::map<std::int64_t, std::vector<Candidate>>& buckets) {
  std::vector<Candidate> result;
  std::size_t size = 0;
  for (const auto& entry : buckets) size += entry.second.size();
  result.reserve(size);
  for (const auto& entry : buckets)
    result.insert(result.end(), entry.second.begin(), entry.second.end());
  std::sort(result.begin(), result.end(), better);
  return result;
}

std::pair<std::vector<Candidate>, StageSummary> scan_first_stage(
    std::int64_t numerator_bound, std::int64_t denominator_bound,
    std::int64_t bucket_width, int cap,
    const std::vector<PrimeTable>& tables, std::int64_t parameter_scale) {
  const auto started = std::chrono::steady_clock::now();
  StageSummary summary;
  summary.stage = 1;
  summary.cap_per_bucket = cap;
  for (const PrimeTable& table : tables) summary.primes.push_back(table.prime);
  std::map<std::int64_t, std::vector<Candidate>> buckets;

  Candidate infinity;
  infinity.numerator = 1;
  infinity.denominator = 0;
  infinity.height = 1;
  infinity = score_later_block(infinity, tables, parameter_scale);
  ++summary.population_scored;
  summary.pareto_insertions += insert_pareto(buckets[0], infinity, cap);

  for (std::int64_t denominator = 1; denominator <= denominator_bound;
       ++denominator) {
    std::vector<int> steps;
    std::vector<int> residues;
    steps.reserve(tables.size());
    residues.reserve(tables.size());
    for (const PrimeTable& table : tables) {
      const int prime = table.prime;
      const int scale = normalized_mod(parameter_scale, prime);
      if (scale == 0)
        throw std::runtime_error(
            "parameter scale is not invertible modulo a table prime");
      if (denominator % prime == 0) {
        steps.push_back(-1);
        residues.push_back(prime);
      } else {
        const int inverse = power_mod(normalized_mod(denominator, prime),
                                      prime - 2, prime);
        steps.push_back(multiply_mod(scale, inverse, prime));
        residues.push_back(multiply_mod(
            multiply_mod(scale, normalized_mod(-numerator_bound, prime), prime),
            inverse, prime));
      }
    }

    for (std::int64_t numerator = -numerator_bound;
         numerator <= numerator_bound; ++numerator) {
      if (std::gcd(numerator < 0 ? -numerator : numerator, denominator) == 1) {
        Candidate candidate;
        candidate.numerator = numerator;
        candidate.denominator = denominator;
        candidate.height = std::max(numerator < 0 ? -numerator : numerator,
                                    denominator);
        std::int64_t block_score = 0;
        for (std::size_t i = 0; i < tables.size(); ++i)
          add_symbol(candidate, tables[i].symbols[static_cast<std::size_t>(
                                    steps[i] < 0 ? tables[i].prime : residues[i])],
                     block_score);
        candidate.block_scores[0] = block_score;
        candidate.completed_blocks = 1;
        ++summary.population_scored;
        const std::int64_t bucket = (candidate.height - 1) / bucket_width;
        summary.pareto_insertions +=
            insert_pareto(buckets[bucket], candidate, cap);
      }
      for (std::size_t i = 0; i < tables.size(); ++i) {
        if (steps[i] >= 0) {
          residues[i] += steps[i];
          if (residues[i] >= tables[i].prime) residues[i] -= tables[i].prime;
        }
      }
    }
  }
  std::vector<Candidate> retained = flatten_and_sort(buckets);
  summary.bucket_count = buckets.size();
  summary.retained_count = retained.size();
  summary.seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started)
                        .count();
  return {std::move(retained), summary};
}

std::pair<std::vector<Candidate>, StageSummary> scan_later_stage(
    const std::vector<Candidate>& input, std::int64_t bucket_width, int cap,
    int stage, const std::vector<PrimeTable>& tables,
    std::int64_t parameter_scale) {
  const auto started = std::chrono::steady_clock::now();
  StageSummary summary;
  summary.stage = stage;
  summary.cap_per_bucket = cap;
  for (const PrimeTable& table : tables) summary.primes.push_back(table.prime);
  std::map<std::int64_t, std::vector<Candidate>> buckets;
  for (const Candidate& unscored : input) {
    Candidate candidate =
        score_later_block(unscored, tables, parameter_scale);
    ++summary.population_scored;
    const std::int64_t bucket = (candidate.height - 1) / bucket_width;
    summary.pareto_insertions += insert_pareto(buckets[bucket], candidate, cap);
  }
  std::vector<Candidate> retained = flatten_and_sort(buckets);
  summary.bucket_count = buckets.size();
  summary.retained_count = retained.size();
  summary.seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started)
                        .count();
  return {std::move(retained), summary};
}

std::string json_escape(const std::string& text) {
  std::ostringstream output;
  for (const char character : text) {
    if (character == '\\' || character == '"') output << '\\';
    output << character;
  }
  return output.str();
}

void write_calibration_candidate(std::ostream& output,
                                 const Candidate& candidate) {
  output << "{\"parameter\": \"" << candidate.parameter()
         << "\", \"projective_pair\": [" << candidate.numerator << ", "
         << candidate.denominator << "], \"projective_height\": "
         << candidate.height << ", \"standardized_block_signals\": [";
  for (int block = 0; block < candidate.completed_blocks; ++block) {
    if (block) output << ", ";
    output << std::setprecision(17)
           << candidate.standardized_block_scores[block];
  }
  output << "], \"worst_block_signal\": " << std::setprecision(17)
         << candidate.minimum_standardized_block_score()
         << ", \"mean_block_signal\": "
         << candidate.mean_standardized_block_score()
         << ", \"raw_block_score_units_1e12\": [";
  for (int block = 0; block < candidate.completed_blocks; ++block) {
    if (block) output << ", ";
    output << candidate.block_scores[block];
  }
  output << "], \"good_prime_count\": " << candidate.good_primes
         << ", \"bad_reduction_prime_count\": " << candidate.bad_primes
         << "}";
}

void write_full_calibration_json(
    const std::string& path, const std::string& table_path,
    const Tables& tables, std::int64_t numerator_bound,
    std::int64_t denominator_bound, std::int64_t parameter_scale,
    const FullCalibrationSummary& summary, int argc, char** argv) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot create output file: " + path);
  bool accepted = true;
  for (std::uint64_t rank : summary.control_ranks)
    accepted = accepted &&
               static_cast<long double>(rank) / summary.population <= 0.01L;
  output << "{\n"
         << "  \"schema\": \"elkies-2026-positive-control-worst-block-nagao-v1\",\n"
         << "  \"status\": \""
         << (accepted ? "PASS_POSITIVE_CONTROL_SCORING_GATE"
                      : "FAIL_DISCARD_POSITIVE_CONTROL_SCORING_METHOD")
         << "\",\n"
         << "  \"proof_boundary\": \"This is a complete bounded heuristic ranking, not a rank or Selmer bound. Singular local fibres are mean-imputed and counted.\",\n"
         << "  \"model_sha256\": \"" << tables.model_sha256 << "\",\n"
         << "  \"table_file\": \"" << json_escape(table_path) << "\",\n"
         << "  \"search\": {\"coordinate\": \"published compact t\", "
         << "\"numerator_interval\": [-" << numerator_bound << ", "
         << numerator_bound << "], \"denominator_interval\": [1, "
         << denominator_bound
         << "], \"primitive_pairs_only\": true, \"includes_infinity\": true, "
         << "\"height\": \"max(abs(a),b)\", \"height_limit\": "
         << std::max(numerator_bound, denominator_bound)
         << ", \"parameter_scale\": " << parameter_scale << "},\n"
         << "  \"scoring\": {\"prime_ensembles\": [";
  for (std::size_t block = 0; block < tables.blocks.size(); ++block) {
    if (block) output << ", ";
    output << '[';
    for (std::size_t index = 0; index < tables.blocks[block].size(); ++index) {
      if (index) output << ", ";
      output << tables.blocks[block][index].prime;
    }
    output << ']';
  }
  output << "], \"ensembles_are_pairwise_disjoint\": true, "
         << "\"per_prime_standardization\": \"center and population-standardize over good fibres of P1(F_p)\", "
         << "\"singular_fibre_policy\": \"mean imputation (standardized contribution zero)\", "
         << "\"block_normalization\": \"sum(z_p)/sqrt(number of primes in block)\", "
         << "\"primary_ranking_key\": \"minimum block signal\", "
         << "\"tie_breaker\": \"mean block signal, good primes, bad primes, height, denominator, numerator\", "
         << "\"positive_control_acceptance_threshold\": \"every disclosed fibre in top 1 percent of complete H<=10000 population\"},\n"
         << "  \"population_count\": " << summary.population << ",\n"
         << "  \"positive_controls\": [\n";
  for (std::size_t index = 0; index < summary.controls.size(); ++index) {
    output << "    {\"score\": ";
    write_calibration_candidate(output, summary.controls[index]);
    const long double fraction =
        static_cast<long double>(summary.control_ranks[index]) /
        summary.population;
    output << ", \"population_rank\": " << summary.control_ranks[index]
           << ", \"population_fraction\": " << std::setprecision(17)
           << static_cast<double>(fraction)
           << ", \"passes_top_one_percent_gate\": "
           << (fraction <= 0.01L ? "true" : "false") << "}"
           << (index + 1 == summary.controls.size() ? "\n" : ",\n");
  }
  output << "  ],\n  \"finalists\": [\n";
  for (std::size_t index = 0; index < summary.finalists.size(); ++index) {
    output << "    ";
    write_calibration_candidate(output, summary.finalists[index]);
    output << (index + 1 == summary.finalists.size() ? "\n" : ",\n");
  }
  output << "  ],\n  \"runtime_seconds\": " << std::setprecision(12)
         << summary.seconds << ",\n  \"reproducing_command\": \"";
  for (int i = 0; i < argc; ++i) {
    if (i) output << ' ';
    output << json_escape(argv[i]);
  }
  output << "\"\n}\n";
}

void write_json(const std::string& path, const std::string& table_path,
                const Tables& tables, std::int64_t numerator_bound,
                std::int64_t denominator_bound, std::int64_t bucket_width,
                const std::vector<int>& keeps,
                const std::vector<StageSummary>& stages,
                const std::vector<Candidate>& survivors, int finalist_count,
                double total_seconds, std::int64_t parameter_scale, int argc,
                char** argv) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot create output file: " + path);
  output << "{\n"
         << "  \"schema\": \"h92-q12o5867-rootless-projective-nagao-cpp-"
         << (parameter_scale == 1 ? "v1" : "skew-v2") << "\",\n"
         << "  \"status\": \"PASS_BOUNDED_HEURISTIC_PROJECTIVE_NAGAO_CPP_SIEVE\",\n"
         << "  \"proof_boundary\": \"Table scores and survival are heuristics, not rank bounds; no section was evaluated.\",\n"
         << "  \"model_sha256\": \"" << tables.model_sha256 << "\",\n"
         << "  \"table_file\": \"" << json_escape(table_path) << "\",\n"
         << "  \"score_scale\": " << tables.score_scale << ",\n"
         << "  \"search\": {\n"
         << "    \"numerator_interval\": [-" << numerator_bound << ", "
         << numerator_bound << "],\n"
         << "    \"denominator_interval\": [1, " << denominator_bound << "],\n"
         << "    \"includes_zero\": true,\n"
         << "    \"includes_infinity\": true,\n"
         << "    \"primitive_pairs_only\": true,\n"
         << "    \"parameter_chart\": \"u=parameter_scale*v\",\n"
         << "    \"parameter_scale\": " << parameter_scale << ",\n"
         << "    \"height_bucket_width\": " << bucket_width << ",\n"
         << "    \"keep_per_bucket\": [";
  for (std::size_t i = 0; i < keeps.size(); ++i) {
    if (i) output << ", ";
    output << keeps[i];
  }
  output << "]\n  },\n  \"stages\": [\n";
  for (std::size_t i = 0; i < stages.size(); ++i) {
    const StageSummary& stage = stages[i];
    output << "    {\"stage\": " << stage.stage << ", \"primes\": [";
    for (std::size_t j = 0; j < stage.primes.size(); ++j) {
      if (j) output << ", ";
      output << stage.primes[j];
    }
    output << "], \"cap_per_height_bucket\": " << stage.cap_per_bucket
           << ", \"population_scored\": " << stage.population_scored
           << ", \"height_bucket_count\": " << stage.bucket_count
           << ", \"pareto_insertions\": " << stage.pareto_insertions
           << ", \"retained_count\": " << stage.retained_count
           << ", \"runtime_seconds\": " << std::setprecision(12)
           << stage.seconds << ", \"parameters_per_second\": "
           << (stage.seconds == 0.0 ? 0.0
                                    : stage.population_scored / stage.seconds)
           << "}" << (i + 1 == stages.size() ? "\n" : ",\n");
  }
  const int written = std::min<int>(finalist_count, survivors.size());
  output << "  ],\n  \"final_survivor_count\": " << survivors.size()
         << ",\n  \"finalists\": [\n";
  for (int i = 0; i < written; ++i) {
    const Candidate& candidate = survivors[static_cast<std::size_t>(i)];
    std::int64_t actual_numerator = candidate.numerator;
    std::int64_t actual_denominator = candidate.denominator;
    if (actual_denominator != 0) {
      actual_numerator *= parameter_scale;
      const std::int64_t common =
          std::gcd(actual_numerator < 0 ? -actual_numerator : actual_numerator,
                   actual_denominator);
      actual_numerator /= common;
      actual_denominator /= common;
      if (actual_denominator < 0) {
        actual_numerator = -actual_numerator;
        actual_denominator = -actual_denominator;
      }
    }
    const std::int64_t actual_height =
        actual_denominator == 0
            ? 1
            : std::max(actual_numerator < 0 ? -actual_numerator
                                            : actual_numerator,
                       actual_denominator);
    output << "    {\"parameter\": \"";
    if (actual_denominator == 0)
      output << "infinity";
    else
      output << actual_numerator << "/" << actual_denominator;
    output << "\", \"projective_pair\": ["
           << actual_numerator << ", " << actual_denominator
           << "], \"projective_height\": " << actual_height
           << ", \"chart_parameter\": \"" << candidate.parameter()
           << "\", \"chart_projective_pair\": [" << candidate.numerator
           << ", " << candidate.denominator << "], \"chart_height\": "
           << candidate.height << ", \"block_score_units_1e12\": [";
    for (int j = 0; j < candidate.completed_blocks; ++j) {
      if (j) output << ", ";
      output << candidate.block_scores[j];
    }
    output << "], \"total_score_units_1e12\": " << candidate.total_score()
           << ", \"total_score\": " << std::fixed << std::setprecision(12)
           << static_cast<double>(candidate.total_score()) / tables.score_scale
           << ", \"good_prime_count\": " << candidate.good_primes
           << ", \"bad_reduction_prime_count\": " << candidate.bad_primes
           << "}" << (i + 1 == written ? "\n" : ",\n");
  }
  output << "  ],\n  \"runtime_seconds\": " << std::setprecision(12)
         << total_seconds << ",\n  \"reproducing_command\": \"";
  for (int i = 0; i < argc; ++i) {
    if (i) output << ' ';
    output << json_escape(argv[i]);
  }
  output << "\"\n}\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 8 && argc != 9 && argc != 10) {
    std::cerr << "usage: " << argv[0]
              << " TABLE NUM DEN BUCKET_WIDTH KEEP1,KEEP2,... FINALISTS OUTPUT.json [PARAMETER_SCALE [CONTROL_A/B,...]]\n";
    return 2;
  }
  try {
    const std::string table_path = argv[1];
    const std::int64_t numerator_bound = std::stoll(argv[2]);
    const std::int64_t denominator_bound = std::stoll(argv[3]);
    const std::int64_t bucket_width = std::stoll(argv[4]);
    const std::vector<int> keeps = parse_keeps(argv[5]);
    const int finalists = std::stoi(argv[6]);
    const std::string output_path = argv[7];
    const std::int64_t parameter_scale = argc >= 9 ? std::stoll(argv[8]) : 1;
    if (numerator_bound < 0 || denominator_bound < 1 || bucket_width < 1 ||
        finalists < 1 || parameter_scale == 0)
      throw std::runtime_error("bounds, bucket width, and finalists are invalid");

    const Tables tables = read_tables(table_path);
    if (tables.a_degree != 8 || tables.b_degree != 12 ||
        tables.score_scale != 1000000000000LL)
      throw std::runtime_error("unexpected model degrees or score scale");
    if (keeps.size() != tables.blocks.size())
      throw std::runtime_error("one keep count is required per table block");

    if (argc == 10) {
      if (tables.blocks.size() < 3)
        throw std::runtime_error(
            "positive-control calibration requires at least three prime ensembles");
      const auto controls = parse_controls(argv[9]);
      const FullCalibrationSummary summary = scan_full_worst_block(
          numerator_bound, denominator_bound, tables, parameter_scale,
          controls, finalists);
      write_full_calibration_json(
          output_path, table_path, tables, numerator_bound, denominator_bound,
          parameter_scale, summary, argc, argv);
      bool accepted = true;
      for (std::uint64_t rank : summary.control_ranks)
        accepted = accepted &&
                   static_cast<long double>(rank) / summary.population <= 0.01L;
      std::cout << (accepted ? "PASS" : "FAIL")
                << " positive_control_population=" << summary.population
                << " control_ranks=";
      for (std::size_t index = 0; index < summary.control_ranks.size(); ++index) {
        if (index) std::cout << ',';
        std::cout << summary.control_ranks[index];
      }
      std::cout << " seconds=" << std::fixed << std::setprecision(3)
                << summary.seconds << " output=" << output_path << "\n";
      return accepted ? 0 : 3;
    }

    const auto started = std::chrono::steady_clock::now();
    auto first = scan_first_stage(numerator_bound, denominator_bound,
                                  bucket_width, keeps[0], tables.blocks[0],
                                  parameter_scale);
    std::vector<Candidate> survivors = std::move(first.first);
    std::vector<StageSummary> stages;
    stages.push_back(first.second);
    for (std::size_t block = 1; block < tables.blocks.size(); ++block) {
      auto next = scan_later_stage(survivors, bucket_width, keeps[block],
                                   static_cast<int>(block + 1),
                                   tables.blocks[block], parameter_scale);
      survivors = std::move(next.first);
      stages.push_back(next.second);
    }
    const double total_seconds = std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - started)
                                     .count();
    write_json(output_path, table_path, tables, numerator_bound,
               denominator_bound, bucket_width, keeps, stages, survivors,
               finalists, total_seconds, parameter_scale, argc, argv);
    std::cout << "PASS population=" << stages.front().population_scored
              << " survivors=" << survivors.size() << " seconds="
              << std::fixed << std::setprecision(3) << total_seconds
              << " rate=" << std::setprecision(0)
              << stages.front().population_scored / stages.front().seconds
              << "/s output=" << output_path << "\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 2;
  }
  return 0;
}
