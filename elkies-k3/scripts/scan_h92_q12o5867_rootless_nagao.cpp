// Fast table-lookup scan for the exact q12/orbit5867 rootless R17 family.
//
// The Python companion exports every point of P^1(F_p).  This hot loop does no
// coefficient or section evaluation: it enumerates signed primitive (a:b),
// including zero and infinity, and exactly reproduces the staged bucketed
// Pareto retention used by search_h92_q12o5867_rootless_nagao.py.
//
// Usage:
//   scan TABLE NUM DEN BUCKET_WIDTH KEEP1,KEEP2,... FINALISTS OUTPUT.json
//        [PARAMETER_SCALE]
//
// The optional scale searches the rational chart u=PARAMETER_SCALE*v while
// retaining v in balanced height buckets.  It must be nonzero and invertible
// modulo every table prime.  This is useful when the exact binary A_8,B_12
// model is badly conditioned in its inherited projective coordinate.

#include <algorithm>
#include <array>
#include <chrono>
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
  if (argc != 8 && argc != 9) {
    std::cerr << "usage: " << argv[0]
              << " TABLE NUM DEN BUCKET_WIDTH KEEP1,KEEP2,... FINALISTS OUTPUT.json [PARAMETER_SCALE]\n";
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
    const std::int64_t parameter_scale = argc == 9 ? std::stoll(argv[8]) : 1;
    if (numerator_bound < 0 || denominator_bound < 1 || bucket_width < 1 ||
        finalists < 1 || parameter_scale == 0)
      throw std::runtime_error("bounds, bucket width, and finalists are invalid");

    const Tables tables = read_tables(table_path);
    if (tables.a_degree != 8 || tables.b_degree != 12 ||
        tables.score_scale != 1000000000000LL)
      throw std::runtime_error("unexpected model degrees or score scale");
    if (keeps.size() != tables.blocks.size())
      throw std::runtime_error("one keep count is required per table block");

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
