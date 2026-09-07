// Closed local-score sweep for four split-infinity six-root Mestre families.
//
// Usage: score_mestre_dsquare_four NUMERATOR_BOUND DENOMINATOR_BOUND KEEP
//        [FAMILY_INDEX [NEAR_SQRT_WINDOW [U_MIN U_MAX]]]
// The output is heuristic ordering only.  Exact conductors and point
// independence are handled by search_mestre_dsquare_four.py.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using i64 = std::int64_t;

struct Family {
  i64 c_num;
  i64 c_den;
  std::array<const char *, 9> a;
  std::array<const char *, 13> b;
};

static const std::array<Family, 4> FAMILIES = {{
    {84878, 3,
     {"-6335506384367705667", "0", "-187965968101634304", "0",
      "16096946597856", "0", "215550720", "0", "-34992"},
     {"5398702003432801270887248574", "0",
      "920235175155986372785584432", "0", "-265569147155331220745484",
      "0", "39076108758158185728", "0", "-1846021419819648", "0",
      "-23279477760", "0", "2519424"}},
    {196250, 3,
     {"-2849935192246345367547", "0", "-453784306113509400", "0",
      "51920932719456", "0", "-538876800", "0", "-34992"},
     {"58522711676291264158258068274986", "0",
      "15989934323715193045999944000", "0",
      "-2686936627228512207246444", "0", "171374816542764794400", "0",
      "-6279655654021248", "0", "58198694400", "0", "2519424"}},
    {39146, 1,
     {"-3364483917221298387", "0", "-335524216685184", "0",
      "139074579936", "0", "-6531840", "0", "-432"},
     {"2368221563410632931499846766", "0", "418411864834246730627568",
      "0", "-221282839923690745836", "0", "32526162664756992", "0",
      "-2557747746432", "0", "78382080", "0", "3456"}},
    {55950, 1,
     {"-1704784199738338039707", "0", "-297408087083885400", "0",
      "31435815302496", "0", "-576201600", "0", "-34992"},
     {"26924764436923313334978849778506", "0",
      "7986062326126687827971716800", "0", "-1177858739105628340438764",
      "0", "91628387901857268000", "0", "-4163605746749568", "0",
      "62229772800", "0", "2519424"}},
}};

struct Candidate {
  int family;
  int numerator;
  int denominator;
  double score;
  int good;
};

static i64 mod_pow(i64 base, i64 exponent, i64 modulus) {
  i64 answer = 1;
  base %= modulus;
  while (exponent) {
    if (exponent & 1) answer = answer * base % modulus;
    base = base * base % modulus;
    exponent >>= 1;
  }
  return answer;
}

static i64 inverse(i64 value, i64 prime) {
  return mod_pow((value % prime + prime) % prime, prime - 2, prime);
}

static i64 decimal_mod(const char *text, int prime) {
  bool negative = *text == '-';
  if (negative) ++text;
  i64 value = 0;
  while (*text) {
    value = (10 * value + (*text - '0')) % prime;
    ++text;
  }
  return negative && value ? prime - value : value;
}

template <std::size_t N>
static i64 evaluate(const std::array<const char *, N> &coefficients,
                    i64 argument, int prime) {
  i64 value = 0;
  for (int index = static_cast<int>(N) - 1; index >= 0; --index) {
    value = (value * argument + decimal_mod(coefficients[index], prime)) % prime;
  }
  return value;
}

static int legendre(i64 value, int prime) {
  value = (value % prime + prime) % prime;
  if (!value) return 0;
  return mod_pow(value, (prime - 1) / 2, prime) == 1 ? 1 : -1;
}

static std::vector<int> primes(int lower, int upper) {
  std::vector<bool> sieve(upper + 1, true);
  sieve[0] = sieve[1] = false;
  for (int p = 2; p * p <= upper; ++p)
    if (sieve[p])
      for (int multiple = p * p; multiple <= upper; multiple += p)
        sieve[multiple] = false;
  std::vector<int> answer;
  for (int p = lower; p <= upper; ++p)
    if (sieve[p]) answer.push_back(p);
  return answer;
}

static Candidate score_candidate(int family_index, int numerator,
                                 int denominator,
                                 const std::vector<int> &prime_list) {
  const Family &family = FAMILIES[family_index];
  double score = 0.0;
  int good = 0;
  for (int prime : prime_list) {
    i64 den = denominator % prime;
    i64 num = numerator % prime;
    i64 c_den = family.c_den % prime;
    if (!den || !num || !c_den) continue;
    i64 u = num * inverse(den, prime) % prime;
    i64 c = (family.c_num % prime) * inverse(c_den, prime) % prime;
    i64 t_den = 2 * u % prime;
    if (!t_den) continue;
    i64 t = (c - u * u) % prime;
    if (t < 0) t += prime;
    t = t * inverse(t_den, prime) % prime;
    i64 coefficient_a = evaluate(family.a, t, prime);
    i64 coefficient_b = evaluate(family.b, t, prime);
    i64 discriminant_core =
        (4 * coefficient_a % prime * coefficient_a % prime * coefficient_a
         + 27 * coefficient_b % prime * coefficient_b) % prime;
    if (!discriminant_core) continue;
    int character_sum = 0;
    for (i64 x = 0; x < prime; ++x) {
      i64 rhs = (x * x % prime * x + coefficient_a * x + coefficient_b) % prime;
      character_sum += legendre(rhs, prime);
    }
    int trace = -character_sum;
    score += static_cast<double>(2 - trace) /
             static_cast<double>(prime + 1 - trace);
    ++good;
  }
  return {family_index, numerator, denominator, score, good};
}

int main(int argc, char **argv) {
  if (argc < 4 || argc > 8 || argc == 7) {
    std::cerr << "usage: " << argv[0]
              << " NUMERATOR_BOUND DENOMINATOR_BOUND KEEP "
                 "[FAMILY_INDEX [NEAR_SQRT_WINDOW [U_MIN U_MAX]]]\n";
    return 2;
  }
  const int numerator_bound = std::stoi(argv[1]);
  const int denominator_bound = std::stoi(argv[2]);
  const int keep = std::stoi(argv[3]);
  const int selected_family = argc >= 5 ? std::stoi(argv[4]) : -1;
  const int near_sqrt_window = argc == 6 ? std::stoi(argv[5]) : 0;
  const int u_min = argc == 8 ? std::stoi(argv[6]) : 0;
  const int u_max = argc == 8 ? std::stoi(argv[7]) : 0;
  if (numerator_bound < 1 || denominator_bound < 1 || keep < 1) return 2;
  if (selected_family < -1 || selected_family >= 4) return 2;
  if (near_sqrt_window < 0 || (near_sqrt_window && selected_family < 0)) return 2;
  if (argc == 8 && (u_min < 0 || u_max < u_min)) return 2;
  const std::vector<int> prime_list = primes(11, 251);
  std::vector<std::pair<int, int>> rationals;
  for (int denominator = 1; denominator <= denominator_bound; ++denominator) {
    int numerator_begin = 1;
    int numerator_end = numerator_bound;
    if (argc == 8) {
      numerator_begin = std::max(numerator_begin, u_min * denominator);
      numerator_end = std::min(numerator_end, u_max * denominator);
    }
    if (near_sqrt_window) {
      const Family &family = FAMILIES[selected_family];
      const double center =
          std::sqrt(static_cast<double>(family.c_num) / family.c_den) * denominator;
      numerator_begin = std::max(1, static_cast<int>(std::floor(center)) - near_sqrt_window);
      numerator_end = std::min(numerator_bound,
                               static_cast<int>(std::ceil(center)) + near_sqrt_window);
    }
    for (int numerator = numerator_begin; numerator <= numerator_end; ++numerator)
      if (std::gcd(numerator, denominator) == 1)
        rationals.emplace_back(numerator, denominator);
  }

  std::cout << "MESTRE_DSQUARE_FOUR_SCORE_V1\n";
  std::cout << "PRIMES";
  for (int prime : prime_list) std::cout << ' ' << prime;
  std::cout << "\n";
  std::cout << "DOMAIN " << numerator_bound << ' ' << denominator_bound << ' '
            << rationals.size() << "\n";
  const int family_begin = selected_family < 0 ? 0 : selected_family;
  const int family_end = selected_family < 0 ? 4 : selected_family + 1;
  for (int family = family_begin; family < family_end; ++family) {
    std::vector<Candidate> candidates(rationals.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (std::size_t index = 0; index < rationals.size(); ++index) {
      candidates[index] = score_candidate(family, rationals[index].first,
                                          rationals[index].second, prime_list);
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate &left,
                                                       const Candidate &right) {
      if (left.score != right.score) return left.score > right.score;
      if (left.good != right.good) return left.good > right.good;
      if (left.denominator != right.denominator)
        return left.denominator < right.denominator;
      return left.numerator < right.numerator;
    });
    const int retained = std::min<int>(keep, candidates.size());
    for (int rank = 0; rank < retained; ++rank) {
      const Candidate &candidate = candidates[rank];
      std::cout.precision(17);
      std::cout << "C " << family << ' ' << rank + 1 << ' '
                << candidate.numerator << ' ' << candidate.denominator << ' '
                << candidate.score << ' ' << candidate.good << "\n";
    }
  }
  std::cout << "DONE\n";
  return 0;
}
