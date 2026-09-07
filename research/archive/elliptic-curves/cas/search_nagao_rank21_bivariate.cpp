// Exhaustive rational-parameter/offset square search in Nagao's rank-21 family.
//
// For the primitive Mestre quartic attached to roots
//
//   (0, 4, 47, 352, 380, 399),
//
// write T=a/b and X=T+k/h.  Its coefficients are
//
//   c0=-23T^6+4306112T^4-5788761232T^2+197804341504,
//   c1=-18124T^4+98440030T^2-37842168672,
//   c2=46T^4-3814698T^2+6945772145,
//   c3=18124T^2-35986906,
//   c4=47342-23T^2.
//
// Hence F=b^6*h^4*H_T(T+k/h) is integral and H_T(X) is a rational
// square exactly when F is an integer square.  This program checks every
// primitive positive (a,b), every reduced k/h, and every declared k, after
// imposing the four automatically discovered clean split-multiplicative
// residue unions
//
//   T mod 5  in {1,4},       T mod 7  in {0,3,4},
//   T mod 11 in {2,9},       T mod 13 in {3,10}.
//
// The six X=T+r displayed sections are omitted.  The other six displayed
// sections and every genuine extra hit are separated by exact downstream
// code.  This bounded square search is not a rank certificate.

#include <gmpxx.h>

#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {

mpz_class evaluate(long a_long, long b_long, long h_long, long k_long,
                   long slope_long) {
  const mpz_class a = a_long;
  const mpz_class b = b_long;
  const mpz_class h = h_long;
  const mpz_class k = k_long;
  const mpz_class slope = slope_long;
  const mpz_class a2 = a * a;
  const mpz_class a4 = a2 * a2;
  const mpz_class a6 = a4 * a2;
  const mpz_class b2 = b * b;
  const mpz_class b4 = b2 * b2;
  const mpz_class b6 = b4 * b2;
  const mpz_class h2 = h * h;
  const mpz_class h3 = h2 * h;
  const mpz_class h4 = h2 * h2;
  const mpz_class z = slope * a * h + b * k;
  const mpz_class z2 = z * z;
  const mpz_class z3 = z2 * z;
  const mpz_class z4 = z2 * z2;

  const mpz_class c0 = -23 * a6 + 4306112 * a4 * b2
      - mpz_class("5788761232") * a2 * b4
      + mpz_class("197804341504") * b6;
  const mpz_class c1 = -18124 * a4 + mpz_class("98440030") * a2 * b2
      - mpz_class("37842168672") * b4;
  const mpz_class c2 = 46 * a4 - mpz_class("3814698") * a2 * b2
      + mpz_class("6945772145") * b4;
  const mpz_class c3 = 18124 * a2 - mpz_class("35986906") * b2;
  const mpz_class c4 = 47342 * b2 - 23 * a2;
  return c0 * h4 + c1 * b * h3 * z + c2 * h2 * z2
      + c3 * b * h * z3 + c4 * z4;
}

bool in_residue_union(long a, long b, long prime,
                      const std::initializer_list<long>& residues) {
  const long denominator = b % prime;
  if (denominator == 0) return false;
  long numerator = a % prime;
  if (numerator < 0) numerator += prime;
  for (long residue : residues) {
    if (numerator == (residue * denominator) % prime) return true;
  }
  return false;
}

bool conductor_friendly(long a, long b) {
  return in_residue_union(a, b, 5, {1, 4})
      && in_residue_union(a, b, 7, {0, 3, 4})
      && in_residue_union(a, b, 11, {2, 9})
      && in_residue_union(a, b, 13, {3, 10});
}

bool displayed_section(long slope, long h, long k) {
  if (slope != 1 && slope != -1) return false;
  if (h != 1) return false;
  return k == 0 || k == 4 || k == 47 || k == 352 || k == 380 || k == 399;
}

long parse_long(const char* text, const char* name) {
  char* end = nullptr;
  const long result = std::strtol(text, &end, 10);
  if (end == text || *end != '\0') {
    throw std::invalid_argument(std::string("invalid ") + name);
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 7 && argc != 8 && argc != 10) {
    std::cerr << "usage: search_nagao_rank21_bivariate "
                 "A_MAX B_MAX H_MAX K_MIN K_MAX A_MIN "
                 "[LOCAL_FILTER [SLOPE_MIN SLOPE_MAX]]\n";
    return 2;
  }
  const long a_max = parse_long(argv[1], "A_MAX");
  const long b_max = parse_long(argv[2], "B_MAX");
  const long h_max = parse_long(argv[3], "H_MAX");
  const long k_min = parse_long(argv[4], "K_MIN");
  const long k_max = parse_long(argv[5], "K_MAX");
  const long a_min = parse_long(argv[6], "A_MIN");
  const bool local_filter = argc == 7 || parse_long(argv[7], "LOCAL_FILTER") != 0;
  const long slope_min = argc == 10 ? parse_long(argv[8], "SLOPE_MIN") : 1;
  const long slope_max = argc == 10 ? parse_long(argv[9], "SLOPE_MAX") : 1;
  if (a_min <= 0 || a_max < a_min || b_max <= 0 || h_max <= 0
      || k_max < k_min || slope_max < slope_min) {
    std::cerr << "require positive A/B/H bounds, A_MIN <= A_MAX, "
                 "K_MIN <= K_MAX\n";
    return 2;
  }

  const auto started = std::chrono::steady_clock::now();
  unsigned long long parameter_count = 0;
  unsigned long long tested = 0;
  unsigned long long hits = 0;
  for (long b = 1; b <= b_max; ++b) {
    for (long a = a_min; a <= a_max; ++a) {
      if (std::gcd(a, b) != 1
          || (local_filter && !conductor_friendly(a, b))) continue;
      ++parameter_count;
      for (long slope = slope_min; slope <= slope_max; ++slope) {
        for (long h = 1; h <= h_max; ++h) {
          mpz_class f0 = evaluate(a, b, h, k_min, slope);
          const mpz_class f1 = evaluate(a, b, h, k_min + 1, slope);
          const mpz_class f2 = evaluate(a, b, h, k_min + 2, slope);
          const mpz_class f3 = evaluate(a, b, h, k_min + 3, slope);
          const mpz_class f4 = evaluate(a, b, h, k_min + 4, slope);
          mpz_class d0 = f0;
          mpz_class d1 = f1 - f0;
          mpz_class d2 = f2 - 2 * f1 + f0;
          mpz_class d3 = f3 - 3 * f2 + 3 * f1 - f0;
          const mpz_class d4 = f4 - 4 * f3 + 6 * f2 - 4 * f1 + f0;
          for (long k = k_min; k <= k_max; ++k) {
            if (std::gcd(std::labs(k), h) == 1) {
              ++tested;
              if (!displayed_section(slope, h, k) && d0 >= 0
                  && mpz_perfect_square_p(d0.get_mpz_t())) {
                mpz_class root;
                mpz_sqrt(root.get_mpz_t(), d0.get_mpz_t());
                std::cout << a << '\t' << b << '\t' << slope << '\t' << h
                          << '\t' << k << '\t' << root << '\n';
                ++hits;
              }
            }
            d0 += d1;
            d1 += d2;
            d2 += d3;
            d3 += d4;
          }
        }
      }
    }
  }
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  std::cerr << "parameters=" << parameter_count << " tested=" << tested
            << " hits=" << hits << " seconds=" << seconds << '\n';
  return 0;
}
