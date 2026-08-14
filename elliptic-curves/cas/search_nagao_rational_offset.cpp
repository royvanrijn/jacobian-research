// Exhaustive rational-u, integral-offset search in Nagao's rank-13 family.
//
// With u=a/b in lowest terms, T=(23550*b^2-a^2)/(2*a*b), and X=T+k/h,
// symbolic expansion of Nagao's quartic gives
//
//   H_T(X) = P(a/b,k)/(4*(a/b)^4).
//
// After homogenizing P to degree eight, an integer square
// W^2=h^4*b^8*P(a/b,k/h) yields Y=W/(2*a^2*b^2*h^2).  This program checks
// every primitive positive (a,b), every reduced k/h, and every declared k.
// The six X=T+r Mestre sections are omitted; all other hits are reported for
// exact downstream classification.  No independence or rank claim is made.

#include <gmpxx.h>

#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {

mpz_class term(const char* coefficient,
               const std::array<mpz_class, 9>& a_power,
               const std::array<mpz_class, 9>& b_power,
               int degree) {
  return mpz_class(coefficient) * a_power[degree] * b_power[8 - degree];
}

std::array<mpz_class, 5> coefficients(long a_long, long b_long) {
  const mpz_class a = a_long;
  const mpz_class b = b_long;
  std::array<mpz_class, 9> ap;
  std::array<mpz_class, 9> bp;
  ap[0] = bp[0] = 1;
  for (int index = 1; index <= 8; ++index) {
    ap[index] = ap[index - 1] * a;
    bp[index] = bp[index - 1] * b;
  }

  std::array<mpz_class, 5> c;
  c[4] = term("9", ap, bp, 6) + term("423900", ap, bp, 4)
       + term("4991422500", ap, bp, 2);
  c[3] = term("-18", ap, bp, 7) + term("-2700", ap, bp, 6)
       + term("-423900", ap, bp, 5) + term("-128436840", ap, bp, 4)
       + term("9982845000", ap, bp, 3)
       + term("-1497426750000", ap, bp, 2)
       + term("235095999750000", ap, bp, 1);
  c[2] = term("9", ap, bp, 8) + term("4050", ap, bp, 7)
       + term("820050", ap, bp, 6) + term("97277760", ap, bp, 5)
       + term("-21781294044", ap, bp, 4)
       + term("-2290891248000", ap, bp, 3)
       + term("454801780125000", ap, bp, 2)
       + term("-52896599943750000", ap, bp, 1)
       + term("2768255397056250000", ap, bp, 0);
  c[1] = term("-1350", ap, bp, 8) + term("-820050", ap, bp, 7)
       + term("-94110480", ap, bp, 6) + term("31110626544", ap, bp, 5)
       + term("4792288971600", ap, bp, 4)
       + term("-732655255111200", ap, bp, 3)
       + term("-52193907484200000", ap, bp, 2)
       + term("10710581921943750000", ap, bp, 1)
       + term("-415238309558437500000", ap, bp, 0);
  c[0] = term("112225", ap, bp, 8) + term("46738530", ap, bp, 7)
       + term("-2973418919", ap, bp, 6)
       + term("-2733206451300", ap, bp, 5)
       + term("32192690719900", ap, bp, 4)
       + term("64367011928115000", ap, bp, 3)
       + term("-1649065566024697500", ap, bp, 2)
       + term("-610446746510853750000", ap, bp, 1)
       + term("34518606881626406250000", ap, bp, 0);
  return c;
}

mpz_class evaluate(const std::array<mpz_class, 5>& c, long k) {
  mpz_class result = c[4];
  for (int index = 3; index >= 0; --index) result = result * k + c[index];
  return result;
}

bool displayed_section(long k) {
  return k == 0 || k == 25 || k == 57 || k == 104 || k == 116 || k == 148;
}

long parse_long(const char* text, const char* name) {
  char* end = nullptr;
  const long result = std::strtol(text, &end, 10);
  if (end == text || *end != '\0') throw std::invalid_argument(std::string("invalid ") + name);
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 7) {
    std::cerr << "usage: search_nagao_rational_offset A_MAX B_MAX H_MAX K_MIN K_MAX A_MIN\n";
    return 2;
  }
  const long a_max = parse_long(argv[1], "A_MAX");
  const long b_max = parse_long(argv[2], "B_MAX");
  const long h_max = parse_long(argv[3], "H_MAX");
  const long k_min = parse_long(argv[4], "K_MIN");
  const long k_max = parse_long(argv[5], "K_MAX");
  const long a_min = parse_long(argv[6], "A_MIN");
  if (a_min <= 0 || a_max < a_min || b_max <= 0 || h_max <= 0 || k_max < k_min) {
    std::cerr << "require positive A/B/H bounds, A_MIN <= A_MAX, K_MIN <= K_MAX\n";
    return 2;
  }

  const auto started = std::chrono::steady_clock::now();
  unsigned long long parameter_count = 0;
  unsigned long long tested = 0;
  unsigned long long hits = 0;
  for (long b = 1; b <= b_max; ++b) {
    for (long a = a_min; a <= a_max; ++a) {
      if (std::gcd(a, b) != 1) continue;
      ++parameter_count;
      const auto c = coefficients(a, b);
      for (long h = 1; h <= h_max; ++h) {
        std::array<mpz_class, 5> scaled;
        mpz_class h_power = 1;
        for (int index = 4; index >= 0; --index) {
          scaled[index] = c[index] * h_power;
          h_power *= h;
        }
        const mpz_class f0 = evaluate(scaled, k_min);
        const mpz_class f1 = evaluate(scaled, k_min + 1);
        const mpz_class f2 = evaluate(scaled, k_min + 2);
        const mpz_class f3 = evaluate(scaled, k_min + 3);
        const mpz_class f4 = evaluate(scaled, k_min + 4);
        mpz_class d0 = f0;
        mpz_class d1 = f1 - f0;
        mpz_class d2 = f2 - 2 * f1 + f0;
        mpz_class d3 = f3 - 3 * f2 + 3 * f1 - f0;
        const mpz_class d4 = f4 - 4 * f3 + 6 * f2 - 4 * f1 + f0;

        for (long k = k_min; k <= k_max; ++k) {
          if (std::gcd(std::labs(k), h) == 1) {
            ++tested;
            if (!(h == 1 && displayed_section(k)) && d0 >= 0
                && mpz_perfect_square_p(d0.get_mpz_t())) {
              mpz_class root;
              mpz_sqrt(root.get_mpz_t(), d0.get_mpz_t());
              std::cout << a << '\t' << b << '\t' << h << '\t' << k
                        << '\t' << root << '\n';
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
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  std::cerr << "parameters=" << parameter_count << " tested=" << tested
            << " hits=" << hits << " seconds=" << seconds << '\n';
  return 0;
}
