// Exhaustive integral-offset search in Nagao's rank-13 base change.
//
// For T=(23550-u^2)/(2u), substitute X=T+k in Nagao's primitive quartic
// H_T(X).  Exact symbolic expansion gives
//
//     (2u)^6 H_T(T+k) = 16 u^2 P(u,k).
//
// Thus an integral solution w^2=P(u,k) gives the exact quartic point
// (X,Y)=(T+k,w/(2u^2)).  This program exhausts a rectangular integer box.
// It reports only nonnegative square roots and omits the six displayed
// Mestre sections k=0,25,57,104,116,148.  It makes no rank claim.

#include <gmpxx.h>

#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

mpz_class evaluate(const std::array<mpz_class, 5>& coefficients, long k) {
  mpz_class answer = coefficients[4];
  for (int index = 3; index >= 0; --index) {
    answer *= k;
    answer += coefficients[index];
  }
  return answer;
}

std::array<mpz_class, 5> coefficients(long u_long) {
  const mpz_class u = u_long;
  std::array<mpz_class, 9> power;
  power[0] = 1;
  for (int index = 1; index <= 8; ++index) power[index] = power[index - 1] * u;

  std::array<mpz_class, 5> c;
  c[4] = 9 * power[6] + 423900 * power[4]
       + mpz_class("4991422500") * power[2];
  c[3] = -18 * power[7] - 2700 * power[6] - 423900 * power[5]
       - 128436840 * power[4] + mpz_class("9982845000") * power[3]
       - mpz_class("1497426750000") * power[2]
       + mpz_class("235095999750000") * power[1];
  c[2] = 9 * power[8] + 4050 * power[7] + 820050 * power[6]
       + 97277760 * power[5] - mpz_class("21781294044") * power[4]
       - mpz_class("2290891248000") * power[3]
       + mpz_class("454801780125000") * power[2]
       - mpz_class("52896599943750000") * power[1]
       + mpz_class("2768255397056250000");
  c[1] = -1350 * power[8] - 820050 * power[7] - 94110480 * power[6]
       + mpz_class("31110626544") * power[5]
       + mpz_class("4792288971600") * power[4]
       - mpz_class("732655255111200") * power[3]
       - mpz_class("52193907484200000") * power[2]
       + mpz_class("10710581921943750000") * power[1]
       - mpz_class("415238309558437500000");
  c[0] = 112225 * power[8] + 46738530 * power[7]
       - mpz_class("2973418919") * power[6]
       - mpz_class("2733206451300") * power[5]
       + mpz_class("32192690719900") * power[4]
       + mpz_class("64367011928115000") * power[3]
       - mpz_class("1649065566024697500") * power[2]
       - mpz_class("610446746510853750000") * power[1]
       + mpz_class("34518606881626406250000");
  return c;
}

bool displayed_section(long k) {
  return k == 0 || k == 25 || k == 57 || k == 104 || k == 116 || k == 148;
}

long parse_long(const char* text, const char* name) {
  char* end = nullptr;
  const long answer = std::strtol(text, &end, 10);
  if (end == text || *end != '\0') throw std::invalid_argument(std::string("invalid ") + name);
  return answer;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: search_nagao_integer_offset U_MIN U_MAX K_MIN K_MAX\n";
    return 2;
  }
  const long u_min = parse_long(argv[1], "U_MIN");
  const long u_max = parse_long(argv[2], "U_MAX");
  const long k_min = parse_long(argv[3], "K_MIN");
  const long k_max = parse_long(argv[4], "K_MAX");
  if (u_min <= 0 || u_max < u_min || k_max < k_min) {
    std::cerr << "require 0 < U_MIN <= U_MAX and K_MIN <= K_MAX\n";
    return 2;
  }

  const auto started = std::chrono::steady_clock::now();
  unsigned long long tested = 0;
  unsigned long long hits = 0;
  for (long u = u_min; u <= u_max; ++u) {
    const auto c = coefficients(u);

    // Forward differences make each new k cost only four large additions.
    const mpz_class f0 = evaluate(c, k_min);
    const mpz_class f1 = evaluate(c, k_min + 1);
    const mpz_class f2 = evaluate(c, k_min + 2);
    const mpz_class f3 = evaluate(c, k_min + 3);
    const mpz_class f4 = evaluate(c, k_min + 4);
    mpz_class d0 = f0;
    mpz_class d1 = f1 - f0;
    mpz_class d2 = f2 - 2 * f1 + f0;
    mpz_class d3 = f3 - 3 * f2 + 3 * f1 - f0;
    const mpz_class d4 = f4 - 4 * f3 + 6 * f2 - 4 * f1 + f0;

    for (long k = k_min; k <= k_max; ++k) {
      ++tested;
      if (!displayed_section(k) && d0 >= 0 && mpz_perfect_square_p(d0.get_mpz_t())) {
        mpz_class root;
        mpz_sqrt(root.get_mpz_t(), d0.get_mpz_t());
        std::cout << u << '\t' << k << '\t' << root << '\n';
        ++hits;
      }
      d0 += d1;
      d1 += d2;
      d2 += d3;
      d3 += d4;
    }
  }
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  std::cerr << "tested=" << tested << " hits=" << hits
            << " seconds=" << seconds << '\n';
  return 0;
}
