// Exact bounded q=80 CM24 seed search over a compile-time prime.
//
// Compile and run with:
//   c++ -O3 -std=c++17 -DPRIME=13 elkies-k3/scripts/search_q80_p3_mod11.cpp -o /tmp/search_q80_p3
//   /tmp/search_q80_p3_mod11

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {
#ifndef PRIME
#define PRIME 11
#endif
constexpr int P = PRIME;
using Poly = std::array<int, 32>;

int mod(int value) {
  value %= P;
  return value < 0 ? value + P : value;
}

int power(int base, int exponent) {
  int result = 1;
  for (; exponent; exponent >>= 1, base = mod(base * base))
    if (exponent & 1) result = mod(result * base);
  return result;
}

std::uint64_t integer_power(std::uint64_t base, int exponent) {
  std::uint64_t result = 1;
  while (exponent--) result *= base;
  return result;
}

int inverse(int value) { return power(value, P - 2); }

Poly multiply(const Poly& left, const Poly& right, int max_degree) {
  Poly result{};
  int left_degree = max_degree, right_degree = max_degree;
  while (left_degree >= 0 && left[left_degree] == 0) --left_degree;
  while (right_degree >= 0 && right[right_degree] == 0) --right_degree;
  for (int i = 0; i <= left_degree; ++i)
    for (int j = 0; j <= right_degree && j + i <= max_degree; ++j)
      result[i + j] = mod(result[i + j] + left[i] * right[j]);
  return result;
}

Poly subtract(const Poly& left, const Poly& right) {
  Poly result{};
  for (int index = 0; index < 32; ++index)
    result[index] = mod(left[index]-right[index]);
  return result;
}

Poly negate(const Poly& value) {
  Poly result{};
  for (int index = 0; index < 32; ++index) result[index] = mod(-value[index]);
  return result;
}

int polynomial_degree(const Poly& value) {
  int degree = 31;
  while (degree >= 0 && value[degree] == 0) --degree;
  return degree;
}

Poly remainder(Poly numerator, const Poly& denominator) {
  const int denominator_degree = polynomial_degree(denominator);
  if (denominator_degree < 0) return numerator;
  const int inverse_leading = inverse(denominator[denominator_degree]);
  int numerator_degree = polynomial_degree(numerator);
  while (numerator_degree >= denominator_degree) {
    const int shift = numerator_degree-denominator_degree;
    const int factor = mod(numerator[numerator_degree]*inverse_leading);
    for (int index = 0; index <= denominator_degree; ++index)
      numerator[index+shift] = mod(
          numerator[index+shift]-factor*denominator[index]);
    numerator_degree = polynomial_degree(numerator);
  }
  return numerator;
}

int gcd_degree(Poly left, Poly right) {
  while (polynomial_degree(right) >= 0) {
    Poly next = remainder(left, right);
    left = right;
    right = next;
  }
  return polynomial_degree(left);
}

int evaluate(const Poly& value, int point, int max_degree);
int derivative_evaluate(
    const Poly& value, int point, int max_degree, int order = 1);

Poly polynomial_gcd(Poly left, Poly right) {
  while (polynomial_degree(right) >= 0) {
    Poly next = remainder(left, right);
    left = right;
    right = next;
  }
  const int degree = polynomial_degree(left);
  if (degree < 0) return left;
  const int scale = inverse(left[degree]);
  for (int index = 0; index <= degree; ++index)
    left[index] = mod(left[index]*scale);
  return left;
}

Poly quotient_exact(Poly numerator, const Poly& denominator) {
  Poly quotient{};
  const int denominator_degree = polynomial_degree(denominator);
  if (denominator_degree < 0) return quotient;
  const int inverse_leading = inverse(denominator[denominator_degree]);
  int numerator_degree = polynomial_degree(numerator);
  while (numerator_degree >= denominator_degree) {
    const int shift = numerator_degree-denominator_degree;
    const int factor = mod(numerator[numerator_degree]*inverse_leading);
    quotient[shift] = mod(quotient[shift]+factor);
    for (int index = 0; index <= denominator_degree; ++index)
      numerator[index+shift] = mod(
          numerator[index+shift]-factor*denominator[index]);
    numerator_degree = polynomial_degree(numerator);
  }
  if (numerator_degree >= 0) std::abort();
  return quotient;
}

int pair_intersection(const Poly& XP, const Poly& YP, const Poly& ZP,
                      const Poly& XQ, const Poly& YQ, const Poly& ZQ) {
  const Poly ZP2 = multiply(ZP, ZP, 31);
  const Poly ZQ2 = multiply(ZQ, ZQ, 31);
  const Poly ZP3 = multiply(ZP2, ZP, 31);
  const Poly ZQ3 = multiply(ZQ2, ZQ, 31);
  const Poly D = subtract(multiply(XP, ZQ2, 31), multiply(XQ, ZP2, 31));
  const Poly S = subtract(multiply(YP, ZQ3, 31), negate(multiply(YQ, ZP3, 31)));
  const Poly H = multiply(D, multiply(ZP, ZQ, 31), 31);
  const Poly D2 = multiply(D, D, 31);
  const Poly coordinate_sum = subtract(
      multiply(XP, ZQ2, 31), negate(multiply(XQ, ZP2, 31)));
  const Poly N = subtract(multiply(S, S, 31), multiply(D2, coordinate_sum, 31));
  const Poly first = polynomial_gcd(H, N);
  const Poly N_over_first = quotient_exact(N, first);
  const Poly cancellation = polynomial_gcd(H, N_over_first);
  const Poly H_reduced = quotient_exact(H, cancellation);
  const Poly cancellation2 = multiply(cancellation, cancellation, 31);
  const Poly N_reduced = quotient_exact(N, cancellation2);
  const int finite = polynomial_degree(H_reduced);
  const int excess = polynomial_degree(N_reduced)-2*finite-4;
  if (excess > 0 && (excess & 1)) std::abort();
  return finite+(excess > 0 ? excess/2 : 0);
}

bool opposite_i4_branches(const Poly& A, const Poly& B, int d,
                          const Poly& P1X, const Poly& P1Y,
                          const Poly& X3, const Poly& Y3, const Poly& Z3) {
  const int z_one = evaluate(Z3, 1, 1);
  if (z_one == 0) return false;
  const int xi1 = derivative_evaluate(P1X, 1, 3);
  const int eta1 = derivative_evaluate(P1Y, 1, 5);
  const int x3_one = evaluate(X3, 1, 6);
  const int xi3 = mod(
      (derivative_evaluate(X3, 1, 6)*z_one-2*x3_one)
      * inverse(mod(z_one*z_one*z_one)));
  const int eta3 = mod(
      derivative_evaluate(Y3, 1, 9)
      * inverse(mod(power(z_one, 3))));
  const int a1 = derivative_evaluate(A, 1, 5);
  const int cone_constant = mod(
      (d*derivative_evaluate(A, 1, 5, 2)
       + derivative_evaluate(B, 1, 8, 2))*inverse(2));
  return mod(
      eta1*eta3 + 3*d*xi1*xi3
      + a1*inverse(2)*(xi1+xi3) + cone_constant) == 0;
}

bool is_square(int value) {
  static const std::array<bool, P> table = [] {
    std::array<bool, P> answer{};
    for (int candidate = 0; candidate < P; ++candidate)
      answer[mod(candidate * candidate)] = true;
    return answer;
  }();
  return table[mod(value)];
}

int evaluate(const Poly& value, int point, int max_degree) {
  int result = 0;
  for (int i = max_degree; i >= 0; --i) result = mod(result * point + value[i]);
  return result;
}

int derivative_evaluate(const Poly& value, int point, int max_degree, int order) {
  int result = 0;
  for (int i = max_degree; i >= order; --i) {
    int coefficient = value[i];
    for (int j = 0; j < order; ++j) coefficient = mod(coefficient * (i - j));
    result = mod(result * point + coefficient);
  }
  return result;
}

int binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  int result = 1;
  for (int i = 1; i <= k; ++i)
    result = mod(result * (n + 1 - i) * inverse(i));
  return result;
}

bool solve4(int augmented[4][5], int answer[4]) {
  for (int column = 0; column < 4; ++column) {
    int pivot = column;
    while (pivot < 4 && augmented[pivot][column] == 0) ++pivot;
    if (pivot == 4) return false;
    for (int j = column; j < 5; ++j)
      std::swap(augmented[column][j], augmented[pivot][j]);
    const int scale = inverse(augmented[column][column]);
    for (int j = column; j < 5; ++j) augmented[column][j] = mod(augmented[column][j] * scale);
    for (int row = 0; row < 4; ++row) {
      if (row == column) continue;
      const int factor = augmented[row][column];
      for (int j = column; j < 5; ++j)
        augmented[row][j] = mod(augmented[row][j] - factor * augmented[column][j]);
    }
  }
  for (int i = 0; i < 4; ++i) answer[i] = augmented[i][4];
  return true;
}

bool square_root(const Poly& value, Poly& root) {
  int degree = 31;
  while (degree >= 0 && value[degree] == 0) --degree;
  if (degree < 0) {
    root.fill(0);
    return true;
  }
  if (degree & 1) return false;
  const int half = degree / 2;
  int leading_root = -1;
  for (int candidate = 0; candidate < P; ++candidate)
    if (mod(candidate * candidate) == value[degree]) {
      leading_root = candidate;
      break;
    }
  if (leading_root < 0) return false;
  root.fill(0);
  root[half] = leading_root;
  const int denominator_inverse = inverse(2 * leading_root);
  for (int index = half - 1; index >= 0; --index) {
    const int target = half + index;
    int known = 0;
    for (int left = index + 1; left < half; ++left) {
      const int right = target - left;
      if (index < right && right < half)
        known = mod(known + root[left] * root[right]);
    }
    root[index] = mod((value[target] - known) * denominator_inverse);
  }
  const Poly check = multiply(root, root, 31);
  return check == value;
}

void print_poly(const Poly& value, int degree) {
  std::cout << "[";
  for (int i = 0; i <= degree; ++i) {
    if (i) std::cout << ",";
    std::cout << value[i];
  }
  std::cout << "]";
}
}  // namespace

int search_p3(const Poly& A, const Poly& B, int d, int rho, int node,
              const Poly& P1X, const Poly& P1Y, const Poly& P2X, const Poly& P2Y) {
  int hits = 0;
  for (int pole = 0; pole < P; ++pole) {
#ifdef POLE_START
    if (pole < POLE_START) continue;
#endif
#ifdef POLE_END
    if (pole >= POLE_END) continue;
#endif
    if (pole == 0 || pole == 1 || pole == rho) continue;
    Poly Z{}; Z[0] = mod(-pole); Z[1] = 1;
    const Poly Z2 = multiply(Z, Z, 2);
    const Poly Z4 = multiply(Z2, Z2, 4);
    const Poly Z6 = multiply(Z4, Z2, 6);
    const int matrix_det = mod(power(rho, 6) - power(rho, 5));
    const int matrix_det_inverse = inverse(matrix_det);

    std::vector<int> nonzero_squares;
    for (int value = 1; value < P; ++value)
      if (is_square(value)) nonzero_squares.push_back(value);
    auto complete_top = [&](Poly& X, int leading) {
      int sum_one_without_x4 = 0, sum_rho_without_x4 = 0;
      for (int i = 0; i < 4; ++i) {
        sum_one_without_x4 = mod(sum_one_without_x4 + X[i]);
        sum_rho_without_x4 = mod(
            sum_rho_without_x4 + X[i] * power(rho, i));
      }
      const int right_one_without_x4 = mod(
          d * power(evaluate(Z, 1, 1), 2) - sum_one_without_x4);
      const int right_rho_without_x4 = mod(
          node * power(evaluate(Z, rho, 1), 2) - sum_rho_without_x4);
      const int x6_without_x4 = mod(
          (right_rho_without_x4-power(rho, 5)*right_one_without_x4)
          * matrix_det_inverse);
      X[4] = mod(rho*(leading-x6_without_x4));
      const int right_one = mod(right_one_without_x4-X[4]);
      const int right_rho = mod(right_rho_without_x4-X[4]*power(rho, 4));
      X[5] = mod((right_one * power(rho, 6) - right_rho) * matrix_det_inverse);
      X[6] = mod((right_rho - power(rho, 5) * right_one) * matrix_det_inverse);
    };
    Poly zero_response{}, x3_response{};
    complete_top(zero_response, 0);
    x3_response[3] = 1;
    complete_top(x3_response, 0);
    const int x3_to_pole = mod(
        evaluate(x3_response, pole, 6)-evaluate(zero_response, pole, 6));
    if (x3_to_pole == 0) return 4;

    const std::uint64_t total = integer_power(P, 2);
    for (const int constant : nonzero_squares)
    for (const int leading : nonzero_squares)
    for (const int at_pole : nonzero_squares)
    for (std::uint64_t code = 0; code < total; ++code) {
      std::uint64_t remaining = code;
      Poly X{};
      X[0] = constant;
      for (int i = 1; i < 3; ++i) {
        X[i] = remaining % P;
        remaining /= P;
      }
      complete_top(X, leading);
      const int value_without_x3 = evaluate(X, pole, 6);
      X[3] = mod((at_pole-value_without_x3)*inverse(x3_to_pole));
      complete_top(X, leading);
      if (X[6] != leading) return 3;
      if (evaluate(X, pole, 6) != at_pole) return 5;
      if (!is_square(X[0]) || !is_square(X[6])) continue;

      bool evaluation_square = true;
      for (int point = 2; point <= 8; ++point) {
        const int x_value = evaluate(X, point, 6);
        const int z_value = evaluate(Z, point, 1);
        const int z2_value = mod(z_value*z_value);
        const int z4_value = mod(z2_value*z2_value);
        const int z6_value = mod(z4_value*z2_value);
        const int right = mod(
            mod(x_value*x_value)*x_value
            + evaluate(A, point, 5)*x_value*z4_value
            + evaluate(B, point, 8)*z6_value);
        if (!is_square(right)) {
          evaluation_square = false;
          break;
        }
      }
      if (!evaluation_square) continue;

      const Poly X2 = multiply(X, X, 12);
      const Poly X3 = multiply(X2, X, 18);
      const Poly AX = multiply(A, X, 11);
      const Poly AXZ4 = multiply(AX, Z4, 15);
      const Poly BZ6 = multiply(B, Z6, 14);
      Poly value{};
      for (int i = 0; i <= 18; ++i)
        value[i] = mod(X3[i] + AXZ4[i] + BZ6[i]);
      Poly Y{};
      if (!square_root(value, Y)) continue;
      if (evaluate(Y, 1, 9) != 0 || evaluate(Y, rho, 9) != 0) continue;
      if (evaluate(Y, pole, 9) == 0) continue;
      Poly one{}; one[0] = 1;
      for (int sign3 : {1, -1}) {
        const Poly oriented_Y = sign3 == 1 ? Y : negate(Y);
        if (!opposite_i4_branches(
                A, B, d, P1X, P1Y, X, oriented_Y, Z))
          continue;
        if (pair_intersection(P1X, P1Y, one, X, oriented_Y, Z) != 3)
          continue;
        if (pair_intersection(P2X, P2Y, one, X, oriented_Y, Z) != 1)
          continue;
        ++hits;
        std::cout << "Q80P3MODP|HIT|prime=" << P << "|d=" << d << "|rho=" << rho << "|node=" << node
                  << "|P1X=";
        print_poly(P1X, 3);
        std::cout << "|P1Y=";
        print_poly(P1Y, 5);
        std::cout << "|P2X=";
        print_poly(P2X, 4);
        std::cout << "|P2Y=";
        print_poly(P2Y, 6);
        std::cout << "|pole=" << pole << "|X=";
        print_poly(X, 6);
        std::cout << "|Y=";
        print_poly(oriented_Y, 9);
        std::cout << "\n";
      }
    }
  }
  return hits;
}

int main() {
  std::uint64_t ambient_tested = 0, cm_tested = 0, p1_tested = 0, p2_tested = 0;
  int hits = 0;
  for (int d = 1; d < P; ++d) {
#ifdef FIX_D
    if (d != mod(FIX_D)) continue;
#endif
    for (int p = 0; p < P; ++p) {
#ifdef P_START
      if (p < P_START) continue;
#endif
#ifdef P_END
      if (p >= P_END) continue;
#endif
#ifdef FIX_P
      if (p != mod(FIX_P)) continue;
#endif
      for (int q = 0; q < P; ++q) {
#ifdef P1_FAMILY
        if (q != mod(18 - 2 * p)) continue;
#endif
#ifdef FIX_Q
        if (q != mod(FIX_Q)) continue;
#endif
        for (int e = 0; e < P; ++e) {
#ifdef P1_FAMILY
          if (e != mod(mod(p - 42) * mod(p - 42) * inverse(36))) continue;
#endif
#ifdef FIX_E
          if (e != mod(FIX_E)) continue;
#endif
          ++ambient_tested;
          const int r = mod(-3*d*d + 3 - p - q);
          Poly A{};
          A[2] = mod(-3); A[3] = p; A[4] = q; A[5] = r;

          // Cubic binomial jet of the chosen square-root branch at T=1.
          Poly ajet{};
          for (int jet = 0; jet < 4; ++jet)
            for (int degree = 2; degree <= 5; ++degree)
              ajet[jet] = mod(ajet[jet] + A[degree] * binomial(degree, jet));
          Poly u{};
          const int denominator_inverse = inverse(mod(-3*d*d));
          for (int jet = 1; jet < 4; ++jet) u[jet] = mod(ajet[jet] * denominator_inverse);
          Poly u2 = multiply(u, u, 3), u3 = multiply(u2, u, 3), branch{};
          branch[0] = 1;
          for (int jet = 0; jet < 4; ++jet)
            branch[jet] = mod(branch[jet] + mod(3*inverse(2))*u[jet]
                              + mod(3*inverse(8))*u2[jet] - inverse(16)*u3[jet]);
          for (int jet = 0; jet < 4; ++jet) branch[jet] = mod(2*d*d*d*branch[jet]);
          int augmented[4][5]{};
          for (int jet = 0; jet < 4; ++jet) {
            for (int column = 0; column < 4; ++column)
              augmented[jet][column] = binomial(4 + column, jet);
            const int fixed = mod(2*binomial(3, jet) + e*binomial(8, jet));
            augmented[jet][4] = mod(branch[jet] - fixed);
          }
          int bcoeff[4]{};
          if (!solve4(augmented, bcoeff)) return 2;
          Poly B{};
          B[3] = 2;
          for (int i = 0; i < 4; ++i) B[4+i] = bcoeff[i];
          B[8] = e;

          Poly A2 = multiply(A, A, 10), A3 = multiply(A2, A, 15);
          Poly B2 = multiply(B, B, 16), discriminant{};
          for (int i = 0; i <= 16; ++i)
            discriminant[i] = mod(4*A3[i] + 27*B2[i]);

          for (int rho = 2; rho < P; ++rho) {
            if (evaluate(discriminant, rho, 16) != 0
                || derivative_evaluate(discriminant, rho, 16) != 0
                || derivative_evaluate(discriminant, rho, 16, 2) == 0) continue;
            int node = -1;
            for (int x = 0; x < P; ++x)
              if (mod(x*x*x + evaluate(A, rho, 5)*x + evaluate(B, rho, 8)) == 0
                  && mod(3*x*x + evaluate(A, rho, 5)) == 0) node = x;
            if (node < 0) continue;
            ++cm_tested;

            std::vector<std::pair<Poly, Poly>> p1s;
            for (int a = 0; a < P; ++a) for (int b = 0; b < P; ++b) {
              Poly X{}; X[1] = a; X[2] = b; X[3] = mod(d-a-b);
              Poly X2 = multiply(X, X, 8), X3 = multiply(X2, X, 12);
              Poly AX = multiply(A, X, 9), value{};
              for (int i = 0; i <= 12; ++i) value[i] = mod(X3[i] + AX[i] + B[i]);
              Poly Y{};
              if (!square_root(value, Y) || evaluate(X, rho, 3) != node
                  || evaluate(Y, 1, 5) != 0 || evaluate(Y, rho, 5) != 0) continue;
              // Exact nonidentity charts.  At I1* the spinor branch follows
              // the double root x/T=1 and has (ord_0 x,ord_0 y)=(1,2).
              // (The simple root -2 is the order-two D5 class.)  At IV*, after
              // (xbar,ybar)=(u^4*x(1/u),u^6*y(1/u)), a nonzero E6 class has
              // valuations at least (2,2), hence deg(x)<=2 and deg(y)=4.  The previous
              // X[3],Y[5] nonvanishing test selected the wrong IV* stratum.
              if (X[1] != 1 || Y[2] == 0 || X[3] != 0 || Y[5] != 0
                  || Y[4] == 0)
                continue;
              p1s.push_back({X, Y});
            }
            if (p1s.empty()) continue;
#ifdef REQUIRE_P1_STANDARD
            {
              std::vector<std::pair<Poly, Poly>> filtered;
              for (const auto& P1 : p1s)
                if (P1.first[1] == 1 && P1.first[2] == 2 && P1.first[3] == 0)
                  filtered.push_back(P1);
              p1s.swap(filtered);
              if (p1s.empty()) continue;
            }
#endif
            ++p1_tested;
#ifdef STOP_AFTER_P1
            for (const auto& P1 : p1s) {
              std::cout << "Q80P1MODP|prime=" << P << "|d,p,q,e="
                        << d << "," << p << "," << q << "," << e
                        << "|rho=" << rho << "|node=" << node << "|X=";
              print_poly(P1.first, 3);
              std::cout << "|Y=";
              print_poly(P1.second, 5);
              std::cout << "\n";
            }
            continue;
#endif

            std::vector<std::pair<Poly, Poly>> p2s;
            std::vector<int> nonzero_squares;
            for (int value = 1; value < P; ++value)
              if (is_square(value)) nonzero_squares.push_back(value);
            const std::uint64_t middle_total = integer_power(P, 3);
            for (const int constant : nonzero_squares)
            for (const int leading : nonzero_squares)
            for (std::uint64_t code = 0; code < middle_total; ++code) {
              std::uint64_t remaining = code; Poly X{};
              X[0] = constant;
              for (int i = 1; i < 4; ++i) { X[i] = remaining % P; remaining /= P; }
              X[4] = leading;
              if (evaluate(X, 1, 4) == d) continue;
              // Several fiber evaluations are cheap necessary square tests.
              // They reject almost all quartics before any convolution.
              bool evaluation_square = true;
              for (int point = 1; point <= 6; ++point) {
                const int x_value = evaluate(X, point, 4);
                const int right = mod(
                    mod(x_value * x_value) * x_value
                    + evaluate(A, point, 5) * x_value
                    + evaluate(B, point, 8));
                if (!is_square(right)) {
                  evaluation_square = false;
                  break;
                }
              }
              if (!evaluation_square) continue;
              Poly X2 = multiply(X, X, 8), X3 = multiply(X2, X, 12);
              Poly AX = multiply(A, X, 9), value{};
              for (int i = 0; i <= 12; ++i) value[i] = mod(X3[i] + AX[i] + B[i]);
              Poly Y{};
              if (!square_root(value, Y)) continue;
              // P2 has trivial component label at the residual I2.  Passing
              // through its singular cubic node would instead give label 1;
              // this missing gate created the earlier spurious standard-P1
              // modular branch.
              if (evaluate(X, rho, 4) == node && evaluate(Y, rho, 6) == 0)
                continue;
              p2s.push_back({X, Y});
            }
            if (p2s.empty()) continue;
            std::vector<std::array<Poly, 4>> paired_sections;
            Poly one{}; one[0] = 1;
            for (const auto& P1 : p1s) for (const auto& P2 : p2s) {
              for (int sign2 : {1, -1}) {
                const Poly oriented_Y2 = sign2 == 1 ? P2.second : negate(P2.second);
                if (pair_intersection(
                        P1.first, P1.second, one,
                        P2.first, oriented_Y2, one) == 2)
                  paired_sections.push_back({P1.first, P1.second, P2.first, oriented_Y2});
              }
            }
            if (paired_sections.empty()) continue;
            ++p2_tested;
#ifdef STOP_AFTER_P2
            for (const auto& pair : paired_sections) {
              std::cout << "Q80P2MODP|prime=" << P << "|d,p,q,e="
                        << d << "," << p << "," << q << "," << e
                        << "|rho=" << rho << "|node=" << node << "|X=";
              print_poly(pair[2], 4);
              std::cout << "|Y=";
              print_poly(pair[3], 6);
              std::cout << "\n";
            }
            continue;
#endif
            for (const auto& pair : paired_sections) {
              std::cout << "Q80P3MODP|CANDIDATE|prime=" << P << "|d,p,q,e=" << d << "," << p << ","
                        << q << "," << e << "|rho=" << rho << "|node=" << node << "\n";
              hits += search_p3(A, B, d, rho, node, pair[0], pair[1], pair[2], pair[3]);
            }
          }
        }
      }
    }
  }
  std::cout << "Q80P3MODP|SUMMARY|prime=" << P << "|ambient=" << ambient_tested << "|cm=" << cm_tested
            << "|with_P1=" << p1_tested << "|with_P2=" << p2_tested << "|hits=" << hits
            << "|status=BOUNDED_EXPERIMENT\n";
  return hits ? 0 : 1;
}
