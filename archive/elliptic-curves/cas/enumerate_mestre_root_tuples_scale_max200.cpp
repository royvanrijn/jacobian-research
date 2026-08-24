// Exact bounded enumerator for primitive affine-normalized Mestre root tuples.
//
// This is a standalone continuation of enumerate_mestre_root_tuples_scale.cpp.
// The older source is a frozen max-root-100 input and is deliberately not
// modified.  Here the six roots are 0<a<b<c<d<e<=200.  We quotient integral
// scaling by gcd(a,b,c,d,e)=1 and reflection in e/2 by retaining the
// lexicographically smaller representative.  The emitted obstruction is
// Mestre's exact homogeneous degree-five obstruction
//
//   -s1^5 + 6*s1^3*s2 - 7*s1^2*s3 - 8*s1*s2^2
//       + 8*s1*s4 + 12*s2*s3 - 24*s5.
//
// Signed 128-bit arithmetic is exact throughout this declared range.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>

using i128 = __int128_t;

static bool lexicographically_normalized(const std::array<int, 6>& roots) {
  const int diameter = roots[5];
  std::array<int, 6> reflected{};
  for (int index = 0; index < 6; ++index) {
    reflected[index] = diameter - roots[5 - index];
  }
  return roots <= reflected;
}

static bool reflection_symmetric(const std::array<int, 6>& roots) {
  const int diameter = roots[5];
  for (int index = 0; index < 3; ++index) {
    if (roots[index] + roots[5 - index] != diameter) return false;
  }
  return true;
}

static i128 obstruction(const std::array<int, 6>& roots) {
  std::array<i128, 6> symmetric{};
  symmetric[0] = 1;
  for (const int root : roots) {
    for (int degree = 5; degree >= 1; --degree) {
      symmetric[degree] += static_cast<i128>(root) * symmetric[degree - 1];
    }
  }
  const i128 s1 = symmetric[1];
  const i128 s2 = symmetric[2];
  const i128 s3 = symmetric[3];
  const i128 s4 = symmetric[4];
  const i128 s5 = symmetric[5];
  return -s1 * s1 * s1 * s1 * s1 + 6 * s1 * s1 * s1 * s2
         - 7 * s1 * s1 * s3 - 8 * s1 * s2 * s2 + 8 * s1 * s4
         + 12 * s2 * s3 - 24 * s5;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: enumerate_mestre_root_tuples_scale_max200 MAX_ROOT\n";
    return 2;
  }
  char* end = nullptr;
  const long parsed = std::strtol(argv[1], &end, 10);
  if (*argv[1] == '\0' || *end != '\0' || parsed < 5 || parsed > 200) {
    std::cerr << "MAX_ROOT must be an integer in [5,200]\n";
    return 2;
  }
  const int max_root = static_cast<int>(parsed);
  std::uint64_t normalized_count = 0;
  std::uint64_t obstruction_count = 0;
  std::uint64_t reflection_count = 0;
  std::uint64_t nonreflection_count = 0;

  std::cout << "MESTRE_ROOT_TUPLES_V2\n";
  for (int e = 5; e <= max_root; ++e) {
    for (int a = 1; a < e; ++a) {
      for (int b = a + 1; b < e; ++b) {
        for (int c = b + 1; c < e; ++c) {
          for (int d = c + 1; d < e; ++d) {
            const std::array<int, 6> roots{0, a, b, c, d, e};
            int common = 0;
            for (int index = 1; index < 6; ++index) {
              common = std::gcd(common, roots[index]);
            }
            if (common != 1 || !lexicographically_normalized(roots)) continue;
            ++normalized_count;
            if (obstruction(roots) != 0) continue;
            ++obstruction_count;
            const bool symmetric = reflection_symmetric(roots);
            if (symmetric) {
              ++reflection_count;
            } else {
              ++nonreflection_count;
            }
            std::cout << "R";
            for (const int root : roots) std::cout << ' ' << root;
            std::cout << ' ' << (symmetric ? 1 : 0) << '\n';
          }
        }
      }
    }
  }
  std::cout << "S " << max_root << ' ' << normalized_count << ' '
            << obstruction_count << ' ' << reflection_count << ' '
            << nonreflection_count << '\n';
  return 0;
}
