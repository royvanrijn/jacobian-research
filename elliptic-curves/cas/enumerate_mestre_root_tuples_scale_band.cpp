// Exact band enumerator for primitive affine-normalized Mestre root tuples.
// It intentionally leaves the frozen max-200 source untouched.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>

using i128 = __int128_t;

static bool normalized(const std::array<int, 6>& roots) {
  std::array<int, 6> reflected{};
  for (int i = 0; i < 6; ++i) reflected[i] = roots[5] - roots[5 - i];
  return roots <= reflected;
}

static bool symmetric(const std::array<int, 6>& roots) {
  return roots[0] + roots[5] == roots[1] + roots[4] &&
         roots[0] + roots[5] == roots[2] + roots[3];
}

static i128 obstruction(const std::array<int, 6>& roots) {
  std::array<i128, 6> elementary{};
  elementary[0] = 1;
  for (int root : roots)
    for (int degree = 5; degree >= 1; --degree)
      elementary[degree] += static_cast<i128>(root) * elementary[degree - 1];
  const i128 s1 = elementary[1], s2 = elementary[2], s3 = elementary[3];
  const i128 s4 = elementary[4], s5 = elementary[5];
  return -s1*s1*s1*s1*s1 + 6*s1*s1*s1*s2 - 7*s1*s1*s3 - 8*s1*s2*s2
         + 8*s1*s4 + 12*s2*s3 - 24*s5;
}

int main(int argc, char** argv) {
  if (argc != 3) return 2;
  char* end = nullptr;
  const long first = std::strtol(argv[1], &end, 10);
  if (*argv[1] == '\0' || *end != '\0') return 2;
  const long last = std::strtol(argv[2], &end, 10);
  if (*argv[2] == '\0' || *end != '\0' || first < 5 || first > last || last > 400) return 2;
  std::uint64_t normalized_count = 0, obstruction_count = 0, symmetric_count = 0;
  std::cout << "MESTRE_ROOT_TUPLES_BAND_V1\n";
  for (int e = static_cast<int>(first); e <= last; ++e)
    for (int a = 1; a < e; ++a)
      for (int b = a + 1; b < e; ++b)
        for (int c = b + 1; c < e; ++c)
          for (int d = c + 1; d < e; ++d) {
            const std::array<int, 6> roots{0, a, b, c, d, e};
            int common = 0;
            for (int i = 1; i < 6; ++i) common = std::gcd(common, roots[i]);
            if (common != 1 || !normalized(roots)) continue;
            ++normalized_count;
            if (obstruction(roots) != 0) continue;
            ++obstruction_count;
            const bool reflected = symmetric(roots);
            symmetric_count += reflected;
            std::cout << "R";
            for (int root : roots) std::cout << ' ' << root;
            std::cout << ' ' << (reflected ? 1 : 0) << '\n';
          }
  std::cout << "S " << first << ' ' << last << ' ' << normalized_count << ' '
            << obstruction_count << ' ' << symmetric_count << '\n';
}
