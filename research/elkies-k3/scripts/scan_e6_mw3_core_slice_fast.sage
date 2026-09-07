from sage.all import *
from pathlib import Path
import argparse
import subprocess
import tempfile


ap = argparse.ArgumentParser(
    description=(
        "Scan E6/P1 (r0,s0,x1) slices in compiled finite-field arithmetic and "
        "emit only cores with exact I4,I4 and two rational exact-I2 fibers."
    )
)
ap.add_argument("--p", type=int, default=31)
ap.add_argument("--r0", type=int)
ap.add_argument("--s0", type=int)
ap.add_argument("--x1", type=int)
ap.add_argument("--all-slices", action="store_true")
ap.add_argument("--slice-start", type=int, default=0)
ap.add_argument("--slice-count", type=int)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--timeout", type=float, default=30.0)
ap.add_argument("--meta", default="artifacts/local/elkies-k3/e6-base.meta.txt")
args = ap.parse_args()

p = args.p
if not is_prime(p) or p in (2, 3, 79):
    raise SystemExit("choose a good prime other than 2,3,79")
if args.all_slices:
    if any(value is not None for value in (args.r0, args.s0, args.x1)):
        raise SystemExit("--all-slices cannot be combined with --r0,--s0,--x1")
    slices = [
        (r0, s0, x1)
        for r0 in range(1, p)
        for s0 in range(1, p)
        for x1 in range(p)
    ]
    if args.slice_start < 0 or args.slice_start >= len(slices):
        raise SystemExit("--slice-start is outside the all-slices range")
    stop = len(slices) if args.slice_count is None else args.slice_start+args.slice_count
    if args.slice_count is not None and args.slice_count <= 0:
        raise SystemExit("--slice-count must be positive")
    slices = slices[args.slice_start:min(stop, len(slices))]
else:
    if args.slice_start != 0 or args.slice_count is not None:
        raise SystemExit("--slice-start/--slice-count require --all-slices")
    if any(value is None for value in (args.r0, args.s0, args.x1)):
        raise SystemExit("give --r0,--s0,--x1, or use --all-slices")
    fixed = tuple(value % p for value in (args.r0, args.s0, args.x1))
    if fixed[0] == 0 or fixed[1] == 0:
        raise SystemExit("r0 and s0 must be nonzero in this triangular chart")
    slices = [fixed]

names = [
    "a1", "a2", "a3", "a4", "b4", "b5", "lam", "mu", "s0", "s1",
    "sl", "sm", "x0", "x1", "y0", "y1", "y2",
]
R = PolynomialRing(GF(p), names)
F = FractionField(R)

metadata = {}
for line in Path(args.meta).read_text().splitlines():
    if "=" in line:
        key, value = line.split("=", 1)
        metadata[key.strip()] = value.strip()
expressions = {
    key: F(metadata[key].replace("^", "**"))
    for key in ("a3", "y2", "b5", "b4")
}


def cpp_polynomial(function_name, polynomial):
    lines = [
        f"static inline int {function_name}(const int *v) {{",
        "  int acc = 0, term = 0;",
    ]
    for exponents, coefficient in R(polynomial).dict().items():
        lines.append(f"  term = {int(coefficient)};")
        for index, exponent in enumerate(exponents):
            if exponent:
                lines.append(f"  term = md((long long)term * pm(v[{index}], {exponent}));")
        lines.append("  acc = md(acc + term);")
    lines.append("  return acc;")
    lines.append("}")
    return "\n".join(lines)


pieces = []
for key, expression in expressions.items():
    pieces.append(cpp_polynomial(f"{key}_num", expression.numerator()))
    pieces.append(cpp_polynomial(f"{key}_den", expression.denominator()))
formula_functions = "\n".join(pieces)
indices = {name: index for index, name in enumerate(names)}
slice_initializer = ",\n".join(
    "  {" + ",".join(map(str, row)) + "}" for row in slices
)

source = f"""
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

static constexpr int P = {p};
static constexpr int SLICE_COUNT = {len(slices)};
static const int SLICES[SLICE_COUNT][3] = {{
{slice_initializer}
}};
static int inverse[P];
static int choose_table[17][17];

static inline int md(long long value) {{
  value %= P;
  return value < 0 ? int(value + P) : int(value);
}}

static inline int pm(int value, int exponent) {{
  int result = 1;
  while (exponent) {{
    if (exponent & 1) result = md((long long)result * value);
    value = md((long long)value * value);
    exponent >>= 1;
  }}
  return result;
}}

{formula_functions}

static inline bool fraction_value(
    int (*numerator)(const int *), int (*denominator)(const int *),
    const int *v, int &result) {{
  int den = denominator(v);
  if (den == 0) return false;
  result = md((long long)numerator(v) * inverse[den]);
  return true;
}}

static inline void convolution(
    const int *left, int left_degree, const int *right, int right_degree,
    int *output) {{
  std::fill(output, output + left_degree + right_degree + 1, 0);
  for (int i = 0; i <= left_degree; ++i)
    for (int j = 0; j <= right_degree; ++j)
      output[i+j] = md(output[i+j] + (long long)left[i]*right[j]);
}}

static inline int hasse(const int *polynomial, int degree, int order, int point) {{
  int result = 0;
  for (int index = order; index <= degree; ++index)
    result = md(result + (long long)polynomial[index]
        * choose_table[index][order] * pm(point, index-order));
  return result;
}}

static std::mutex output_mutex;
static std::atomic<unsigned long long> tested{{0}}, chart_points{{0}}, hits{{0}};

static void worker(int worker_id, int worker_count) {{
  int v[17] = {{0}};
  for (int task = worker_id; task < SLICE_COUNT*P; task += worker_count) {{
    int slice_index = task/P;
    int a1 = task%P;
    const int R0 = SLICES[slice_index][0];
    const int S0 = SLICES[slice_index][1];
    const int X1 = SLICES[slice_index][2];
    const int x0 = md((long long)R0*R0 - 2*S0);
    const int y0 = md((long long)R0*(md((long long)R0*R0 - 3*S0)));
    v[{indices['s0']}] = S0;
    v[{indices['x0']}] = x0;
    v[{indices['x1']}] = X1;
    v[{indices['y0']}] = y0;
    v[{indices['a1']}] = a1;
    v[{indices['y1']}] = md((long long)(a1
        + 3*md((long long)R0*R0-S0)*X1) * inverse[md(2*R0)]);
    for (int a2 = 0; a2 < P; ++a2) {{
      v[{indices['a2']}] = a2;
      for (int a4 = 0; a4 < P; ++a4) {{
        v[{indices['a4']}] = a4;
        for (int s1 = 1; s1 < P; ++s1) {{
          ++tested;
          v[{indices['s1']}] = s1;
          int a3, y2, b5, b4;
          if (!fraction_value(a3_num, a3_den, v, a3)) continue;
          v[{indices['a3']}] = a3;
          if (!fraction_value(y2_num, y2_den, v, y2)) continue;
          v[{indices['y2']}] = y2;
          if (!fraction_value(b5_num, b5_den, v, b5)) continue;
          v[{indices['b5']}] = b5;
          if (!fraction_value(b4_num, b4_den, v, b4)) continue;
          v[{indices['b4']}] = b4;
          ++chart_points;

          int a[6], b[9];
          a[0] = md(-3LL*S0*S0);
          a[1] = a1; a[2] = a2; a[3] = a3; a[4] = a4;
          a[5] = md(-3LL*s1*s1-a[0]-a1-a2-a3-a4);
          b[0] = md(2LL*S0*S0*S0);
          b[1] = md(-(long long)S0*a1);
          b[2] = md((long long)a1*a1*inverse[md(12*S0)]-(long long)S0*a2);
          b[3] = md(
              (long long)md((long long)a1*a1*a1
                + 36LL*a1*a2*S0*S0
                - 216LL*a3*S0*S0*S0*S0)
              * inverse[md(216LL*S0*S0*S0)]);
          b[4] = b4; b[5] = b5;
          int c_at_one = md(2LL*s1*s1*s1
              - b[0]-b[1]-b[2]-b[3]-b[4]-b[5]-1);
          int weighted_b = md(b[1]+2LL*b[2]+3LL*b[3]+4LL*b[4]+5LL*b[5]);
          int derivative_a = md(a[1]+2LL*a[2]+3LL*a[3]+4LL*a[4]+5LL*a[5]);
          b[6] = md(weighted_b+7LL*c_at_one+8LL+(long long)s1*derivative_a);
          b[7] = md(c_at_one-b[6]);
          b[8] = 1;

          int a2_poly[11], a3_poly[16], b2_poly[17], delta[17];
          convolution(a, 5, a, 5, a2_poly);
          convolution(a2_poly, 10, a, 5, a3_poly);
          convolution(b, 8, b, 8, b2_poly);
          for (int index = 0; index <= 16; ++index) {{
            int av = index <= 15 ? a3_poly[index] : 0;
            delta[index] = md(4LL*av + 27LL*b2_poly[index]);
          }}
          if (delta[16] == 0) continue;
          bool exact_i4_zero = true;
          for (int order = 0; order < 4; ++order)
            if (delta[order] != 0) exact_i4_zero = false;
          if (!exact_i4_zero || delta[4] == 0) continue;
          bool exact_i4_one = true;
          for (int order = 0; order < 4; ++order)
            if (hasse(delta, 16, order, 1) != 0) exact_i4_one = false;
          if (!exact_i4_one || hasse(delta, 16, 4, 1) == 0) continue;

          int roots[P], root_count = 0;
          bool higher_repeated = false;
          for (int point = 2; point < P; ++point) {{
            int h0 = hasse(delta, 16, 0, point);
            int h1 = hasse(delta, 16, 1, point);
            if (h0 != 0 || h1 != 0) continue;
            if (hasse(delta, 16, 2, point) == 0) {{
              higher_repeated = true;
              break;
            }}
            // ord(Delta)=2 with A=0 is additive type II, not multiplicative I2.
            if (hasse(a, 5, 0, point) == 0) {{
              higher_repeated = true;
              break;
            }}
            roots[root_count++] = point;
          }}
          if (higher_repeated || root_count != 2) continue;
          ++hits;
          std::lock_guard<std::mutex> guard(output_mutex);
          std::cout << "E6FASTCORE|p=" << P << "|r0=" << R0
                    << "|s0=" << S0 << "|x1=" << X1
                    << "|a1=" << a1 << "|a2=" << a2 << "|a4=" << a4
                    << "|s1=" << s1 << "|roots=" << roots[0] << "," << roots[1]
                    << "|A=";
          for (int i = 0; i < 6; ++i) std::cout << (i ? "," : "") << a[i];
          std::cout << "|B=";
          for (int i = 0; i < 9; ++i) std::cout << (i ? "," : "") << b[i];
          std::cout << "|X1=" << x0 << "," << X1 << "," << md(s1-x0-X1)
                    << "|Y1=" << y0 << "," << v[{indices['y1']}] << "," << y2
                    << "," << md(-1-y0-v[{indices['y1']}]-y2) << ",1"
                    << std::endl;
        }}
      }}
    }}
  }}
}}

int main(int argc, char **argv) {{
  int threads = argc > 1 ? std::max(1, std::atoi(argv[1])) : 8;
  for (int value = 1; value < P; ++value)
    for (int candidate = 1; candidate < P; ++candidate)
      if (md((long long)value*candidate) == 1) {{ inverse[value] = candidate; break; }}
  for (int n = 0; n <= 16; ++n) {{
    choose_table[n][0] = choose_table[n][n] = 1;
    for (int k = 1; k < n; ++k)
      choose_table[n][k] = md(choose_table[n-1][k-1]+choose_table[n-1][k]);
  }}
  std::vector<std::thread> pool;
  for (int index = 0; index < threads; ++index)
    pool.emplace_back(worker, index, threads);
  for (auto &thread : pool) thread.join();
  std::cout << "E6FASTSUMMARY|p=" << P << "|slices=" << SLICE_COUNT
            << "|tested=" << tested.load() << "|chart_points=" << chart_points.load()
            << "|hits=" << hits.load() << "|threads=" << threads << std::endl;
}}
"""

with tempfile.TemporaryDirectory(prefix="e6-mw3-core-") as temporary:
    binary = Path(temporary)/"scan"
    compilation = subprocess.run(
        ["c++", "-x", "c++", "-", "-O3", "-std=c++17", "-pthread", "-o", str(binary)],
        input=source, text=True, capture_output=True, check=False,
    )
    if compilation.returncode:
        raise RuntimeError("C++ compilation failed:\n"+compilation.stderr)
    try:
        result = subprocess.run(
            [str(binary), str(args.threads)], text=True, capture_output=True,
            timeout=args.timeout, check=False,
        )
    except subprocess.TimeoutExpired as error:
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode(errors="replace")
        print(partial, end="")
        print(f"E6FAST|status=TIMEOUT|seconds={args.timeout}", flush=True)
        raise SystemExit(124)
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    if result.returncode:
        raise SystemExit(result.returncode)
