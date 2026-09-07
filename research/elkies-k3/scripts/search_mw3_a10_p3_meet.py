#!/usr/bin/env python3
"""Exhaust the fixed-surface A10/MW3 P3 system by two-sided recursion.

The input is the msolve-format coefficient system produced by
``build_mw3_a10_p3_fixed.sage``.  Its equations are S_2,...,S_16 for

    Y^2 - X^3 - A*X*(t-r)^4 - B*(t-r)^6 = 0.

Rather than expand six rational substitutions or scan all twelve variables,
the search meets the polynomial-square recursion from both ends:

* S_2,S_3,S_4 determine v0,v1,v2 from q0,q1,q2;
* S_16,S_15,S_14 determine v6,v5,v4 from q3,q2,q1;
* S_13 determines v3 and S_5 is the first join constraint.

All arithmetic is over GF(31).  A small C++ kernel is generated in memory,
compiled through stdin, and run with a hard timeout.  The kernel also checks
the two open conditions omitted from the coefficient ideal: the numerator
must not cancel at the proposed pole, and the section must not meet the node
at the nominal identity-component fiber t=1.  Every reported hit must still
be replayed in Sage.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from sage.all import GF, PolynomialRing


def load_system(path: Path):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    names = [name.strip() for name in lines[0].split(",") if name.strip()]
    prime = int(lines[1])
    ring = PolynomialRing(GF(prime), names, order="degrevlex")
    equations = [
        ring((line[:-1] if line.endswith(",") else line).replace("^", "**"))
        for line in lines[2:]
    ]
    if names != ["r", "q0", "q1", "q2", "q3"] + [f"v{i}" for i in range(7)]:
        raise RuntimeError(f"unexpected variable order: {names}")
    if prime != 31 or len(equations) != 15:
        raise RuntimeError("expected the 15 equations S_2,...,S_16 over GF(31)")
    return names, equations


def cpp_equation(index, polynomial):
    terms = []
    for exponents, coefficient in polynomial.dict().items():
        factors = [f"(long long){int(coefficient)}"]
        for variable, exponent in enumerate(exponents):
            if exponent:
                factors.append(f"pw[v[{variable}]][{exponent}]")
        terms.append("*".join(factors))
    expression = "+".join(terms) if terms else "0"
    return (
        f"static inline int e{index}(const int *v) {{"
        f" return md((long long)({expression})); }}"
    )


def build_cpp(equations, lam, nodes, sinf):
    s0, s1, sl = nodes
    functions = "\n".join(cpp_equation(i, equation) for i, equation in enumerate(equations))
    switches = "\n".join(f"case {i}: return e{i}(v);" for i in range(len(equations)))
    return f"""
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

static constexpr int P = 31;
static constexpr int LAM = {lam};
static constexpr int S0 = {s0};
static constexpr int S1 = {s1};
static constexpr int SL = {sl};
static constexpr int SINF = {sinf};
static int pw[P][7];
static int invv[P];

static inline int md(long long value) {{
    value %= P;
    return value < 0 ? int(value + P) : int(value);
}}

{functions}

static inline int equation(int index, const int *v) {{
    switch (index) {{
        {switches}
        default: return -1;
    }}
}}

static inline int solve_linear(int equation_index, int variable, int *v) {{
    int old = v[variable];
    v[variable] = 0;
    int constant = equation(equation_index, v);
    v[variable] = 1;
    int coefficient = md(equation(equation_index, v) - constant);
    v[variable] = old;
    if (coefficient == 0) return constant == 0 ? -2 : -1;
    return md(-(long long)constant * invv[coefficient]);
}}

static inline bool exact_open_conditions(const int *v) {{
    int r = v[0];
    int inverse_lambda = invv[LAM];
    int lambda_minus_r = md(LAM - r);
    int c_at_r = md(
        (long long)S0 * pw[r][2] * md(r - LAM) * md(-inverse_lambda)
        + (long long)SL * pw[lambda_minus_r][2] * r * inverse_lambda
    );
    int q_at_r = md(
        v[1] + (long long)v[2]*r + (long long)v[3]*pw[r][2]
        + (long long)v[4]*pw[r][3] + (long long)SINF*pw[r][4]
    );
    int x_at_r = md(c_at_r + (long long)r*md(r - LAM)*q_at_r);
    if (x_at_r == 0) return false;

    int c_at_one = md(
        (long long)S0 * pw[r][2] * md(1 - LAM) * md(-inverse_lambda)
        + (long long)SL * pw[lambda_minus_r][2] * inverse_lambda
    );
    int q_at_one = md(v[1] + v[2] + v[3] + v[4] + SINF);
    int x_at_one = md(c_at_one + (long long)md(1 - LAM)*q_at_one);
    int y_at_one = md((long long)md(1 - LAM)
        * (v[5] + v[6] + v[7] + v[8] + v[9] + v[10] + v[11]));
    int node_x_at_one = md((long long)S1 * pw[md(1-r)][2]);
    return x_at_one != node_x_at_one || y_at_one != 0;
}}

struct Bottom {{ int q0, v0, v1, v2; }};
struct Top {{ int q3, v4, v5, v6; }};

static std::mutex output_mutex;
static std::atomic<unsigned long long> joins{{0}}, middle{{0}}, raw_hits{{0}},
    open_rejections{{0}}, hits{{0}};

static void worker(int worker_id, int worker_count) {{
    int v[12] = {{0}};
    for (int r = worker_id; r < P; r += worker_count) {{
        if (r == 0 || r == 1 || r == LAM) continue;
        v[0] = r;
        std::array<std::vector<Bottom>, P*P> bottom;
        std::array<std::vector<Top>, P*P> top;

        // Low recursion: S2, S3, S4.
        for (int q0 = 0; q0 < P; ++q0) {{
            v[1] = q0;
            for (int v0 = 0; v0 < P; ++v0) {{
                v[5] = v0;
                if (e0(v) != 0) continue;
                for (int q1 = 0; q1 < P; ++q1) {{
                    v[2] = q1;
                    int v1 = solve_linear(1, 6, v);
                    if (v1 == -1) continue;
                    int v1_begin = v1 == -2 ? 0 : v1;
                    int v1_end = v1 == -2 ? P : v1 + 1;
                    for (int vv1 = v1_begin; vv1 < v1_end; ++vv1) {{
                        v[6] = vv1;
                        for (int q2 = 0; q2 < P; ++q2) {{
                            v[3] = q2;
                            int v2 = solve_linear(2, 7, v);
                            if (v2 == -1) continue;
                            int v2_begin = v2 == -2 ? 0 : v2;
                            int v2_end = v2 == -2 ? P : v2 + 1;
                            for (int vv2 = v2_begin; vv2 < v2_end; ++vv2) {{
                                v[7] = vv2;
                                bottom[q1*P + q2].push_back({{q0, v0, vv1, vv2}});
                            }}
                        }}
                    }}
                }}
            }}
        }}

        // High recursion: S16, S15, S14.  Exact component 10~1 has v6 != 0.
        for (int q3 = 0; q3 < P; ++q3) {{
            v[4] = q3;
            for (int v6 = 1; v6 < P; ++v6) {{
                v[11] = v6;
                if (e14(v) != 0) continue;
                for (int q2 = 0; q2 < P; ++q2) {{
                    v[3] = q2;
                    int v5 = solve_linear(13, 10, v);
                    if (v5 < 0) continue;
                    v[10] = v5;
                    for (int q1 = 0; q1 < P; ++q1) {{
                        v[2] = q1;
                        int v4 = solve_linear(12, 9, v);
                        if (v4 < 0) continue;
                        top[q1*P + q2].push_back({{q3, v4, v5, v6}});
                    }}
                }}
            }}
        }}

        // Join on q1,q2.  S13 determines v3; S5 is the first constraint.
        for (int q1 = 0; q1 < P; ++q1) {{
            v[2] = q1;
            for (int q2 = 0; q2 < P; ++q2) {{
                v[3] = q2;
                int key = q1*P + q2;
                for (const auto &lo : bottom[key]) {{
                    v[1] = lo.q0; v[5] = lo.v0; v[6] = lo.v1; v[7] = lo.v2;
                    for (const auto &hi : top[key]) {{
                        ++joins;
                        v[4] = hi.q3; v[9] = hi.v4; v[10] = hi.v5; v[11] = hi.v6;
                        int v3 = solve_linear(11, 8, v);
                        if (v3 < 0) continue;
                        v[8] = v3;
                        if (e3(v) != 0) continue;
                        ++middle;
                        bool good = true;
                        for (int equation_index = 4; equation_index <= 10; ++equation_index) {{
                            if (equation(equation_index, v) != 0) {{ good = false; break; }}
                        }}
                        if (!good) continue;
                        ++raw_hits;
                        if (!exact_open_conditions(v)) {{ ++open_rejections; continue; }}
                        ++hits;
                        std::lock_guard<std::mutex> guard(output_mutex);
                        std::cout << "MW3A10P3_HIT";
                        const char *names[12] = {{"r","q0","q1","q2","q3","v0","v1","v2","v3","v4","v5","v6"}};
                        for (int i = 0; i < 12; ++i) std::cout << "|" << names[i] << "=" << v[i];
                        std::cout << std::endl;
                    }}
                }}
            }}
        }}
    }}
}}

int main(int argc, char **argv) {{
    int threads = argc > 1 ? std::max(1, std::atoi(argv[1])) : 8;
    for (int value = 0; value < P; ++value) {{
        pw[value][0] = 1;
        for (int exponent = 1; exponent <= 6; ++exponent)
            pw[value][exponent] = md((long long)pw[value][exponent-1] * value);
    }}
    for (int value = 1; value < P; ++value)
        for (int candidate = 1; candidate < P; ++candidate)
            if (md((long long)value*candidate) == 1) {{ invv[value] = candidate; break; }}

    std::vector<std::thread> pool;
    for (int i = 0; i < threads; ++i) pool.emplace_back(worker, i, threads);
    for (auto &thread : pool) thread.join();
    std::cout << "MW3A10P3_SUMMARY|joins=" << joins.load()
              << "|middle=" << middle.load() << "|raw_hits=" << raw_hits.load()
              << "|open_rejections=" << open_rejections.load()
              << "|hits=" << hits.load()
              << "|threads=" << threads << std::endl;
    return 0;
}}
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--lambda", dest="lam", type=int, default=23)
    parser.add_argument("--nodes", default="3,21,10", help="nodes at 0,1,lambda")
    parser.add_argument("--sinf", type=int, default=29)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()

    _, equations = load_system(args.input)
    nodes = tuple(int(value) % 31 for value in args.nodes.split(","))
    if len(nodes) != 3:
        raise SystemExit("--nodes must contain exactly three comma-separated values")
    source = build_cpp(equations, args.lam % 31, nodes, args.sinf % 31)
    with tempfile.TemporaryDirectory(prefix="mw3-a10-p3-") as temporary:
        binary = Path(temporary) / "search"
        compile_result = subprocess.run(
            ["c++", "-x", "c++", "-", "-O3", "-std=c++17", "-pthread", "-o", str(binary)],
            input=source,
            text=True,
            capture_output=True,
            check=False,
        )
        if compile_result.returncode:
            raise RuntimeError("C++ compilation failed:\n" + compile_result.stderr)
        try:
            result = subprocess.run(
                [str(binary), str(args.threads)],
                text=True,
                capture_output=True,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            partial = error.stdout or ""
            if isinstance(partial, bytes):
                partial = partial.decode(errors="replace")
            print(partial, end="")
            print(f"MW3A10P3|status=TIMEOUT|seconds={args.timeout}", flush=True)
            raise SystemExit(124)
        print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="")
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
