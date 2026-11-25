# Electronic supplement

This supplement contains the exact checks for *The Strong Factorial
Conjecture for Univariate Binomials*. The
[positive-exponent boundary appendix](positive-boundary/README.md)
contains the written contour arguments cited in Appendix C and their
exact input and replay instructions.

## Finite coprimality and squarefreeness

Run these commands from this directory with a C++20 compiler that
supports OpenMP. The programs use integer arithmetic, with Euclidean
division over prime fields in the modular checks.

    check_dir=$(mktemp -d)
    for source in verify_*.cpp; do
      name=$(basename "$source" .cpp)
      g++ -std=c++20 -O3 -Wall -Wextra -fopenmp "$source" -o "$check_dir/$name" || exit
      OMP_NUM_THREADS=8 "$check_dir/$name" || exit
    done

| Program | Role in the proof | Verified cases |
| --- | --- | ---: |
| [verify_initial_terminals.cpp](verify_initial_terminals.cpp) | Initial orders with positive lower exponent | 217,952 |
| [verify_ordinary_binomial_residual.cpp](verify_ordinary_binomial_residual.cpp) | Central range after the prime-support and Padé reductions | 74,362 |
| [verify_small_support.cpp](verify_small_support.cpp) | Small positive lower exponents and adjacent supports below the Padé cutoff | 18,904 |
| [verify_constant_support.cpp](verify_constant_support.cpp) | Constant-plus-power squarefreeness in Appendix B | 202,021 |
| [verify_positive_boundary.cpp](verify_positive_boundary.cpp) | The three finite domains at large relative gaps in Appendix C | 326,146 |

Every program reports zero unresolved cases. The manuscript
specifies the enumeration bounds and proves that the applicable filters
cover every possibility.

The central verifier's table stores one 16-bit entry per odd argument
and covers arguments through 1,156,520,634. A replay with g++ 13.3.0,
eight threads and undefined-behavior instrumentation took 11.63 seconds
and used 1,133,864 KiB of peak resident memory on Ubuntu 24.04.
The decimal square-root estimate used by the sieve is corrected in
both directions by exact integer squares.

An additional replay checks the 64 polynomial pairs left after the
arithmetic sieve in a different prime field, using FLINT's gcd.
Install python-flint and run this command after the C++ build above.

    set -o pipefail
    "$check_dir/verify_ordinary_binomial_residual" --list-gcd-pairs | python3 crosscheck_central_gcds.py

This check uses the pairs exported by the C++ sieve. The argument in
Section 3 of the article justifies the coverage of that sieve.

## Padé inequalities

Run the three scripts with Python 3 and its standard library.

    python3 verify_pade_intervals.py
    python3 verify_hard_pade_endpoints.py
    python3 verify_large_tail.py

The first script checks the prime sums, the directed angular bounds,
and scalar inequalities for the intermediate orders. The second checks
the two-height comparisons at the endpoints of the terminal and gap
intervals. The third checks the endpoint inequalities for
$N\ge e^{15}$. The scripts enclose logarithms, exponentials, and
trigonometric values by rational or directed interval arithmetic.
The article proves how these finite comparisons extend across the
stated parameter ranges.
