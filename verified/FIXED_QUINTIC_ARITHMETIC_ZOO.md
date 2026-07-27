# The arithmetic zoo of one explicit quintic Keller map

The canonical theorem and proofs are consolidated in
[`FIXED_QUINTIC_MODULI_DOMINANCE.md`](../FIXED_QUINTIC_MODULI_DOMINANCE.md).
This companion is the finite certificate ledger.

For the single determinant-`-2` degree-five Keller map with inverse pencil

\[
E_{\Pi,B,C}(S)=\Pi^5S^5-5\Pi S^3-2BS^2+4S-2C,
\]

use the centered monic presentation

\[
\widetilde E_{\Pi,B,C}(T)
=\Pi^5E_{\Pi,B,C}(\Pi^{-2}T),
\]

and the following ledger.  A prime annotation such as \(11:(5)\) records
the squarefree factor-degree partition modulo that prime.  The real column
is the number of real roots.

<!-- BEGIN GENERATED FIXED QUINTIC LEDGER -->
| purpose | target `(Pi,B,C)` | `H_proj` | type | witness primes / exact certificate | real | local certificate |
|---|---|---:|---|---|---:|---|
| split | `(1,0,0)` | 1 | `Q^5` | T(T-1)(T+1)(T-2)(T+2) | 5 | — |
| signature | `(1,-1,-1)` | 1 | `S_5` | 11:(5); 7:(4,1); 3:(3,2) | 1 | — |
| signature | `(-1,-1,-1)` | 1 | `S_5` | 5:(5); 43:(2,1,1,1); nonsquare discriminant | 3 | — |
| signature | `(1,0,-1/2)` | 2 | `S_5` | 2:(5); 7:(4,1); 19:(3,2) | 5 | — |
| alternating | `(1,0,-2/5)` | 5 | `A_5` | 3:(5); 23:(3,1,1); discriminant=232^2 | 5 | — |
| cyclic | `(1,0,-7/10)` | 10 | `C_5` | 2:(5); explicit order-five automorphism | 5 | — |
| dihedral | `(2/5,-21/10,2)` | 21 | `D_5` | square discriminant; pair-sum resolvent split 5+5 (both 3:(5)); 11:(2,2,1) | 5 | — |
| Frobenius | `(1/2,3/2,2/5)` | 15 | `F_20` | 29:(5); Dummit resolvent root -13/2; nonsquare discriminant | 5 | — |
| product | `(1,-3/2,-9/2)` | 9 | `K_2 x K_3` | (T^2+T+1)(T^3-T^2-5T+9); cubic irreducible modulo 5 | 1 | — |
| Hasse failure | `(5,-144/5,-188/3125)` | 90000 | `irreducible 2+3` | common quadratic resolvent Q(sqrt(-31)); cubic irreducible modulo 5 | 1 | 2: quadratic splits; 31: cubic simple root 15; all other finite primes: unramified common-resolvent argument |
<!-- END GENERATED FIXED QUINTIC LEDGER -->

Thus the table contains all three quintic real signatures, all five
transitive subgroups of \(S_5\), split and quadratic-times-cubic algebras,
and an everywhere locally soluble fiber with no rational point.  Every row
has \(\Pi\ne0\) and nonzero discriminant, so reconstruction identifies the
entire Keller fiber with the displayed etale algebra.

The five transitive-group headline rows have common projective-height bound
\(21\).  Their exact certificates are oracle-free; the assertion that the
first individual heights are \(1,5,10,15,21\) is the separate bounded
PARI/GP computation recorded in the
[height-21 witness card](UNIVERSAL_QUINTIC_CALCULATOR.md).

For the Hasse row,

\[
\widetilde E(T)=(T^2-8T+47)(T^3+8T^2+12T+8).
\]

The factor discriminants are \(-31\cdot2^2\) and \(-31\cdot8^2\).
Consequently only \(2\) and \(31\) require special treatment.  At \(2\),
\(-31\equiv1\pmod8\), so the quadratic splits over \(\mathbb Q_2\).  At
\(31\), the cubic has the simple root \(15\).  At every other finite prime,
the cyclic decomposition group in the common \(S_3\) splitting field fixes
a root of one factor.  The cubic supplies a real root, neither factor has a
rational root, and hence this is a Hasse failure.  Its primitive projective
target \([3125:15625:-90000:-188]\) has height \(90000\); no global
height-minimality is claimed.

The seven unramified partitions of five already occur modulo \(7\):

| partition | \((\bar\Pi,\bar B,\bar C)\) |
|---|---|
| \((5)\) | \((1,0,1)\) |
| \((4,1)\) | \((1,0,3)\) |
| \((3,2)\) | \((1,0,2)\) |
| \((3,1,1)\) | \((1,1,3)\) |
| \((2,2,1)\) | \((3,2,0)\) |
| \((2,1,1,1)\) | \((1,2,6)\) |
| \((1,1,1,1,1)\) | \((1,0,0)\) |

The coefficient map

\[
(c_2,c_3,c_4)=(-2\Pi B,4\Pi^3,-2\Pi^5C)
\]

has Jacobian \(-48\Pi^8\).  This is the concise geometric dominance
certificate; it is separate from the finite arithmetic rows above.

For provenance, the bounded trace search that found the earlier transported
\(F_{20}\) row is

```bash
.venv/bin/python scripts/search_fixed_quintic_trace_points.py \
  --u -10 --v 20 --bound 18
```

The unified checker generates both this table and
[`fixed_quintic_certificate_ledger.json`](../artifacts/generated-results/fixed_quintic_certificate_ledger.json),
then runs the four underlying exact checkers:

```bash
.venv/bin/python scripts/verify_fixed_quintic_certificate_ledger.py
```

To refresh the generated JSON and Markdown after an intentional row change,
run

```bash
.venv/bin/python scripts/verify_fixed_quintic_certificate_ledger.py --write
```

The underlying commands are

```bash
.venv/bin/python scripts/verify_fixed_quintic_moduli_dominance.py
.venv/bin/python scripts/verify_fixed_quintic_arithmetic_zoo.py
.venv/bin/python scripts/verify_universal_quintic_calculator.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_thirty_one.py
```

The bounded PARI/GP height search is reproduced by

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_targets.py
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py
```

The first command's default box and the limitation to bounded search evidence
are stated in the canonical note.  The second varies the squarefree common
quadratic discriminant and found the \(\mathbb Q(\sqrt{-7})\) row recorded
in the canonical note.
The two earlier Hasse rows, with common resolvents
\(\mathbb Q(\sqrt{-3})\) and \(\mathbb Q(\sqrt{-7})\), remain useful
independent regressions.  Their exact audits are

```bash
.venv/bin/python scripts/verify_fixed_quintic_arithmetic_zoo.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_seven.py
```

Infinitely many Hasse failures in this fixed pencil remain open.  The
canonical note records an exact rational parametrization of the
common-resolvent threefold.  In proportional conic-parameter slices it finds
elliptic curves of PARI rank two and one.  The accompanying bounded search

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_curves.py
```

tests the remaining cube, irreducibility, and small-prime local-root
conditions; in its default box the only survivors at \(2,3,5\) are
presentations of the known Hasse target.  This is experimental evidence, not
an infinitude proof.  The canonical note also proves a useful negative result:
the standard pure-cubic family
\(\mathbb Q(\sqrt{-3})\times\mathbb Q(\sqrt[3]{m})\) cannot enter the
normalized trace chart, because it would require the conic
\(5v^2-9u^2=15\), which has no \(\mathbb Q_5\)-point.

For the clean \(\mathbb Q(\sqrt{-31})\) row, the first exact curve search
excludes every affine-linear base curve through the certified point and
every degree-at-most-two curve on the coordinate-fixed slices \(A=-8\),
\(R=2\), and \(\Pi=5\).  Its general bounded quadratic continuation tests
15024 genuine coefficient tuples in the box \([-2,2]^6\) and finds no
square pullback:

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_rational_curves.py
```

The exact obstruction and bounded experiment are recorded in
[`fixed_quintic_hasse_curve_search.json`](../artifacts/generated-results/fixed_quintic_hasse_curve_search.json).
Rational curves of higher complexity and infinitely many Hasse failures
remain open.
