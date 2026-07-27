# Fixed-quintic arithmetic-zoo certificate

The canonical theorem and proof are consolidated in
[`FIXED_QUINTIC_MODULI_DOMINANCE.md`](../FIXED_QUINTIC_MODULI_DOMINANCE.md#61-one-explicit-arithmetic-zoo).
This companion records the arithmetic certificate and its direct
reproduction command.

For the single determinant-`-2` degree-five Keller map with inverse pencil

\[
E_{\Pi,B,C}(S)=\Pi^5S^5-5\Pi S^3-2BS^2+4S-2C,
\]

the rows are:

| type | target \((\Pi,B,C)\) | inverse polynomial |
|---|---|---|
| \(\mathbb Q^5\) | \((1,0,0)\) | \(S(S-1)(S+1)(S-2)(S+2)\) |
| irreducible \(S_5\) | \((1,0,-1/2)\) | \(S^5-5S^3+4S+1\) |
| irreducible \(A_5\) | \((1,-4/3,6)\) | \(S^5-5S^3+\frac83S^2+4S-12\) |
| \(K_2\times K_3\) | \((1,-3/2,-9/2)\) | \((S^2+S+1)(S^3-S^2-5S+9)\) |
| Hasse failure | \((4,-335/27,4807/20736)\) | \((192S^2-72S+19)(55296S^3+20736S^2+1224S-253)/10368\) |

The Hasse row has common quadratic resolvent
\(\mathbb Q(\sqrt{-3})\).  Its primitive projective target coordinates are

\[
[20736:82944:-257280:4807],
\]

of height \(257280\).  The proof checks the only possibly bad primes
\(2,3,7,19\) explicitly.

Varying the shared quadratic field gives a second, smaller Hasse row:

\[
\left(\Pi,B,C\right)
=\left(-7,\frac{387}{14},\frac{400}{2401}\right),
\]

with normalized inverse polynomial

\[
(T^2-4T+32)(T^3+4T^2-21T+175).
\]

The two factor discriminants are \(-7\cdot4^2\) and
\(-7\cdot(5\cdot79)^2\).  Exact local witnesses at \(2,5,7,79\) prove
everywhere local solubility.  Its primitive target
\([4802:-33614:132741:800]\) has height \(132741\).

A wider search gives the still smaller row

\[
\left(5,-\frac{144}{5},-\frac{188}{3125}\right),
\qquad
(T^2-8T+47)(T^3+8T^2+12T+8).
\]

Its factor discriminants are \(-31\cdot2^2\) and \(-31\cdot8^2\).
The only exceptional primes are \(2\) and \(31\): the quadratic splits over
\(\mathbb Q_2\), and the cubic has the simple root \(15\) modulo \(31\).
The primitive target \([3125:15625:-90000:-188]\) has height \(90000\).

The same map also contains every transitive quintic Galois group.  In terms
of

\[
\widetilde E_{\Pi,B,C}(T)
=\Pi^5E_{\Pi,B,C}(\Pi^{-2}T),
\]

the three additional solvable rows are:

| group | target \((\Pi,B,C)\) | \(\widetilde E_{\Pi,B,C}(T)\) |
|---|---|---|
| \(C_5\) | \((1,-15/11,331/242)\) | \(T^5-5T^3+\frac{30}{11}T^2+4T-\frac{331}{121}\) |
| \(D_5\) | \((5/2,-27/8,-738/3125)\) | \(T^5-5T^3+\frac{135}{8}T^2+\frac{125}{2}T+\frac{369}{8}\) |
| \(F_{20}\) | \((31/5,5229/310,9618099/114516604)\) | \(T^5-5T^3-\frac{5229}{25}T^2+\frac{119164}{125}T-\frac{9618099}{6250}\) |

Together with the \(A_5\) and \(S_5\) rows above, these give all five
transitive subgroups of \(S_5\).  The exact certificates are respectively:
an explicit order-five cyclotomic automorphism; a two-quintic pair-sum
resolvent plus a \((2,2,1)\) modular factor pattern; and a rational root of
Cayley's sextic resolvent plus a nonsquare discriminant.

The bounded search that found the \(F_{20}\) transport is

```bash
.venv/bin/python scripts/search_fixed_quintic_trace_points.py \
  --u -10 --v 20 --bound 18
```

Run

```bash
.venv/bin/python scripts/verify_fixed_quintic_arithmetic_zoo.py
```

The checker expands the fixed map, verifies its constant Jacobian, every
target substitution, factorization and discriminant, all five transitive
Galois-group certificates, squarefreeness, the common quadratic resolvent,
and every exceptional Hensel witness.

The bounded PARI/GP height search is reproduced by

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_targets.py
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py
```

The first command's default box and the limitation to bounded search evidence
are stated in the canonical note.  The second varies the squarefree common
quadratic discriminant and found the \(\mathbb Q(\sqrt{-7})\) row above.
Its exact independent audit is

```bash
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_seven.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_thirty_one.py
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
