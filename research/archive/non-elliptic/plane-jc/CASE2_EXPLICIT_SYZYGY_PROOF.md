# The \((72,108)\) Case-2 contradiction by resultants and syzygy

This note replaces the statement “Singular returns the unit ideal” in two
ways.  First it gives a compact projective-resultant proof.  Second it records
the direct characteristic-zero syzygy as an independent fallback.  Both
concern the primary Case-2 coefficient ideal, not the later residue-cover
prefilters.

## Statement

Let

\[
  H(U)\in \mathbf Q[U]
\]

be the irreducible degree-\(35\) polynomial in
`exact_replay/firstblock_Q_exact.sing`, and put

\[
  K=\mathbf Q[U]/(H),\qquad S=K[r,s,h].
\]

After the exact triangular solution of \((J4)\), followed by the exact linear
solution of \((J3)\) and \((J2)\), let

\[
 R_0,R_1,R_7,R_9\in S
\]

be the four residual equations stored in
`exact_replay/case2_exact_certificate.pkl`.  Then there are explicit
polynomials

\[
 T_0,T_1,T_7,T_9\in S
\]

such that

\[
 \boxed{\,1=T_0R_0+T_1R_1+T_7R_7+T_9R_9\,}.                 \tag{1}
\]

Consequently the four residual equations have no common zero over any
extension of \(K\), so the \((72,108)\) Case-2 Laurent system is impossible.

## Compact projective-resultant proof

The support shapes are more informative than a generic four-polynomial
system.  On the chart \(s\ne0\), write

\[
 r=\lambda x,\qquad s=\lambda,\qquad q=\lambda^2.
\]

After division by the common factor \(\lambda\), the two cubics have the
form

\[
 R_i/\lambda=A_i(x)h+B_i(x)q+C_i(x),\qquad i=0,1,
\]

where \(\deg A_i,\deg C_i\le1\) and \(\deg B_i\le3\).  Put

\[
\begin{aligned}
D&=A_0B_1-A_1B_0,\\
\mathcal H&=B_0C_1-B_1C_0,\\
Q&=A_1C_0-A_0C_1.
\end{aligned}
\]

If \(D=0\), consistency requires \(Q=0\).  If \(D\ne0\), Cramer's rule gives

\[
 h=\mathcal H/D,\qquad q=Q/D.
\]

Each quartic \(R_7,R_9\) becomes quadratic in \(h,q\).  Substitution and
multiplication by \(D^2\) therefore gives two univariate polynomials
\(N_7(x),N_9(x)\), both of degree \(8\).

All coefficients lie in a localization of the order
\(\mathbf Q[U]/(H(U))\).  The prime \(101\) misses every denominator in the
minimal polynomial and the four residuals, and

\[
 H(55)=0\pmod {101}.
\]

Thus \(U\mapsto55\) defines a good reduction to \(\mathbf F_{101}\).  Direct
Euclidean identities there give

\[
\gcd(\bar D,\bar Q)=1,\qquad
\gcd(\bar N_7,\bar N_9)=1.                              \tag{2}
\]

The same identities hold on the \(r\ne0\) chart.  At \(r=s=0\), the two
quadratics \(R_7(0,0,h)\) and \(R_9(0,0,h)\) are also coprime modulo \(101\).

Each equality in (2) is an explicit Bézout identity over
\(\mathbf F_{101}[x]\), hence the corresponding resultant has nonzero
reduction.  Its characteristic-zero antecedent is therefore nonzero.  The
singular Cramer branch, the nonsingular branch, the missing projective
direction, and the origin are all excluded.  This proves that the four
\(R_i\) have no common zero.

The checker
[`cas/verify_case2_resultant_proof.py`](cas/verify_case2_resultant_proof.py)
constructs every displayed polynomial from the pinned characteristic-zero
residuals and verifies the extended-gcd identities.  The complete numerical
degree ledger is

\[
\begin{array}{c|c|c|c}
\text{locus}&\deg D,\deg Q&\deg N_7,\deg N_9&\deg\gcd\\ \hline
s\ne0&(4,2)&(8,8)&0\\
r\ne0&(4,2)&(8,8)&0\\
r=s=0&-& (2,2)&0.
\end{array}
\]

## Direct characteristic-zero syzygy

The four equations have the following intrinsic degree shapes:

\[
\begin{array}{c|c|c|c}
 &\deg&\deg_h&\text{number of terms}\\ \hline
R_0,R_1&3&1&8\\
R_7,R_9&4&2&14 .
\end{array}
\]

Thus it is natural to search for \(T_0,T_1\) of degree at most \(5\) and
\(T_7,T_9\) of degree at most \(4\).  Equating coefficients in (1) through
total degree \(8\) is a linear problem over the already proved number field
\(K\).  There are

\[
  \binom{8+3}{3}=165
\]

coefficient rows and

\[
  2\binom{5+3}{3}+2\binom{4+3}{3}=182
\]

candidate multiplier columns.  Exact elimination has rank \(151\) and
places the constant vector \(1\) in their image.  One resulting solution has
respectively \(23,20,14,11\) nonzero multiplier terms.  This is the saved
syzygy (1).

This observation is not being used as a probabilistic rank test.  Every entry
is a rational polynomial in \(U\) of degree below \(35\); addition and
multiplication are performed in \(\mathbf Q[U]/(H)\), and the \(165\)
coefficients of the left side of (1) are checked exactly.

## Direct verification of the syzygy

The first-block calculation proves that \(H\) is irreducible, hence \(K\) is
a field.  The triangular and linear substitutions used to obtain the four
\(R_i\) are exact equivalences and use only nonzero elements of \(K\).

The serialized certificate has SHA-256

```text
cfbc3c39d7a28013671144f43ef76f0498542eaf6d562dd624bba3311194e4aa
```

and contains:

1. all coefficients of \(H\);
2. all coefficients of \(R_0,R_1,R_7,R_9\);
3. all coefficients of \(T_0,T_1,T_7,T_9\);
4. the Macaulay degree, columns, and pivot columns used to find them.

The independent checker
[`cas/verify_case2_syzygy_independent.py`](cas/verify_case2_syzygy_independent.py)
pins the hash, verifies irreducibility of the stored \(H\), checks the degree
and support data above, and multiplies out (1) without importing the
equation generator, `exact_core`, or Singular.  It obtains the zero
polynomial from

\[
 T_0R_0+T_1R_1+T_7R_7+T_9R_9-1.
\]

If a common zero \(p\) of the four residuals existed, evaluation of (1) at
\(p\) would give \(1=0\).  This contradiction proves the claim.

## Reproduction

From the repository root, with the pinned FLINT environment:

```bash
.venv/bin/python plane-jc/cas/verify_case2_resultant_proof.py
.venv/bin/python plane-jc/cas/verify_case2_syzygy_independent.py
```

The expected final lines are

```text
CASE2_PROJECTIVE_RESULTANT_PASS
CASE2_EXPLICIT_SYZYGY_PASS
```

The argument uses no finite-field specialization and no inference from a
standard-basis printout.
