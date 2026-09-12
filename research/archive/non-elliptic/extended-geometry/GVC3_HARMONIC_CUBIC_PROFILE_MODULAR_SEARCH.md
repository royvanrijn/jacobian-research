# Full harmonic-cubic obstruction below \(\Delta^6\)

## 1. Scope and status

Put

\[
 \rho=xy+t^2,
 \qquad
 \Delta=4\partial_x\partial_y+\partial_t^2,
\]

and define

\[
 E_4=\rho(xy-2t^2)-x^2t^2,
 \qquad
 O_3=\rho y-3xt^2.
\tag{1.1}
\]

This calculation treats the complete degree-six family

\[
 \boxed{P=\alpha\rho E_4+H_3O_3,\qquad H_3\in\mathcal H_3.}
\tag{1.2}
\]

On the sphere \(\rho=1\), it is

\[
 F=\alpha E+H_3O,
 \quad
 E=xy-2t^2-x^2t^2,
 \quad
 O=y-3xt^2.
\tag{1.3}
\]

Write \(\mu_m=\mathcal R(F^m)\) for normalized spherical Reynolds
projection.  Since \(P\) is homogeneous of degree six, \(\mu_m\) differs by
a nonzero scalar from \(\Delta^{3m}(P^m)\).  Thus its pure-moment zero locus
is exactly the relevant \(\Delta^3\) zero locus.

The radial-linear part \(\rho H_1\) of a general cubic profile is covered by
the exact independent-parity quartic obstruction.  Family (1.2) is the
complementary genuine harmonic-cubic repair and contains all seven
coordinates of \(\mathcal H_3\).  It is not a one-profile cusp ansatz.

> **Theorem 1.1 — full harmonic-cubic repair obstruction.**  In
> characteristic zero, if the first eight pure Reynolds moments of (1.2)
> vanish, then \(\alpha=0\) and \(H_3\) lies on one of two one-sided harmonic
> planes.  Both planes satisfy the GVC conclusion for every fixed multiplier.
> Consequently (1.2) contains no homogeneous GVC witness for \(\Delta^3\).

The proof uses exact characteristic-zero unit bases for every chart and an
exact radical certificate on the residual boundary.  Singular calculations
over \(\mathbb F_{101},\mathbb F_{103},\mathbb F_{107}\) are retained only as
the modular discovery replay.  No conclusion is inferred from good primes.

This theorem closes the complete genuine harmonic-cubic repair of Long's
profile.  It is not a classification of arbitrary ternary sextics, of
noncoherent vectors in other harmonic summands, or of repeated profiles
inside one summand, and it does not prove that \(\Delta^6\) is globally
minimal.

## 2. Harmonic coordinates and the first moment

Use the weight basis

\[
\begin{array}{c|c}
3&x^3\\
2&x^2t\\
1&x^2y-4xt^2\\
0&3xyt-2t^3\\
-1&xy^2-4yt^2\\
-2&y^2t\\
-3&y^3
\end{array}
\tag{2.1}
\]

of \(\mathcal H_3=\ker(\Delta:\mathbb C[x,y,t]_3\to
\mathbb C[x,y,t]_1)\).  Write the corresponding coefficients as

\[
 (p_3,p_2,p_1,p_0,n_1,n_2,n_3).
\]

Exact Reynolds extraction gives

\[
 \boxed{\mu_1=\frac{16}{35}n_1.}
\tag{2.2}
\]

Thus every pure-zero point has \(n_1=0\).  The search therefore has the
seven remaining channels

\[
 (\alpha,p_3,p_2,p_1,p_0,n_2,n_3).
\tag{2.3}
\]

## 3. Invariant moment compiler

On \(\rho=1\), substitute

\[
 y=\frac{1-t^2}{x}.
\tag{3.1}
\]

If \(G=\sum c_{r,s}x^rt^s\), normalized spherical Reynolds projection is

\[
 \mathcal R(G)=
 \sum_{\substack{s\ge0\\s\text{ even}}}\frac{c_{0,s}}{s+1}.
\tag{3.2}
\]

After (2.2), the seven coefficient channels of \(F\) are the following
Laurent polynomials:

\[
\begin{array}{c|l}
\alpha&1-3t^2-t^2x^2\\
p_3&x^2(1-t^2)-3t^2x^4\\
p_2&x(t-t^3)-3t^3x^3\\
p_1&1-6t^2+5t^4+x^2(15t^4-3t^2)\\
p_0&x(15t^5-9t^3)+x^{-1}(5t^5-8t^3+3t)\\
n_2&x^{-1}(-3t^7+6t^5-3t^3)
      +x^{-3}(-t^7+3t^5-3t^3+t)\\
n_3&x^{-2}(3t^8-9t^6+9t^4-3t^2)
      +x^{-4}(t^8-4t^6+6t^4-4t^2+1).
\end{array}
\tag{3.3}
\]

For an occupation vector \(e=(e_0,\ldots,e_6)\) with \(|e|=m\), the
compiler multiplies the seven small Laurent powers in (3.3), applies (3.2),
and multiplies by \(\binom m e\).  It never expands \(P^m\) in all
coefficients of \(x,y,t\).  The resulting numbers of nonzero occupation
terms are

\[
\begin{array}{c|rrrrrrr}
m&2&3&4&5&6&7&8\\ \hline
\#\operatorname{supp}\mu_m&10&27&74&162&331&616&1093.
\end{array}
\tag{3.4}
\]

This is the invariant moment system used below.

## 4. The projective cover on \(\alpha\ne0\)

Normalize \(\alpha=1\).  A primitive integer multiple of the second moment
is

\[
\begin{aligned}
M_2={}&720n_2p_0-1872n_2p_2-2592n_3p_1-4992n_3p_3+936n_3\\
&-4524p_0^2+2288p_0p_2+4576p_1^2+10296p_1+9009.
\end{aligned}
\tag{4.1}
\]

Its two pivot coefficients are

\[
 D=936-2592p_1-4992p_3,
 \qquad
 K=720p_0-1872p_2.
\tag{4.2}
\]

The projective chart is covered by

\[
 D\ne0,
 \qquad
 D=0, K\ne0,
 \qquad
 D=K=0.
\tag{4.3}
\]

For modular discovery on \(D\ne0\), equation (4.1) eliminates \(n_3\); each
later moment is cleared by the smallest required power of \(D\), and
\(zD-1\) records the localization.  The exact replay instead keeps the
original moments and adjoins \(zD-1\) directly.  On the second stratum,
\(D=0\) and \(zK-1\) are adjoined.  On the last stratum, both linear
equations are imposed.  The results are:

\[
\begin{array}{c|c|c|ccc}
\text{stratum}&\text{moments used}&\mathbb Q&101&103&107\\ \hline
D\ne0&2,\ldots,8&(1)&(1)&(1)&(1)\\
D=0, K\ne0&2,\ldots,7&(1)&(1)&(1)&(1)\\
D=K=0&2,\ldots,7&(1)&(1)&(1)&(1).
\end{array}
\tag{4.4}
\]

Here \((1)\) means that the ideal is the unit ideal.  In the
characteristic-zero column, msolve returns the literal reduced basis
\([1]\).  Thus the three strata give an exhaustive exact exclusion of
\(\alpha\ne0\); the finite-field columns document discovery and independent
replay.

## 5. The \(\alpha=0\) boundary

Let \(I_7\) be the ideal generated by \(\mu_2,\ldots,\mu_7\) after
\(\alpha=0\), and put

\[
\begin{aligned}
J_+&=(p_1,p_0,n_2,n_3),\\
J_-&=(p_3,p_2,p_1,p_0),\\
J&=J_+\cap J_-\\
 &=(p_1,p_0,p_3n_2,p_3n_3,p_2n_2,p_2n_3).
\end{aligned}
\tag{5.1}
\]

The occupation compiler verifies exactly that every compiled moment
vanishes on both planes, so \(I_7\subseteq J\) over \(\mathbb Q\).  For each
of the six displayed generators \(g\) of \(J\), the calculation adjoins
\(zg-1\).  Exact characteristic-zero F4 elimination with msolve returns the
literal reduced basis \([1]\) in every case:

\[
 I_7+(zg-1)=(1).
\tag{5.2}
\]

Since \(J\) is radical, (5.2) proves over \(\mathbb Q\) that

\[
 \boxed{\sqrt{I_7}=J.}
\tag{5.3}
\]

This proves the boundary classification.  Modular Singular replay gives the
same six unit tests at \(101,103,107\); the unreduced boundary Gröbner basis
has dimension two and 569 elements at each prime.  Equation (5.3), rather
than those large nonreduced bases, is the geometric output.

On \(J_+\), only \(p_3,p_2\) remain and every Laurent monomial in (3.3) has
strictly positive phase.  On \(J_-\), only \(n_2,n_3\) remain and every
monomial has strictly negative phase.  If \(Q\) has degree \(q\), phase
weight therefore gives

\[
 \Delta^{3m}(QP^m)=0\qquad(m>2q)
\tag{5.4}
\]

on both terminal planes.  This all-order terminal argument and the
exhaustion of the \(\alpha=0\) boundary are both exact.

## 6. Candidate gate and conclusion

No nonterminal characteristic-zero point survives the eighth pure moment.
Consequently there is no candidate whose multiplier channel or moment
sequence reaches the promotion gate: no hypergeometric formula, rational
diagonal representation, or P-recursive recurrence is being claimed for a
new witness.  The residual planes have the stronger all-order phase cutoff
(5.4) and hence cannot be witnesses.  This proves Theorem 1.1.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/research_gvc3_harmonic_cubic_profile.py \
  --cases alpha1_n3:8 alpha1_d_k:7 alpha1_dk:7 \
          alpha0:7 alpha0_radical:7 \
  --primes 101 103 107 \
  --timeout 900 --exact-all --msolve-threads 4
```

The command writes
`artifacts/generated-results/gvc3_harmonic_cubic_profile_modular.json`.
It requires Singular and msolve.  The artifact records the nine exact
characteristic-zero msolve unit bases, input hashes, term counts, dimensions,
basis hashes, and all modular discovery replays.
