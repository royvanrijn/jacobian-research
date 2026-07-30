# Programme 6: rank-stratified moment theory

## 1. Objective and status

Let
\[
 V_d=\operatorname {Sym}^d(k^2)^*\otimes
      \operatorname {Sym}^d(k^2)
     \simeq\operatorname {End}(\operatorname {Sym}^d(k^2)).
\]
In the monomial bases, write \(C_F\) for the \((d+1)\)-by-\((d+1)\)
coefficient matrix of \(F\in V_d\).  The organising invariant is
\[
 \rho(F)=\operatorname {rank}C_F. \tag{1.1}
\]
Because \(V_d\) has two tensor factors, this ordinary matrix rank is
exactly tensor rank.  It is invariant under the diagonal
\(\operatorname {SL}_2\)-action and does not depend on the chosen bases.

Let \(Z_d^\infty\) be the common zero set of all pure contraction moments
\(\mu_m\), and let \({\cal N}_d\) be the one-sided nullcone.  Define
\[
 r_d=\min\{\rho(F):F\in Z_d^\infty\setminus{\cal N}_d\}, \tag{1.2}
\]
with \(r_d=\infty\) if no such point exists.  Thus \(r_d\) measures
moment--nullcone failure in balanced degree \(d\).  Promotion to an SIC
counterexample additionally records one fixed multiplier with a nonzero
mixed tail; the known witnesses have this stronger property.

The certified all-order ledger is
\[
\begin{array}{c|c}
d&\text{current knowledge}\\ \hline
2&r_2=\infty,\\
3&r_3=\infty\ \text{or}\ 2\le r_3\le4,\\
4&2\le r_4\le5,\\
d\ge5&2\le r_d\le d+1.
\end{array} \tag{1.3}
\]
The lower bound is the split-symbol theorem: rank one is safe in every
degree.  The degree-four upper bound is the known full-rank witness.  The
radially propagated witnesses give the final row, without asserting that
their displayed coefficient ranks are minimal.

Ranks two through four remain unclassified in degree four.  No
finite-prefix calculation below changes that all-order statement.

For SIC discovery, the split-symbol theorem changes the opening
completely.  Rank one is not an experimental frontier: the entire Segre
cone is already safe in every degree, without a root-partition or finite
moment hypothesis.  The bidegree-\((4,4)\) programme therefore starts at
rank two:

1. parameterize an exact-rank-\(r\) determinantal chart, beginning with
   \(r=2\);
2. impose enough pure moments to reach the expected quotient dimension
   and extract exact closed points or positive-dimensional components;
3. specialize the factor-period representation on those exact
   points/components and use holonomic creative telescoping to derive a
   checkable all-order recurrence; and
4. test one fixed mixed multiplier on every all-order pure survivor.

The known sixteen-term witness has rank five, while exact ranks two,
three, and four remain open.  A rank-two counterexample would therefore
replace the present witness by two separated channels.  Broad support
searches and generic binary-symbol orbit closure are not the primary
search variables for this question.

## 2. Why this is the common SIC/GVC/GMC filtration

A rank-one point has the separated form
\[
 F=A(\xi)P(z), \tag{2.1}
\]
so its pure moments are the homogeneous GVC expressions
\[
 \mu_m(F)=A(\partial)^mP^m. \tag{2.2}
\]
Rank \(r\) replaces one separated channel by
\[
 F=\sum_{q=1}^r A_q(\xi)P_q(z). \tag{2.3}
\]
This is the first determinantal enlargement of homogeneous GVC inside
balanced SIC.  Circular Gaussian expectation identifies the same
contraction moments with balanced complex Gaussian moments.  Any transfer
to a real GMC polynomial must nevertheless track complexification,
polarization, and rank change explicitly; rank is not declared invariant
under an unrelated reduction.

Thus the present programme studies the common balanced moment core.  SIC,
GVC, and GMC consequences are attached only after the relevant
multiplier or reduction has been verified.

## 3. Determinantal strata

Put
\[
 X_{d,r}=\{F\in V_d:\rho(F)\le r\}.
\]
The exact-rank stratum is
\[
 X_{d,r}^{\circ}=X_{d,r}\setminus X_{d,r-1}. \tag{3.1}
\]
The factor chart
\[
 C=UW,\qquad
 U\in\operatorname {Mat}_{d+1,r},\quad
 W\in\operatorname {Mat}_{r,d+1} \tag{3.2}
\]
parametrizes \(X_{d,r}^{\circ}\) when both factors have rank \(r\); this
is enforced on affine charts by inverting one \(r\)-by-\(r\) minor of
each factor.  It has the generic \(\operatorname {GL}_r\)-gauge.  Hence
\[
 \dim X_{d,r}=r(2(d+1)-r). \tag{3.3}
\]
For \(d\ge2\), the generic diagonal \(\operatorname {SL}_2\)-stabilizer is
finite on the cases used below, so the quotient dimension is
\[
 q_{d,r}=r(2(d+1)-r)-3. \tag{3.4}
\]

The \(m\)-th moment is
\[
 \mu_m(C)=
 \sum_{I=0}^{dm}(dm-I)!\,I!\,
 [x^Iy^I]\left(\sum_{i,j=0}^dc_{ij}x^iy^j\right)^m. \tag{3.5}
\]
After (3.2), it is a gauge-invariant polynomial of degree \(2m\) in the
factor variables and degree \(m\) on \(X_{d,r}\).

## 4. Exact Hilbert test

The degree-\(n\) coordinate ring of \(X_{d,r}\) has the Cauchy
decomposition
\[
 k[X_{d,r}]_n=
 \bigoplus_{\substack{\lambda\vdash n\\\ell(\lambda)\le r}}
 S_\lambda(\operatorname {Sym}^d)^*
 \otimes S_\lambda(\operatorname {Sym}^d). \tag{4.1}
\]
For an \(\operatorname {SL}_2\)-character with even weights, invariant
multiplicity is the weight-zero multiplicity minus the weight-two
multiplicity.  This gives the exact invariant Hilbert series
\[
 H_{d,r}(t)=
 H_{k[X_{d,r}]^{\operatorname {SL}_2}}(t). \tag{4.2}
\]

Suppose \(q_{d,r}\) homogeneous invariants of degrees
\(e_1,\ldots,e_{q_{d,r}}\) were a homogeneous system of parameters.
Cohen--Macaulayness would force
\[
 H_{d,r}(t)\prod_i(1-t^{e_i}) \tag{4.3}
\]
to have nonnegative coefficients.  A negative coefficient therefore
proves that their common zero fiber contains a semistable point.  It does
not make that point an all-order moment zero.

The checker evaluates (4.1) by the principal-specialization formula for
\(S_\lambda(\operatorname {Sym}^d)\), rather than by sampling or numerical
integration.

## 5. Low-degree census

The following table is exact.  `Jac` is the Jacobian rank of the
consecutive moments \(\mu_1,\ldots,\mu_q\), evaluated modulo the good
prime \(1000003\).  Equality `Jac=q` proves algebraic independence in
characteristic zero.  The last column gives the first negative
coefficient of (4.3), or records only that none occurs through degree
\(85\).

| \(d\) | \(r\) | \(q_{d,r}\) | Jac | consecutive numerator |
|---:|---:|---:|---:|---|
| 2 | 1 | 2 | 2 | nonnegative through 85 |
| 2 | 2 | 5 | 5 | \([t^5]=-1\) |
| 2 | 3 | 6 | 6 | nonnegative through 85 |
| 3 | 1 | 4 | 4 | nonnegative through 85 |
| 3 | 2 | 9 | 9 | \([t^{29}]=-58\) |
| 3 | 3 | 12 | 12 | nonnegative through 85 |
| 3 | 4 | 13 | 13 | \([t^{63}]=-2186\) |
| 4 | 1 | 6 | 6 | nonnegative through 85 |
| 4 | 2 | 13 | 13 | \([t^{69}]=-5266\) |
| 4 | 3 | 18 | 18 | nonnegative through 85 |
| 4 | 4 | 21 | 21 | nonnegative through 85 |

The degree-three rank-two row is new:

> **Proposition 5.1.** The common zero fiber of
> \(\mu_1,\ldots,\mu_9\) on \(X_{3,2}\) contains a semistable point of
> exact coefficient rank two.

Indeed, the negative coefficient and Cohen--Macaulayness give a
semistable point of rank at most two.  The balanced cubic rank-one theorem
says that already \(\mu_1,\ldots,\mu_4\) cut the rank-one Segre cone down
to the nullcone.  The semistable point therefore cannot have rank one.

This is an existential finite-prefix point.  It has no recorded
coordinates or residue field, and it need not satisfy \(\mu_m=0\) for
\(m\ge10\).

## 6. Least one-degree repairs

For every row with a negative consecutive numerator, the checker replaces
only the last moment and tests increasing orders.  The first
Hilbert-compatible, full-Jacobian systems are
\[
\begin{array}{c|c}
(d,r)&\text{first single replacement}\\ \hline
(2,2)&1,2,3,4,6,\\
(3,2)&1,2,3,4,5,6,7,8,12,\\
(3,4)&1,2,\ldots,12,14,\\
(4,2)&1,2,\ldots,12,14.
\end{array} \tag{6.1}
\]
In each case the candidate numerator is nonnegative through degree \(85\),
its last observed nonzero term is at the Gorenstein-predicted top degree,
and it is palindromic through that degree.  These are necessary Hilbert
tests, not proofs of a homogeneous parameter system or zero-fiber
classifications.  In particular, orders ten and eleven do not repair the
degree-three rank-two system; order twelve is the first single
replacement that passes the tested conditions.

## 7. Research gates

The split-symbol theorem removes rank one from the all-order SIC search.
Finite rank-one moment--nullcone calculations remain useful for the
stronger problem of classifying truncated fibers, but they are not a
prerequisite for beginning the rank-two recurrence programme.

The efficient order is:

1. **Degree four, exact rank two.** Work on full-rank factor charts in
   (3.2), impose a dimension-sized moment system, and compute exact
   determinantal components or closed points.  The existential
   thirteen-moment point on \(X_{4,2}\) is motivation, not a substitute
   for this extraction: it has no coordinates and has not been proved to
   lie in \(X_{4,2}^{\circ}\).
2. **Rank-two recurrence extraction.** For each exact point or component,
   specialize the constant-term/beta period, derive a scalar recurrence
   by creative telescoping, include endpoint certificates and
   singular-step checks, and evaluate enough exact initial and bridge
   moments to decide the all-order tail.
3. **Mixed SIC test.** On every all-order pure survivor, test a fixed
   low-bidegree multiplier and prove an infinite nonzero mixed tail.  A
   pure recurrence alone proves only moment-zero status.
4. **Degree three companion laboratory.** Extract a point or component
   from the existential exact-rank-two nine-moment fiber, evaluate
   \(\mu_{10},\mu_{11},\mu_{12}\), and connect it to the existing
   rank-two holonomic probes.
5. **Ranks three and four.** Move upward only after a global rank-two
   exclusion or after the rank-two components have been classified.
   Consecutive moments pass the Hilbert test through degree \(85\) here,
   so corrected degree sets and exact component geometry should precede
   numerical solving.
6. **All degree.** Compute \(H_{d,r}\), generic moment transcendence
   degree, and the first Hilbert-compatible degree set for fixed
   \(r=2,3,4\) as functions of \(d\).  The structural target is a uniform
   low-rank theorem, not a growing list of bounded searches.

Promotion rules remain strict: a bounded moment survivor is an
experiment or finite-prefix theorem; an all-order pure identity is a
moment-zero result; and an SIC counterexample additionally needs one
fixed multiplier with infinitely many nonzero mixed moments.

## 8. Reproduction

Run

```bash
python3 scripts/verify_rank_stratified_moment_census.py
```

The checker writes
`artifacts/generated-results/rank_stratified_moment_census.json`.  It
reproduces the previously certified degree-four rank-two coefficient
\(-5266\) and the ambient degree-three coefficient \(-2186\) as
regressions.  The degree cutoff \(85\) contains the
Gorenstein-predicted top degree of every repaired system in (6.1).
Absence of a negative coefficient or of a later nonzero coefficient
means only absence through this cutoff.
