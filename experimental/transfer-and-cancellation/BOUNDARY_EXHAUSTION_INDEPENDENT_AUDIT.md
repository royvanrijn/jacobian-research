# Independent audit of C04--C24 boundary exhaustion

This note independently audits the statement that the displayed C04 and C24
branches are **all** irreducible components of the boundary in the canonical
finite normalization.  It starts from the defining inverse equations and
reconstruction formulas; it does not use the component lists asserted in the
boundary-exhaustion certificate or the resolvent-signature table.

The audit has four targets:

1. identify the normalized inverse incidence with the normalization of the
   target in `k(A^3)`;
2. prove that the complement of the displayed candidate divisors reconstructs
   entirely inside the affine source;
3. classify every normalization prime over the generic points of those
   divisors; and
4. saturate the local field degree by the resulting `sum e f` count.

Throughout, the coefficient field has characteristic zero.  Statements about
individual split roots are geometric; over the original field they descend by
factoring the residual root polynomial.

## 1. Two clean-room lemmas

### Lemma 1.1 (normalization identification)

Let `Y` be normal integral with function field `K`, and suppose an integral
finite `Y`-scheme `V` has function field `L=k(A^3)`.  Then the normalization
of `V` is canonically

\[
 \operatorname{Norm}_Y(L).
\]

If the rational inverse coordinates on this normalization are regular on an
open `R` and give inverse morphisms with the source marking, then
`R` is exactly the distinguished affine-source open.

**Proof.**  Both normalizations are obtained by taking, on every affine open
of `Y`, the integral closure of its coordinate ring in the same field `L`.
The second assertion follows because the two rational compositions are
identities in `L`; wherever all coordinates are regular they are morphism
identities.  QED

### Lemma 1.2 (DVR saturation)

Let `R` be a target DVR, let `L/K` be finite separable of degree `d`, and let
`S` be the integral closure of `R` in `L`.  If constructed distinct primes
`q_i` of `S` satisfy

\[
 \sum_i e_i f_i=d,                                         \tag{1}
\]

then they are every prime of `S`.

**Proof.**  The torsion-free finite `R`-module `S` is free of rank `d`.
The Artinian ring `S/pi S` is the product of its localizations at the primes
over `pi`, and those local dimensions are `e_i f_i`.  Every missing prime
would add a positive summand to (1).  QED

These lemmas refer to the canonical normalization and its DVRs, not to a
chosen compactification.

## 2. C04 from the normalized marked-root incidence

Let

\[
 H(W)=hW^2(W-1)\prod_{j=1}^{n-3}(W-\rho_j)
\]

with the roots distinct and nonzero, and impose the usual open genericity
conditions on the endpoint construction.  The inverse incidence is

\[
 E(W)=H(W)-BCW+cAC^2=0.                                   \tag{2}
\]

Its leading `W`-coefficient is a unit, so (2) is finite flat of degree `n`
over target affine space.  It is integral: as a polynomial linear in `A`, its
coefficients `cC^2` and `H(W)-BCW` are coprime.  On `C E_W !=0`, the formulas

\[
 \gamma={BC-H'(W)\over c},\quad x={C\over\gamma},\quad
 u={W\over\gamma},\quad y={u-1\over x},                    \tag{3}
\]

together with the displayed formula for `z`, invert the source marking.
Thus its function field is `k(A^3)`, and Lemma 1.1 identifies its
normalization with the canonical finite normalization.

The independent C04 audit already derives (2)--(3), the constant Jacobian,
and the equality of the affine source with the simultaneous regularity locus
without importing the weighted-model implementation.  We now inspect every
possible polar prime.

### 2.1 Support away from `C=0`

On `D(C Delta_H)`, all `n` roots are simple.  On the incidence,

\[
 E_W=H'(W)-BC=-c\gamma,                                   \tag{4}
\]

so `C`, `gamma`, `x`, and all denominators in (3) are units or cancel by the
inverse identities.  Every normalized point therefore lies in the affine
source.  Hence every boundary image is contained in

\[
 V(\Delta_H)\cup V(C).                                    \tag{5}
\]

At the generic point of the irreducible discriminant, its normalization is

\[
 W\longmapsto
 (BC,cAC^2)=(H'(W),,WH'(W)-H(W)).                         \tag{6}
\]

The two coordinate degrees are `n-1,n`, so (6) is birational.  The generic
root is exactly double.  Its normalized prime has `(e,f)=(2,1)` and is
boundary by (4), whereas the other `n-2` simple roots reconstruct affinely.
The total contribution is `2+(n-2)=n`.

### 2.2 Every branch over `C=0`

At `C=0`, equation (2) becomes `H(W)=0`; there are no other generic fiber
components.

- At `W=1`, `H'(1)=-c`, hence `gamma=1`.  The étale expansion of (2) in `C`
  makes every expression in (3), including the endpoint-cancelled `z`
  numerator, regular.  This is one affine `(1,1)` prime.
- At the double root `W=0`, put `W=C R`.  Dividing (2) by `C^2` gives
  
  \[
   h_0R^2-BR+cA+O(C)=0,qquad h_0\ne0.                    \tag{7}
  \]
  
  Direct substitution into (3) gives nonnegative valuation for all source
  coordinates.  The generic quadratic in (7) is irreducible over `k(A,B)`,
  so this is one affine prime with `(e,f)=(1,2)`.
- At every additional simple root `rho_j`, equation (2) gives an étale
  normalization branch.  Its generic value in (3) is
  
  \[
   W-\gamma=\rho_j+{H'(\rho_j)\over c}.                   \tag{8}
  \]
  
  The genericity conditions make (8) nonzero, so `y=(W-gamma)/C` has a
  simple pole.  These are `n-3` distinct boundary primes, each `(1,1)`.

Their total contribution over `C=0` is

\[
 2+1+(n-3)=n.                                             \tag{9}
\]

Lemma 1.2 and (5), first over the discriminant and then over `C=0`, prove
that no other canonical boundary prime exists.  Notice that the conclusion
comes from the normalization and degree saturation; it is not inferred from
a denominator list.

## 3. C24 from two monic inverse charts

Put

\[
 N=r(m+1)+1,qquad D(T)=1-T(Q-PT)^m,
\]

and

\[
 \Psi(T)=C\int_0^T D(t)^r\,dt-R.                          \tag{10}
\]

Over `D(P)`, this is a degree-`N` polynomial with unit leading coefficient.
It is irreducible over `k(P,Q,R)`: it is primitive and linear in the
indeterminate `R`.  The exact reconstruction is

\[
 y=Q-TP,quad A=D(T)^{-1},\quad x=T D(T)^{-1},              \tag{11}
\]

with the corresponding polynomial expression for `z` whenever `D(T)!=0`.
Thus the normalized `T`-incidence is the canonical normalization over
`D(P)` by Lemma 1.1.

For the missing divisor `P=0`, put `U=PT`.  Multiplication by `P^(r+1)` gives
the second monic degree-`N` chart

\[
 C J(U,P)-RP^{r+1}=0,qquad
 J(U,P)=\int_0^U\{P-V(Q-V)^m\}^r\,dV.                     \tag{12}
\]

Its function field is again `k(A^3)` because `T=U/P` generically.  Its
normalization is therefore the same canonical normalization near `P=0`.

### 3.1 Support and discriminant prime

Since `Psi_T=C D(T)^r`, a root off the discriminant has `D(T)!=0`.  On
`D(P Delta_(m,r))`, every root is therefore reconstructed by (11), proving

\[
 \operatorname{Supp}(\partial_F)longrightarrow
 V(\Delta_{m,r})\cup V(P).                                \tag{13}
\]

The critical equation `D(T)=0` is integral.  With `Y=Q-PT`, it has

\[
 T=Y^{-m},qquad P=(Q-Y)Y^m.                               \tag{14}
\]

Generic critical values are distinct—for example this is visible after the
specialization `Q=0`, where they differ by `(m+1)`-st roots of unity—so the
critical divisor maps birationally to the discriminant.  Transversely,
`D(T)` has a simple zero and integrating `D(T)^r` gives local order `r+1`.
Thus there is one boundary prime `(e,f)=(r+1,1)`.  It cannot lie in the
étale affine source.  Every other root has `D!=0` and reconstructs by (11),
giving total degree

\[
 (r+1)+N-(r+1)=N.                                         \tag{15}
\]

Lemma 1.2 excludes another prime over the discriminant.

### 3.2 Every prime over `P=0`

At `P=0`, equation (12) factors as

\[
 J(U,0)=(-1)^r U^{r+1}Q^{mr}K_{m,r}(U/Q),                 \tag{16}
\]

where

\[
 K_{m,r}(w)=w^{-(r+1)}\int_0^w v^r(1-v)^{mr}\,dv.          \tag{17}
\]

The degree of `K_(m,r)` is `mr`.  It is squarefree: if a nonzero `a` were a
common zero of `K` and `K'`, then the integral in (17) and its derivative
`a^r(1-a)^(mr)` would vanish.  Hence `a=1`, but the beta integral at `1` is
nonzero.

The `U=0` cluster is analyzed without guessing from (16): set `U=PS` and
divide (12) by `P^(r+1)`.  The exact scaling identity

\[
 J(PS,P)=P^{r+1}\int_0^S\{1-t(Q-Pt)^m\}^r\,dt             \tag{18}
\]

gives at `P=0` an irreducible degree-`r+1` polynomial, linear in the generic
target coordinate `R`.  It is the affine divisor `B=0`, with
`(e,f)=(1,r+1)`.

Among the `mr` nonzero simple roots of (17), the cancellation parameter `q`
selects

\[
 w_q=-{q\over1-q},qquad K_{m,r}(w_q)=0.                   \tag{19}
\]

On the affine source this is exactly `A=0`; it gives one affine `(1,1)`
prime.  There can be no third affine prime over `P=0`, because
`P=AB` and the affine inverse image of `P=0` has only the prime divisors
`A=0` and `B=0`.  Therefore the remaining `mr-1` simple roots are distinct
boundary primes `(1,1)`.  The total is

\[
 (r+1)+1+(mr-1)=N.                                        \tag{20}
\]

Lemma 1.2 excludes every further prime over `P=0`; (13) excludes every other
target image.  This proves the complete C24 boundary list in the canonical
normalization.

## 4. Arithmetic descent and marking safety

The residual polynomials in (7) and (17) are explicit over the coefficient
field.  Factoring them groups the geometric roots into arithmetic primes;
the factor degrees are their residue degrees in these charts.  Applying
Lemma 1.2 over the original field gives the same saturated total, so no new
arithmetic prime can appear during descent.

For every noncubic C04 or C24 member, the canonical target diagram therefore
has exactly two vertices.  One is uniquely marked by receiving a ramified
boundary prime; the other receives only unramified boundary primes.  The
marking used in the C04--C24 obstruction is consequently safe against an
unseen component or a larger automorphism group of a partial diagram.

## 5. Independent executable

Run

```bash
python3 scripts/audit_boundary_exhaustion_independent.py
```

The checker imports only Python's standard library and no project module.  It
independently verifies, for a grid of `m,r`, the derivative and degree of
(10), the chart identity (18), the factorization (16), squarefreeness of
(17), the parameter-root identity (19), and both exact degree saturations.
It also invokes no assertion from the existing boundary certificate.  The
all-parameter conclusions come from the coefficient-independent arguments
above; the grid is a regression against transcription errors.

For the C04 polynomiality, incidence, reconstruction, and constant-Jacobian
identities, the already independent standard-library checker
`scripts/audit_c04_independent.py` supplies the complementary clean-room
regression.
