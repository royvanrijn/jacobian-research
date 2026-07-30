# Frobenius and \(p\)-curvature of the SIC2C4 moment recurrence

## 1. Status and conclusion

This note tests whether the exact valuation phase diagram of the SIC2C4
mixed moments can be recovered from \(p\)-curvature.  It separates two
notions that should not be conflated:

1. the classical differential \(p\)-curvature of the Picard--Fuchs
   operator for the normalized angular period; and
2. the recurrence-operator \(p\)-curvature, namely the \(p\)-step action
   of the shift \(\tau:m\mapsto m+1\).

The experiment gives a negative answer to the proposed replacement in its
literal form.

* Every first radial lift \(R^k\widetilde F\) has the same normalized
  angular Picard--Fuchs operator.  Its differential \(p\)-curvature is
  nonzero nilpotent of rank one at every odd prime through \(101\).
  Neither the operator nor this curvature sees \(k\), while the valuation
  phase does.
* The minimal full-moment recurrence has order one.  At every sufficiently
  large prime its recurrence \(p\)-curvature is
  \[
     d^d Z^d,\qquad Z=m^p-m,
  \]
  independent of the seed power \(r\).  In particular, the two degree-eight
  families \(R^4\widetilde F\) and \(\widetilde F^2\) have identical
  recurrence \(p\)-curvature.  The first has negative valuation jumps and
  prime-power re-entry; the second has neither.
* The exact bridge is the **uncancelled local singularity ledger** of the
  coprime first-order recurrence.  Its numerator and denominator factors
  give the consecutive valuation increment.  Passing to \(p\)-curvature
  takes their shift norm and cancels precisely the zero/pole information
  needed for re-entry.

Thus \(p\)-curvature supplies a useful global consistency check, but it
does not replace the digit-sum or factor-valuation mechanism.  The
first-order recurrence does replace most family-specific factorial
manipulation: after it is derived once, exact \(p\)-adic propagation is a
local valuation automaton on its finitely many linear factors.

The recurrence formulas and their good-prime shift curvatures below are
exact.  The differential \(p\)-curvature statement is an exact bounded
calculation through \(p=101\), not an all-prime theorem.

## 2. Integral moments and the minimal recurrence

Use the denominator-cleared seed \(\widetilde F=2F\).  For
\[
 d=4r+k,\qquad r\geq1,\quad k\geq0,\qquad
 \widetilde G_{r,k}=R^k\widetilde F^{\,r},
\]
the characteristic-\(p\) note proves
\[
 M_{d,r}(m)
 :=\mathcal E_2(Z\widetilde G_{r,k}^{\,m})
 =4^{rm}\frac{(dm+2)!((rm)!)^2}{(2rm+1)!}.              \tag{2.1}
\]
Direct division at consecutive orders gives
\[
\boxed{\displaystyle
 \frac{M_{d,r}(m+1)}{M_{d,r}(m)}
 =
 4^r
 \frac{
  \prod_{j=3}^{d+2}(dm+j)
  \left(\prod_{j=1}^{r}(rm+j)\right)^2
 }{
  \prod_{j=2}^{2r+1}(2rm+j)
 }.}                                                     \tag{2.2}
\]
Cancel the numerator and denominator over \(\mathbb Q[m]\), and clear
their common rational content.  This produces coprime primitive
\(A_{d,r},B_{d,r}\in\mathbb Z[m]\) with
\[
 \boxed{
 A_{d,r}(m)M_{d,r}(m+1)-B_{d,r}(m)M_{d,r}(m)=0.}          \tag{2.3}
\]
This recurrence is minimal in order: the sequence is nonzero at every
order in characteristic zero, so no nonzero order-zero operator can
annihilate it.  Coprimality makes (2.3) the primitive minimal
order-one representative, up to sign.

The first rows are
\[
\begin{array}{c|c|c}
(d,r)&A_{d,r}(m)&B_{d,r}(m)\\ \hline
(4,1)&1&
16(m+1)^2(4m+3)(4m+5)\\
(5,1)&2m+3&
10(m+1)^2(5m+3)(5m+4)(5m+6)(5m+7)\\
(8,2)&1&
1024(m+1)^2(2m+1)^2
(8m+3)(8m+5)(8m+7)(8m+9)\\
(9,2)&(4m+3)(4m+5)&
648(m+1)^2(2m+1)(3m+1)(3m+2)\\
&&{}\qquad\cdot(9m+4)(9m+5)(9m+7)(9m+8)(9m+10)(9m+11).
\end{array}                                               \tag{2.4}
\]
The generated artifact records nine such factorizations through
\((d,r)=(13,3)\).

For the radial-power subfamily \(d=4r\), (2.1) is \(A_{rm}\), where
\[
 A_s=4^s\frac{(4s+2)!(s!)^2}{(2s+1)!}.
\]
The already proved quotient
\[
 \frac{A_{s+1}}{A_s}
 =16(4s+3)(4s+5)(s+1)^2                         \tag{2.5}
\]
shows directly that \(A_{4r,r}=1\): advancing \(m\) multiplies \(r\)
polynomial quotients.  Hence every \(p\)-adic valuation is nondecreasing
on this subfamily.

## 3. The normalized Picard--Fuchs operator

Dividing (2.1) by its radial factorial isolates the angular period
\[
 b_r(m)=\frac{(rm)!}{(2rm+1)!!},\qquad
 M_{d,r}(m)=(dm+2)!\,2^{rm}b_r(m).                       \tag{3.1}
\]
The decisive separation is already visible here:
\[
 b_r(m)\ \text{depends on }r,\qquad
 (dm+2)!\ \text{depends on }d=4r+k.                      \tag{3.2}
\]
In particular, every \(R^k\widetilde F\) has \(r=1\), so every first
radial lift has the same angular period up to the harmless dilation
\(x\mapsto2x\).

For \(r=1\),
\[
 (2m+3)b_1(m+1)=(m+1)b_1(m).                             \tag{3.3}
\]
If
\[
 y(x)=\sum_{m\geq0}b_1(m)x^m,
\]
then its minimal homogeneous differential operator is
\[
 \boxed{
 L_{\mathrm{PF}}
 =\theta(2\theta+1)-x(\theta+1)^2,\qquad
 \theta=x\frac d{dx}.}                                   \tag{3.4}
\]
Indeed, the first-order coefficient recurrence gives the inhomogeneous
equation
\[
 ((2\theta+1)-x(\theta+1))y=1,
\]
and applying \(\theta\) yields (3.4).  Equivalently,
\[
 L_{\mathrm{PF}}
 =x^2(2-x)D^2+3x(1-x)D-x,                                \tag{3.5}
\]
so its singular support is \(x=0,2,\infty\).  This is the beta-period
operator behind the SIC2C4 angular calculation.

The factorial multiplier in (3.1) is not a gauge transformation of this
connection.  It is an inverse-Borel/Laplace-type operation on
coefficients and changes the Gevrey growth.  Treating the full integral
moment series and the normalized period as the same Picard--Fuchs module
therefore discards the part that creates factorial divisibility.

## 4. Differential \(p\)-curvature experiment

For each odd prime, the checker reduces the companion connection of
(3.5) over \(\mathbb F_p(x)\).  If its connection matrix is \(D+C\), it
computes the zero-order matrix of \((D+C)^p-D^p\) by
\[
 C_1=C,\qquad C_{n+1}=C_n'+CC_n.                          \tag{4.1}
\]
The result for every
\[
 p=3,5,7,\ldots,101
\]
is:

* the \(p\)-curvature is nonzero;
* its trace and determinant vanish;
* its square is zero, so it has generic rank one and nilpotency index two;
* every denominator factor is supported at \(x=0\) or \(x=2\).

This is compatible with global nilpotence of a period connection, but the
bounded calculation is not used as an all-prime proof.  More importantly
for the proposed bridge, exactly the same calculation applies to every
\(R^k\widetilde F\).  The \(d\)-dependent cutoff
\[
 dm+2<p                                                     \tag{4.2}
\]
and the \(d\)-dependent prime-power phase diagram cannot be reconstructed
from an invariant that is literally unchanged when \(k\) changes.

## 5. Recurrence \(p\)-curvature

There is also a genuine \(p\)-curvature for recurrence operators.  Work
over \(\mathbb F_p(m)[\tau]\) with
\[
 \tau f(m)=f(m+1)\tau .
\]
For the rank-one operator
\[
 A(m)\tau-B(m),
\]
the \(p\)-curvature is the \(\mathbb F_p(m)\)-linear action of
\(\tau^p\), hence the multiplier
\[
 \Psi_p(m)
 =\prod_{i=0}^{p-1}\frac{B(m+i)}{A(m+i)}.                 \tag{5.1}
\]
Its characteristic polynomial is \(T-\Psi_p\).  This is the recurrence
notion used by Zhou and van Hoeij; the constants of the shift are
\(\mathbb F_p(Z)\), where
\[
 Z=m^p-m.
\]

For a good linear factor \(am+b\), with \(a\ne0\pmod p\),
\[
\begin{aligned}
 \prod_{i=0}^{p-1}(a(m+i)+b)
 &=a^p\prod_{i\in\mathbb F_p}(m+i+b/a)\\
 &=a(m^p-m)=aZ.                                           \tag{5.2}
\end{aligned}
\]
Take \(p>d+2\), so every slope in (2.2) is a unit.  Applying (5.2) to
the coprime numerator and denominator shows that their separate shift
norms are scalar powers of \(Z\).  Their exponent difference is
\[
 \deg B_{d,r}-\deg A_{d,r}=d,                              \tag{5.3}
\]
and their leading-coefficient ratio is
\[
 \lim_{m\to\infty}
 \frac{M_{d,r}(m+1)}{m^dM_{d,r}(m)}
 =d^d.                                                     \tag{5.4}
\]
Therefore
\[
\boxed{
 \Psi_{p,d,r}=d^dZ^d,\qquad
 \chi_{p,d,r}(T)=T-d^dZ^d
 \quad(p>d+2).}                                           \tag{5.5}
\]
This is an exact all-good-prime formula.  Direct products were also
computed for the \(d=4,5\) rows at every good prime through \(43\), and
at two good primes for each remaining displayed row.

The angular recurrence has equal numerator and denominator degrees.
The same norm calculation gives
\[
 \Psi^{\rm angular}_{p,r}=2^{-r}                           \tag{5.6}
\]
for the unscaled \(b_r\).  Both (5.5) and (5.6) are generically nonzero,
so these rank-one recurrence curvatures are not nilpotent.

Formula (5.5) is too coarse for the valuation phase.  Every good linear
factor lies in the single shift orbit detected by \(Z\); numerator and
denominator factors cancel down to their degree difference.  The
curvature even identifies the following arithmetically different rows:
\[
\begin{array}{c|c|c|c}
(d,r)&A_{d,r}&\Psi_{p,d,r}&\text{valuation behavior}\\ \hline
(8,1)&2m+3&8^8Z^8&\text{negative jumps and re-entry}\\
(8,2)&1&8^8Z^8&\text{monotone; no re-entry}.
\end{array}                                               \tag{5.7}
\]
Thus neither nilpotency nor the reduced singular factor of the
\(p\)-curvature classifies the exact phase.

## 6. What does recover the phase

Taking valuations in the primitive recurrence (2.3) gives the local rule
\[
\boxed{
 v_p(M_{d,r}(m+1))-v_p(M_{d,r}(m))
 =v_p(B_{d,r}(m))-v_p(A_{d,r}(m)).}                        \tag{6.1}
\]
Because \(A\) and \(B\) are explicitly factored into affine-linear
terms, (6.1) is a finite-state base-\(p\) transition rule.  It is the
factorized version of the signed block formula already proved in
(4.14e) of
[`TWO_PAIR_SIC_CHARACTERISTIC_P.md`](TWO_PAIR_SIC_CHARACTERISTIC_P.md).

The known re-entry becomes transparent.  For \((d,r)=(5,1)\),
\[
 A_{5,1}(m)=2m+3.
\]
At \(p=11,m=4\),
\[
 v_{11}(B_{5,1}(4))=0,\qquad
 v_{11}(A_{5,1}(4))=v_{11}(11)=1.
\]
Hence (6.1) gives the negative jump
\[
 v_{11}(M_{5,1}(5))-v_{11}(M_{5,1}(4))=-1,
\]
and indeed
\[
 v_{11}(M_{5,1}(4))=2,\qquad
 v_{11}(M_{5,1}(5))=1.                                   \tag{6.2}
\]
The order-four moment vanishes modulo \(11^2\), while the order-five
moment reappears.  In the shift norm (5.5), the norm of \(2m+3\) is one
factor of \(Z\), which cancels against a numerator factor of \(Z\).
The \(p\)-curvature has forgotten the pole that caused (6.2).

As a bounded correlation audit, the checker used every odd prime through
\(101\), orders through \(2p\), and moduli \(p,p^2,p^3\).  The counts are
\[
\begin{array}{c|r|r}
(d,r)&\text{negative consecutive jumps}&
\text{vanish--reappear transitions}\\ \hline
(4,1)&0&0\\
(5,1)&44&22\\
(6,1)&46&23\\
(7,1)&40&20\\
(8,1)&44&22\\
(8,2)&0&0\\
(9,2)&75&18\\
(12,3)&0&0\\
(13,3)&93&8.
\end{array}                                               \tag{6.3}
\]
These are regression counts, not asymptotic densities.  Their role is to
show that the denominator-free radial rows and the denominator-bearing
nonradial rows behave differently even when (5.5) is identical.

## 7. Revised arithmetic mechanism

For this family, the reusable mechanism is:

1. normalize the contraction moment into a period part and a factorial
   part;
2. derive and primitively cancel the first-order recurrence;
3. retain the numerator and forward-denominator factors separately;
4. use (6.1) as the exact \(p\)-adic transition rule;
5. use differential or recurrence \(p\)-curvature only as a global
   module/singularity consistency check.

This is genuinely more general than reapplying Legendre's formula to a
closed factorial quotient: a creative-telescoping recurrence can feed the
same local valuation automaton even when no compact factorial formula is
known.  What cannot be discarded is the integral lattice and its
individual recurrence singularities.  Classical \(p\)-curvature is
insensitive to the factorial inverse-Borel step, while recurrence
\(p\)-curvature norms away the zero/pole balance.

For higher-order recurrences, the analogous next experiment should keep a
Smith/Dieudonne-style integral lattice for the companion matrix and track
the valuations of each step matrix before taking the \(p\)-fold product.
Slopes of that integral Frobenius object may classify phase transitions.
The generic characteristic polynomial of \(p\)-curvature alone will not.

## 8. Placement in the existing programme

The present computation should be consumed rather than widened by a larger
prime table.  The degree-eight control is exact and already separates
generic curvature from the valuation phase, so extending the same bounded
differential calculation beyond \(p=101\) would not change the decision.

The bridge is now used in three places:

1. the characteristic-\(p\) phase diagram takes the primitive recurrence
   ledger as its canonical digit-automaton input;
2. the general holonomic workflow requires primitive integral
   period/raw-moment forms and local factor/Smith ledgers after creative
   telescoping; and
3. the bidegree-\((3,3)\) and rank-two bidegree-\((4,4)\) programmes have
   this arithmetic postprocessing specified after their current
   characteristic-zero recurrence gates.

The next expansion is therefore conditional and concrete: once either
rank-two programme produces a certified higher-order recurrence, compute
the integral companion-lattice step invariants and compare them with its
generic \(p\)-curvature.  Until then, more SIC2C4 primes or more radial rows
would be redundant calibration rather than progress on an open gate.

## Reproduction

Run

```bash
.venv/bin/python scripts/research_two_pair_sic_frobenius_curvature.py
```

The script:

* derives and verifies the primitive recurrences at nine propagated rows;
* computes direct recurrence \(p\)-curvature products at the recorded good
  primes and compares them with (5.5);
* computes the differential \(p\)-curvature of (3.4) at every odd prime
  through \(101\);
* verifies its nonzero square-zero form and pole support;
* audits the local valuation recurrence and bounded prime-power re-entry
  table; and
* writes
  [`two_pair_sic_frobenius_curvature.json`](../artifacts/generated-results/two_pair_sic_frobenius_curvature.json).

The recurrence \(p\)-curvature convention and its relation to true
singularities are described in Y. Zhou and M. van Hoeij,
[*Desingularization and \(p\)-Curvature of Recurrence
Operators*](https://arxiv.org/abs/2202.08931).  General global-nilpotence
context for differential connections is supplied, for example, by
M. Dettweiler and S. Reiter,
[*On globally nilpotent differential
equations*](https://arxiv.org/abs/math/0605383).
