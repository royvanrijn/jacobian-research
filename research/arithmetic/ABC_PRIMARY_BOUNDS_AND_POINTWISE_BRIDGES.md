# Primary exceptional-set bounds and pointwise transfers

The fixed-budget count, the small-budget limiting coefficient,
and pointwise bounds on the largest letter are different
quantities. The following records fix the source scope of
recent comparisons and give two deductions with their
own proof dependencies.

## A closed-budget consequence of the BBLT program

Use \(N_\lambda(X)\) for ordered positive primitive triples
with \(c\le X\) and \(\operatorname{rad}(abc)<c^\lambda\).
BBLT v2 states the bound \(X^{(23\lambda+3)/40+\varepsilon}\)
for \(0<\lambda\le2\) in Theorem 1.2, and
\(X^{3/5+\varepsilon}\) for fixed \(0<\lambda<1\) in
Theorem 1.3. The 22 inequalities of Proposition 6.1
have the closed budget constraint \(L\le1\).
Use the anatomic reduction and the arithmetic inputs
to its non-budget rows before enforcing that constraint.
[BBLT v2, Theorems 1.2–1.3 and Proposition 6.1](https://arxiv.org/html/2410.12234v2)

**PBP1.** The same normalized program gives
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{\,3/5+(8/5)(\lambda-1)_++\varepsilon}
                   \qquad(0<\lambda\le2).
\]
For \(\lambda\le1\) one may take the constant from
\(\lambda=1\) by monotonicity. For \(1<\lambda<43/41\)
this exponent is strictly below the formula in
BBLT Theorem 1.2.

To verify the endpoint, suppose a fixed budget \(\lambda\)
admits a count exponent strictly above the proposed value.
Pass to anatomic boxes as in Proposition 2.1. For each
truncation tolerance first fix the resulting number of
valuation layers and the smaller count allowances. Then
choose arbitrarily large violating heights. Keep only
the projected variables \(a,b,c,D,a_1,\ldots,c_6\)
of Proposition 6.1, taking a fixed count exponent \(D\)
strictly between the proposed bound and the violating
exponent. The variables are bounded, the radical
budget satisfies \(L\le\lambda+o(1)\), and the fixed
arithmetic inputs imply (6.2)–(6.10) up to errors
tending to zero with the truncation and count tolerances.
Those errors arise without the assumption \(L\le1\);
in the source proof that assumption is only used for
(6.1). Each minimum row is a finite union of closed
linear alternatives. Compactness therefore gives a
limit satisfying all the projected inequalities exactly
except that (6.1) is replaced by \(L\le\lambda\).
This uses an inner arbitrarily large height choice at
every fixed truncation.

If \(\lambda\le1\), apply Proposition 6.1 directly.
If \(\lambda>1\), set \(\theta=1/\lambda\), scale all
radical variables by \(\theta\), and replace \(D\)
by \(D'=\theta D-(1-\theta)\). A negative \(D'\)
already gives the required bound. Rows (6.2)–(6.4)
and (6.10) pass because the height inequalities
contain the fixed constant \(1\) and the decrease
in \(D'\) only relaxes those rows. In every geometry
row, write the two alternatives as
\(D\le L-S\) or \(D\le L-1+T\), with \(S,T\)
linear in the radical variables. The scaled first
alternative follows from \(D'\le\theta D\).
The second follows from
\(D'\le\theta L-1+\theta T\).
Hence all hypotheses of Proposition 6.1 hold for
the scaled variables, giving \(D'\le3/5\).
Solving for \(D\) gives \(D\le8\lambda/5-1\).
Comparing this to \((23\lambda+3)/40\) gives
the indicated cutoff \(43/41\). This is a deduction
from the source program. Later results in the
repository give stronger fixed-budget counts.
[Arithmetic projection and closed program, BBLT v2, §6](https://arxiv.org/html/2410.12234v2#S6)
· [Stronger fixed-budget counts](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#exceptional-set-bounds-from-the-two-tall-terms)

## A two-parameter inversion of Pasten's pointwise estimate

**PBP2.** Let \(R=\operatorname{rad}(abc)\),
\(m=\min(a,b)\) and \(q=\min(P^+(a),P^+(b),P^+(c))\),
with \(P^+(1)=1\), for a primitive positive sum.
There is an absolute \(K>0\) such that for large \(R\),
\[
   Y:=\max\{\log(c/m),(\log c)/q\}
             \le \exp\!\left(K\sqrt{\log R\log\log R}\right).
\]
For \(Y\) large and \(t=K^{-2}(\log Y)^2\), this implies
\(\log R\ge t/\log t\). In particular, for any fixed
\(0<\beta\le1\),
\[
  q\le(\log c)^{1-\beta}\quad\text{or}\quad
  \log(c/m)\ge(\log c)^\beta
  \quad\Longrightarrow\quad
  \log R\ \ge\
  \left(\frac{\beta^2}{2K^2}-o(1)\right)
               \frac{(\log\log c)^2}{\log\log\log c}.
\]
Choose \(\eta=\log(c/m)/\log c\) in Theorem 1.4(1)
after interchanging the summands if needed.
Theorem 1.4(2) gives the other term in the maximum.
For \(x=\log R\), the two estimates yield
\(x\log x\ge t\). If \(t\) is large,
\((t/\log t)\log(t/\log t)\le t\), and
\(x\log x\) increases for \(x\ge1\).
The disjunction makes \(Y\ge(\log c)^\beta\).
Use this lower bound for \(t\) in the increasing
function \(t/\log t\) to obtain the displayed formula.
The proof of Pasten's Corollary 1.5 already uses the
small-\(q\) mechanism for a largest-prime bound.
The displayed version also prices the real ratio
and gives a radical cutoff in either sector.

For example, at a shrinking radical budget
\(\log R=o((\log\log c)^2/\log\log\log c)\),
the family evades both sectors for every fixed \(\beta>0\).
This is a pointwise restriction at that scale. At
fixed \(\lambda>0\), the bound \(\log R\le\lambda\log c\)
does not imply the required small-\(R\) hypothesis.
Thus neither this inversion nor the Stewart–Yu
pointwise estimate supplies the missing occupied-cell
inequality at fixed radical budget.
[Pasten, Theorem 1.4 and proof of Corollary 1.5](https://arxiv.org/html/2312.03566v1#S1)
· [Stewart–Yu, Theorem 1](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0012-7094-01-10815-6.pdf)

## Current fixed-budget source claims and arithmetic inputs

Hongmin Li's radical-budget paper 17, version 4,
states \(54/91\) on \(0<\lambda\le1\) in Theorem 1.1.
The exact constant \(5588/9417\) is recorded in its
Proposition 5.1 and the discussion of that theorem.
Theorem 1.2 gives the Farkas-envelope bound for
\(0.70\le\lambda\le1.10\), and Corollary 1.3(i)
states \(3/5\) up to \(247/245\). [Primary rStation v4](https://rstation.ai/problems/math/papers/17/v4)

The squares and twisted-energy paper 28, version 4,
states \(3622/6139\) at budget one in Theorem 1.1.
Theorem 1.2 gives its tabulated
\(589903/1000008\) at that budget as part of the
five-point certificate. Corollary 1.3 states the
piecewise-affine extension for \(0<\lambda\le2\).
These statements identify the source versions
and scope. The external certificate trees are
not used as proof premises here. [Primary rStation v4](https://rstation.ai/problems/math/papers/28/v4)

The source arithmetic mechanism used locally is
Zeng's fixed-degree weighted power-difference bound.
It counts after the numerical coefficients and
the degree are fixed. Browning–Verzobio's
coefficient-uniform generalized Fermat bound is a
separate input when each original letter has been
reduced to one exact powered root. The local valuation
grid transfer relies on that exact representation.
The hybrid splice transfer uses a different step:
the two original operand caps pay for primitive
outputs once the third radical is bounded, and
known common support restores the original divided
values. The new pair estimates show that an upper
bound of this kind on typical output radical would
have to supply fresh arithmetic. [Zeng, Lemma 3.1 and §§4–5](https://hexagonmath.org/pdf/2611.00621v1)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)
· [Exact-grid transfer](ABC_FINE_VALUATION_DIVISIBILITY.md#consequence-for-a-genuinely-fine-valuation-grid)
· [Physical shifts and direct splices](ABC_HYBRID_SAME_ROLE_OVERLAP_AND_FRESH_RADICALS.md#primitive-comparison-and-physical-shift-outputs)


### Short powerful windows and coprime lifts

Narumi–Tachiya give large occupancy between consecutive
powers by simultaneous approximation in prescribed
power classes. A local coprimality sieve and a
large-divisor tail from Roth give a positive density
of centres with any fixed number of pairwise coprime
one-sided occupants, all coprime to the centre.
A growing prime-power congruence can be imposed on the
actual short displacements with an explicit count cost.
These constructions saturate the limiting exponent
with **two** original radical caps. The entire
displacement is still outside the small-radical test.
[Narumi–Tachiya, Theorem 1 and Lemma 2](https://arxiv.org/html/2611.20656v1#S2)
· [Coprime short-power windows (CPW1–CPW2)](ABC_COPRIME_SHORT_POWER_WINDOWS.md)


### Quantitative Subspace Theorem on a moving divisor

On the complementary upper-bound side, use the actual
full prime support of the shift and apply
Evertse–Ferretti with the bounded list of algebraic
linear forms selected by the two exact powered letters.
The number of places can grow. Root-choice and
exponent-vector catalogues cost only an exponential in
the number of primes, which is height-subpower.
A cheap divisor of fixed positive numerical weight \(w\)
on the third letter gives a simpler finite-height count.
Fix the actual coefficient and divisor records and recover
the primitive slope as soon as \(kw>2\). The count thereby
allows any threshold \(K(L)\to\infty\), or a qualifying fixed
lower order, and the divisor may contain primes dividing \(k\).
The sector saving follows when the two records together
cost less than one third by a fixed margin. The Subspace
Theorem route above gives the separate unbounded-height
cardinality bound. Extracting affordable power representations
and divisors across arbitrary balanced profiles remains open.
[Evertse–Ferretti, Theorem 3.1](https://annals.math.princeton.edu/wp-content/uploads/annals-v177-n2-p04-p.pdf#page=11)
· [Elementary finite-height count (CPW5)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#an-adaptive-common-power-sector-with-a-third-letter-lock)
