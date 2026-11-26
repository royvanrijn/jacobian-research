# Growing families with private prime packets

A near-full primitive abc family contains a polylogarithmic number
of rows whose prime supports can be made disjoint after removing
arbitrarily small radical mass per row. The same conclusion holds
among the neighbors of one fixed complete letter at the half-radical
fibre exponent. One selection works against all subsequent
polynomial-height prime pools.

## Pricing empirical prime frequencies

Use \(n=X^L\) in the shell \(X/2<c\le X\), with \(L\) fixed
before the height limit. Let \(\mathcal F\) be any nonempty
family of primitive triples in the shell with
\(\operatorname{rad}(abc)<c^L\). For a uniform row let
\(\pi_p=\mathbb P_{\mathcal F}(p\mid abc)\). Price \(p\) at
\[
 w_X(p)=\min\bigl(\log_Xp,-\log_X\pi_p\bigr),
\]
using \(+\infty\) when \(\pi_p=0\). The partition sum in (EB3)
holds uniformly in the family, since
\[
 \sum_{p\le X}X^{-w_X(p)}
   \le\sum_{p\le X}\frac1p+\sum_{p\le X}\pi_p
   \ll \log\log X+\frac{\log X}{\log\log X}=o(\log X).
\]
The second sum is the expected number of prime divisors of
\(abc\le X^3\). For every fixed \(\chi,\eta>0\), at all
sufficiently small fixed \(L\) the number of rows with
\[
 \sum_{p\mid abc}\max\{0,\log(p\pi_p)\}\ge\chi L\log X
\]
is at most
\[
             O_{L,\chi,\eta,\varepsilon}
                    \bigl(X^{((1-\chi)/3+\eta)L+\varepsilon}\bigr).
                                                               \tag{GP1}
\]
For \(0\le\chi\le1\), their effective radical cap is
\((1-\chi)L\), so this is (EB4). For \(\chi>1\) the set is
empty. The prices can be formed anew on any subfamily.
[Effective prime prices and the one-third count (EB3–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

## Growing selections and uniform pool avoidance

Fix \(C>0\) and \(0<\xi<1\). There are \(\delta>0\) and
\(L_0>0\) such that for every fixed \(0<L<L_0\) and all
sufficiently large \(X\) the following conclusions hold.

- If \(|\mathcal F|\ge n^{1/3-\delta}\) under the hypotheses
  above, at least \(\tfrac12\lfloor(\log X)^C\rfloor\) distinct
  rows can be selected. For each selected row \(E\) and
  role \(W\in\{a,b,c\}\), there is a packet \(D_{E,W}\mid W\)
  consisting of whole prime powers of \(W\), with
  \[
         \log_n\operatorname{rad}(W/D_{E,W})\le\xi .
                                                               \tag{GP2}
  \]
  All packet supports across roles and selected rows are disjoint.
- If \(\mathcal F\) instead consists of primitive shell triples
  with one complete letter \(T\) fixed, the assertion holds on the
  two free roles provided
  \(|\mathcal F|\ge n^{1/6-\delta}\) and each free radical is at
  most \(n^{1/3+\delta}\), with \(\operatorname{rad}(T)\le X^{2L}\).
  The constants and height threshold are uniform in its value and role.

For either version, any positive integer \(P\le X^B\) and any \(\zeta>0\)
satisfy the simultaneous avoidance bound
\[
 \#\left\{E:
       \log_n\operatorname{rad}\gcd\left(P,\prod_WD_{E,W}\right)
                              \ge\zeta\right\}
                         \le\frac{B}{\zeta L}.                 \tag{GP3}
\]
Indeed the products of the packet radicals are pairwise coprime,
so the sum of all overlap masses is at most \(\log_n P\).
This includes \(P\) chosen after the rows or depending on
additional queries.

### High recurrence at the global and fixed-letter exponents

Choose \(a>C+2\). Define \(\pi_p\) using all eligible roles
in a row, meaning all three in the global case and just the
two free ones in the fixed-letter case. The pool
\[
              \mathcal P=\{p:\pi_p\ge(\log X)^{-a}\}
\]
has cardinality \(O((\log X)^{a+1})\) by the bound on
\(\sum_p\pi_p\). Put \(s=1-1/(2(a+1))\). Pricing these primes
at \(s\log_Xp\) and the others at \(\log_Xp\) satisfies (EB3):
the decreasing sum of \(p^{-s}\) over a set of that cardinality
is \(O((\log X)^{(a+1)(1-s)})=O((\log X)^{1/2})\).
In the global case (EB4) bounds the rows with total mass at
least \(\xi/2\) in this pool by
\(O_L(n^{1/3-(1-s)\xi/6+\eta+\varepsilon/L})\).

For the fixed-letter case, the Euler-product enumeration of
(EB3) gives the following uniform extension of the fixed-letter
count. Suppose the two free letters have raw radical caps
\(X^{O(L)}\) and effective caps \(X^{L\ell_U},X^{L\ell_V}\) in
the sense that their price sums are at most \(L\ell_U,L\ell_V\).
For prices satisfying (EB3), the count is
\[
 O_{L,\eta,\varepsilon,g}
       \bigl(X^{( \min(\ell_U,\ell_V)/2+\eta)L+\varepsilon}\bigr)
                                                               \tag{GP4}
\]
at sufficiently small fixed \(L\), uniformly in the fixed value
as long as the full raw radical budget is \(X^{O(L)}\). The
uniform constants and budget threshold depend on the fixed raw
cap bound.
To prove this, take the anatomic truncation and coefficient
costs smaller than the allowance. The fixed value has
divisor-many possible root tuples, which can be fixed at subpower
cost. Keep its actual raw height caps and set its dictionary costs
to zero. Bin the price of each exact free root with sufficiently
small mesh at the fixed truncation. The Euler-product count
bounds its dictionary exponent by the upper endpoint of that
bin, or by its raw exponent if smaller. Exact free roots on
each side have disjoint prime supports and their combined price
is at most the cap on the complete letter. Equation (GD3) with
singleton groups gives (GP4), including the truncation and bin
allowances. This uses a bound on the raw radical of the fixed
value only to justify the total raw root cap.

[Grouped dictionary count (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count) · [Fixed-letter half-radical bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

In the fixed-letter case, rows with combined pool mass at least \(\xi/2\) have at least \(\xi/4\) there in
one free letter. Its effective cap is at most
\(1/3+\delta-(1-s)\xi/4\). Equation (GP4), applied to each
choice of offending role, bounds their count by
\(O_L(n^{1/6+\delta/2-(1-s)\xi/8+\eta+\varepsilon/L})\).
In each of the two versions choose \(\delta,\eta\) small and
then \(\varepsilon/L\) small compared with \((1-s)\xi\).
The discarded rows are a power-smaller part of the family.
The small-budget threshold depends only on the fixed
parameters.

### Drawing the private packets

Sample \(m=\lfloor(\log X)^C\rfloor\) distinct rows uniformly
without replacement among those left after the recurrence
discard. For any chosen row and any of its primes outside
\(\mathcal P\), the conditional chance that another sampled
row uses that prime is \(O(m(\log X)^{-a})\), by a union bound
and the removal of only \(o(|\mathcal F|)\) rows. The expected
total normalized radical mass shared with another sampled
row, summed over the sample and outside \(\mathcal P\), is
\(O(m^2(\log X)^{-a})=o(m)\). Here the normalized mass per row
is bounded in both versions. Some sample has this bound.
Delete its \(o(m)\) rows where this loss exceeds \(\xi/2\),
and strip from each remaining row the entire prime powers on
\(\mathcal P\) and at every prime shared with another member
of the original sample. This proves (GP2) and (GP3). In the
fixed-letter version the packets and the sampling concern only
the free letters, which are coprime to the fixed value.

A variant gives \(m=\lfloor n^\kappa\rfloor\) in the global
family when only its primes \(p\ge n^b\) are to be preserved,
for fixed \(0<\kappa<a<b\) and \(\kappa<1/3-\delta\).
Use the recurrence threshold \(n^{-a}\). On this prime range,
\(\log(p\pi_p)\ge(1-a/b)\log p\) at a recurrent prime, so
(GP1) discards any prescribed positive combined mass in
these primes at a fixed power saving. The sampling loss
per row is then \(O(mn^{-a})=o(1)\). This version preserves
the designated prime-size range up to the stated loss and
makes no assumption on how much mass that range carries.

The earlier fixed-rectangle and translation-fibre arguments
control overlap in bounded configurations or against prescribed
pools. Here the selection grows with \(X\), and (GP3) is a
uniform bound for every later pool. A band of scaled degrees
at most \(S\) loses at most \(S\xi\) in numerical weight per
selected row. If each row retains positive numerical weight,
the LCM of the packets can have height \(X^{\Theta((\log X)^C)}\).
It cannot serve as a polynomial-height coefficient or pool in
existing pointwise counts. A fixed-letter selection gives a
large star, with no estimate here on the number of common
neighbors of its growing set of leaves.
[Fixed-rectangle overlap (FG2)](ABC_DIFFUSE_FACTOR_GRAPH.md#rectangles-with-controlled-prime-overlap) · [Fixed translation-fibre overlap (TC5)](ABC_ADDITIVE_TRANSLATION_CUBES.md#joint-radical-size-on-the-translated-forms) · [Exclusive packets on bounded tuples](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#exclusive-prime-packets-on-bounded-tuples)
