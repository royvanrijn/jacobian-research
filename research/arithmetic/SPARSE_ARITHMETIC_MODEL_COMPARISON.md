# Replacing sparse indicators by arithmetic models

This theorem isolates the shared analytic implication used by subgroup
prescriptions, fixed prime-coset self-sums and the fixed-shift
Connors–Keating argument. Its assumptions concern exact multiplicative
invariance, one-variable distribution and a positive sieve upper bound.
They do not assume centered binary cancellation, target-grid energy or a
leading pair asymptotic.

The inputs are precisely the existing ZF and G imports in
[the analytic input record](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md),
together with Walsh's elementary bilinear multiplier and
additive-to-Dirichlet conversion. No additional prime-distribution
theorem or graph construction is used.

The complete standalone statement is Theorem `thm:model-comparison` in
[comparison.tex](../../results/pointwise-binary-splitting-semigroups/comparison.tex).
Its proof is in
[graph.tex](../../results/pointwise-binary-splitting-semigroups/graph.tex),
[fixed-shift-comparison.tex](../../results/pointwise-binary-splitting-semigroups/fixed-shift-comparison.tex)
and [assembly.tex](../../results/pointwise-binary-splitting-semigroups/assembly.tex).
Those sources give the operator identities, exact normalizations and
precision order; the conditions and their independent proof roles are
recorded below.

## 1. Independently checkable hypotheses

Fix an even modulus \(L_0\), a reduced class \(c\bmod L_0\), and
fixed functions \(0\le f_1,f_2\le1\). Put
\[
 L=\log X,\quad y=\exp(L^{2/5}),\quad
 S=\{p:p\le y\},\quad E=\exp(\sqrt L).
\]
Let \(D_{i,X,K}\) be models indexed by a sufficiently large **fixed**
depth \(K\), and let \(\xi_i=f_i-D_{i,X,K}\). The entire model,
cutoff and small-prime record stay tied to the outer scale \(X\)
during every later dilation.

The two line families are:

* the reflected positive line \(Au+Bv=X\), where \(X\) is an integer
  target and \(A,B\) are fixed coprime positive \(L_0\)-supported
  coefficients;
* the fixed-shift line \((u,v)=(n,n+h)\), with fixed \(h\ne0\) and
  \(X\le n\le2X\).

Write \(t=Au\), \(J_X=[0,X]\) in the first family, and
\(t=n\), \(J_X=[X,2X]\) in the second. The displacement is
\(k=X\) or \(k=h\), respectively. Let a positive scale satisfy
\(XL^{-\kappa_0}\le\mathcal S_X\le XL^{\kappa_0}\), with
\(\kappa_0\) fixed.

**(C1) Boundedness.** For each fixed depth and sufficiently large \(X\),
\(0\le D_{i,X,K}\le1\) on \(1\le n\le X^3\). The bound is
independent of the depth; its threshold may depend on that depth.

**(C2) A common invariant prime class.** For every prime
\(p\le y\), \(p\equiv c\pmod{L_0}\), both endpoints satisfy
\[
 f_i(pn)=f_i(n),\qquad D_{i,X,K}(pn)=D_{i,X,K}(n),
\]
including repeated factors. Any residue filter is part of the function
and must satisfy the same invariance.

**(C3) Recordwise one-variable centering.** For every requested fixed
period/smoothness exponent \(B\) and accuracy \(D\), every sufficiently
large fixed depth gives
\[
 \left|\sum_{v\ S\text{-rough}}\xi_i(rv)a(v)W(rv/Y)\right|
      \ll_{B,D,K}(Y/r)L^{-D}.
\]
This is uniform for \(S\)-smooth records \(r\le E\),
\(X^{1/2}\le Y\le X^2\), \(|a|\le1\) periodic of period at most
\(L^B\), and compactly supported smooth rescaled weights with fixed
logarithmic derivatives bounded by fixed logarithmic powers.
It includes polynomial-logarithmic Mellin phases. It permits arbitrary
bounded weights on the **fixed record**, not arbitrary weights in the
rough variable.

**(C4) Positive record upper envelope.** Select disjoint small and big
prime banks in the ranges
\([\exp(L^{1/8}),\exp(L^{1/4})]\) and
\([\exp(L^{1/3}),\exp(L^{2/5})]\), in the class \(c\bmod L_0\),
deleting primes dividing \(k\). Each bank has reciprocal mass in
\([1/2,2]\), with at most \((\log L)/(16\varphi(L_0))\) banks of
either type. For their union \(\mathcal B\), use the independent law
\[
 \mu_p(0,0)=1-2/p,\qquad
 \mu_p(e,0)=\mu_p(0,e)=(1-1/p)p^{-e}\quad(e\ge1).
\]
There is no mass on two positive valuations. For every nonnegative
record function \(F\le L^B\), every interval \(J\subset J_X\),
and all four actual/model endpoint choices, require the upper bound
\[
 \sum_{t\in J}F(u_{\mathcal B},v_{\mathcal B})w_1(u)w_2(v)
 \ll (|J|/X)\mathcal S_X\,\mathbb E_\mu F
       +X^{1/2}L^{C_B}+XL^{C_B}e^{-c_0L^{1/10}}.
\]
This hypothesis is verified by stratum CRT root counts and an upper
sieve in the applications. Its endpoint/logarithmic exponents are
independent of the final fixed model depth. Constants and thresholds
may depend on that depth.

## 2. Shared comparison theorem

**Theorem SparseArithmeticModelComparison.** Under the existing scoped
ZF and G inputs, conditions (C1)–(C4) imply, for every sufficiently deep
fixed model, uniformly for intervals \(J\subset J_X\),
\[
 \sum_{t\in J}\{f_1(u)f_2(v)-D_{1,X,K}(u)D_{2,X,K}(v)\}
       \ll_K\mathcal S_X(\log X)^{-1/(8192\varphi(L_0))}.
\]
The reflected estimate is uniform in the target and its valuations.
The shift \(h\) in the second family stays fixed. There is no asserted
uniformity throughout intermediate growing sublinear shifts.

The model main term is not part of this theorem. Its evaluation is a
separate arithmetic calculation, including any secondary-branch error.

## 3. Why these assumptions imply cancellation

The proof uses mass-one banks in one reduced prime class. Characters
whose primitive quotient has conductor dividing \(L_0\) agree on a
product of a fixed number of bank labels up to the unit phase
\(\theta(c)^{\#\mathrm{slots}}\). Both products in each kernel have
the same slot type, so the character quotient frame reproduces the
required characters with these phases retained. Its Gram norm and signed
probe cancellation are unchanged.

G then supplies the physical root moment and endpoint pairing. Deleting
bank divisors of the displacement makes the residue change of variable
invertible at every bank prime. The shared-label identity leaves two
independent big-prime products. Invariance makes the coefficient left
after extraction of a unique small-bank prime independent of that prime.

For a reflected target, the major denominators are \(Xr\), with
\(r\) polynomial-logarithmic and physical arc width
\((\log X)^B/(XH_1H_2)\). In every gcd profile the cofactor/modulus
ratio is \(Y/(qP)\), which exceeds every logarithmic power. The same
ZF prime estimate controls nonexceptional characters; (C3) centers the
fixed-conductor characters after freezing the smooth record. This derives
target-grid energy. Discrete Parseval includes all grid points and
removes a spurious target factor. When endpoints are not units at the
fixed modulus, exceptional characters are identified by their primitive
conductors; their record factor is the unit phase \(\eta(A/d)\).

For fixed \(h\), the physical windows have width
\((\log X)^B/(|h|H_1H_2)\), which is wider. The proof does **not**
apply the reflected energy estimate at those widths. Walsh's existing
additive-to-Dirichlet conversion gives a low/high Mellin-frequency split.
(C3) controls low frequencies. The same ZF smooth-prime bound controls
high frequencies, while the cofactor length \(Y/P\) exceeds the
polynomial-height range \(X(\log X)^C\) by every logarithmic factor.
The ordinary Dirichlet mean-value estimate therefore pays the cofactor
energy. This is the separately displayed fixed-shift frequency argument.

Finally (C4) discards rare bank records at the correct sparse scale.
The exact Boolean reconstruction costs at most \(\log X\), and its
normalization introduces no divisor multiplicities. Shrinking cutoffs
and the same positive interval bound give the stated logarithmic rate.
Choose graph precision and comparison coefficients first, energy precision
next, and the fixed model depth last. Only constants and thresholds depend
on that last depth.

## 4. Applications and stopping boundary

| Application | Arithmetic checks |
|---|---|
| Fixed subgroup pairs, with fixed residue filters | [Subgroup model section](../../results/pointwise-binary-splitting-semigroups/model.tex): identity-class invariance, projected rough means, direct record sieve and explicit model product. |
| Unfiltered self-sums for one fixed prime coset | [Prime-coset theorem](PRIME_COSET_SELF_SUMS.md): every negative/complex branch retained, positive leading domination, one allowed prime class, and finite-image criterion. |
| Prime factors three modulo four | [Two-branch check](PRIME_THREE_MOD_FOUR_SUMS.md), a specialization of the coset theorem. |
| Connors–Keating fixed shifts | [Indicator check](CONNORS_KEATING_INDICATOR_CORRELATIONS.md): unrestricted even inert valuations, dyadic states and a prime-square upper sieve, followed by the separate arithmetic model evaluation. |

Arbitrary pairs of prime prescriptions are not covered merely by
possessing Euler products. If their allowed prime sets are disjoint,
the common invariant bank class in (C2) can disappear. A replacement
requiring a different graph, a new energy estimate, or another prime
distribution theorem is outside this theorem.

## 5. Classical replacement: a proved range and an open obligation

Walsh v2 explicitly anticipates sums-of-two-squares indicator correlations
in Section 1.2.4 (printed pp.5–6), citing Connors–Keating and Smilansky.
The general fixed-shift application idea is his announced direction.
The completed arithmetic models, retained character and valuation biases,
local factors, application checks and uniform reflected-target argument
are the contributions supplied here. Walsh also states on printed p.2
that his fixed-shift divisor theorem could use classical estimates in
place of the recent zero-free input. His written Lemma5.1 and
Section6.4 still use that input; the statement is not a supplied classical
replacement proof for this indicator comparison.

### 5.1. Trace the conductor and height separately

For fixed `h`, the shared graph comparison uses only classical
Siegel–Walfisz in its character frame, the combinatorial comparison and
the retained graph import. Its Fourier major arcs have
`q <= |h| L^B0`. Imposing the bank class changes the character modulus
only to `lcm(q,L0)`, still polynomial logarithmic. The fixed-shift proof
does not use the reflected target-grid energy estimate or its
target-sized conductors.

The arithmetic subgroup and quadratic-indicator checks of (C3) already
use a classical near-one contour. Positive envelopes, local model
correlations and fixed-class bank Mertens estimates are classical too.
The current complex prime-coset (C3) proof has a separate explicit ZF use;
that reflected application has not had its dependencies removed by this
fixed-shift audit.

The decisive remaining use is the smooth prime polynomial at Mellin
heights
\[
 T_0=xL^{B_0+B_1},\qquad
 e^{cL^{1/8}}\le P\le e^{CL^{1/4}},\qquad L=\log x.
\]
Small `q` does not make this a small-height estimate.

### 5.2. ClassicalModerateHeightPrimeEstimate

**Proved partial statement.** Fix `B,c,C>0`, set
`Z=exp(L^(1/16))`, and let `w` be smooth supported on `[1/2,2]`,
with every required fixed derivative bounded by a fixed power of `L`.
Uniformly for `q<=L^B`, characters `theta mod q`, `|t|<=2Z` and the
prime scales above,
\[
 \sum_p\theta(p)w(p/P)p^{it}
 =\delta_\theta\int_0^\infty
       \frac{w(u/P)u^{it}}{\log u}\,du
       +O\big(P L^{C_w}e^{-c_1L^{1/16}}\big).
 \tag{C.1}
\]
Constants and thresholds may be ineffective. This statement uses
classical estimates only, with no graph or recent zero-free input.

**Proof.** Use the classical Dirichlet zero-free region and Siegel's
theorem in its **all-character, ineffective** form. Exact statements and
their proofs are in Tao's
[254A Notes2](https://terrytao.wordpress.com/2014/12/09/254a-notes-2-complex-analytic-multiplicative-number-theory/),
Proposition51 and Theorem62(ii); Exercise45(i)–(ii) gives the local
zero-count/logarithmic-derivative formulas. These are classical inputs,
not consequences of Walsh's unsupplied alternative.

Put `delta=c2 L^(-1/16)`, with `c2` sufficiently small. For all the
characters in question, the rectangle
`Re s >= 1-2delta`, `|Im s| <= 3Z^2` is zero-free apart from the
allowed principal pole. Indeed `log(q(3Z^2+2))=O(L^(1/16))`.
A possible real exceptional zero is outside this rectangle too:
choose `epsilon<1/(32(B+1))`; Siegel gives
`1-beta >> q^(-epsilon) >> L^(-1/32)`, greater than `2delta` for
sufficiently large `x`. The latter threshold is ineffective.

On `Re s=1-delta`, regularize the principal logarithmic derivative by
subtracting its pole. There are `O(log(q(|Im s|+3)))` zeros in a fixed
local box, and every zero is at distance at least `delta` in the
relevant smaller rectangle. The local formula therefore bounds the
regularized logarithmic derivative by a fixed power of `L`.
The same holds on the horizontal edges used below.

Insert the von Mangoldt weight and the smooth test `w(u/P)/log u`.
Mellin inversion has weight translated in height by `t`; the Mellin
seminorms are polynomial logarithmic after this translation, independent
of the size of `t`. Truncate at imaginary part `Z^2` and move the
vertical line to `1-delta`. The pole produces the displayed principal
integral. The new vertical integral costs at most
`P^(1-delta) L^C`, and
`delta log P >= c3 L^(1/16)`. Fixed sufficiently many integrations by
parts in the smooth Mellin transform control both horizontal segments
and the discarded initial-line tail, because `|t|<=2Z` and the edge
height is `Z^2`. Their contribution is
`O(P L^C exp(-c4 L^(1/16)))`. Prime powers cost
`O(P^(1/2) L^C)` and are smaller. This proves (C.1).

Finite residue-class and fixed-bank restrictions are imposed by
orthogonality at the polynomial-logarithmic effective modulus.
Choosing `B` to cover that modulus only changes the fixed constants.
Integration by parts in the principal integral then supplies every
requested inverse power of `L` once `|t|` exceeds a sufficiently large
fixed power of `L`, throughout this moderate-height range. ∎

The factor two is a buffer: for outer `|t|<=Z`, every retained shift
`sigma+nu` is polynomial logarithmic, hence `|t+sigma+nu|<=2Z` for
sufficiently large `x`. Thus the partial lemma actually covers all the
prime arguments in this intermediate portion of the energy integral.

This proves a classical replacement for the **intermediate-height
piece**, not the fixed-shift correlation theorem without ZF. It does
not assert the strong error `P^sigma1` of the recent half-plane.

### 5.3. Open high-frequency obligation

The unsupplied classical step is on `Z<|t|<=T0`. With the actual
prime/cofactor polynomials from the manuscript, a sufficient target is
\[
 \int_{Z<|t|\le T_0}
   |P_r(t+\sigma+\nu)R_{r'}(t+\sigma+\nu)|^2\,dt
       \ll_D L^{-D},
 \tag{C.2}
\]
uniformly over the smooth dyadic pieces, residue classes, retained
Mellin parameters and **actual** cofactor coefficients, after paying
the fixed logarithmic partition costs. This is an explicitly open
weighted-energy estimate. It is not an admissibility hypothesis smuggled
into a new proved comparison theorem. A proof could also work directly
with the aggregate comparison rather than each individual factor.

At `t` of order `x`, the elementary classical zero-free width is of
order `L^(-1)`. With `log P<=CL^(1/4)`, its contour saving is only
`exp(-c log(P)/L)`, which does not tend to zero. Even a
Vinogradov–Korobov-width contour at the top height would not make
`(1-sigma)log P` grow at these bank scales. This disproves the proposed
**direct substitution**, not the possibility of a different classical
argument. The long-cofactor mean-value estimate bounds its total energy
by a polynomial logarithm; it does not show that this energy avoids
large values of the prime polynomial. Applying ordinary mean value
alone therefore does not prove (C.2).

An inspected alternative is Soundararajan,
[arXiv:1606.08021v1](https://arxiv.org/abs/1606.08021v1), Section4,
printed pp.13–15, using Lemmas2.5–2.6. His large-values bound gives a
small-prime exceptional set of size `x^o(1)` at these scales, but the
restricted mean-value bound for a general cofactor still has a constant
diagonal term. His next step extracts a second prime pool with
`exp(L^(4/5)) <= Q <= exp(L^(9/10))`, where classical prime cancellation
is available, and reuses restricted-set mean value. That pool lies
above our `y=exp(L^(2/5))`. Even when an actual indicator has invariance
under those split primes, its present density model does not satisfy
the same exact invariance there. A suitable second-pool decomposition
for the actual centered model, the omitted-integer error, the required
logarithmic precision and the preservation of low-frequency centering
have not been proved here. This is a concrete candidate route, not a
completed application of his theorem or a new imported hypothesis.

A new classical large-values/weighted-energy argument would need to
establish this missing step and preserve the actual record coefficients
and parameter order. The existing fixed-shift, Connors–Keating and
Smilansky entries consequently retain ZF and G. No reflected dependency
has been removed, and a fixed-shift success would still leave the
growing-target conductor issue there.

The written theorem and application checks have internal mathematical
review. This does not assert external human review, formal verification,
historical priority, or independent revalidation of the complete imported
preprints. Finite character checks test normalization only, not the
infinite analytic conclusion.
