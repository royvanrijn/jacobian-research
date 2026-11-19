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
[comparison.tex](../../unpublished/pointwise-binary-splitting-semigroups/comparison.tex).
Its proof is in
[graph.tex](../../unpublished/pointwise-binary-splitting-semigroups/graph.tex),
[fixed-shift-comparison.tex](../../unpublished/pointwise-binary-splitting-semigroups/fixed-shift-comparison.tex)
and [assembly.tex](../../unpublished/pointwise-binary-splitting-semigroups/assembly.tex).
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
| Fixed subgroup pairs, with fixed residue filters | [Subgroup model section](../../unpublished/pointwise-binary-splitting-semigroups/model.tex): identity-class invariance, projected rough means, direct record sieve and explicit model product. |
| Unfiltered self-sums for one fixed prime coset | [Prime-coset theorem](PRIME_COSET_SELF_SUMS.md): every negative/complex branch retained, positive leading domination, one allowed prime class, and finite-image criterion. |
| Prime factors three modulo four | [Two-branch check](PRIME_THREE_MOD_FOUR_SUMS.md), a specialization of the coset theorem. |
| Connors–Keating fixed shifts | [Indicator check](CONNORS_KEATING_INDICATOR_CORRELATIONS.md): unrestricted even inert valuations, dyadic states and a prime-square upper sieve, followed by the separate arithmetic model evaluation. |

Arbitrary pairs of prime prescriptions are not covered merely by
possessing Euler products. If their allowed prime sets are disjoint,
the common invariant bank class in (C2) can disappear. A replacement
requiring a different graph, a new energy estimate, or another prime
distribution theorem is outside this theorem.

The written theorem and application checks have internal mathematical
review. This does not assert external human review, formal verification,
historical priority, or independent revalidation of the complete imported
preprints. Finite character checks test normalization only, not the
infinite analytic conclusion.
