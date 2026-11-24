# Inverse powers: exact reductions and the remaining contact theorem

Internal research note, 21 November 2026. The general inverse-power
contact bound in Section 1 is **not proved here**. The proved results are
the field-degree statement, the exact cyclic-descent criterion and its
contact corollary, the unicritical case, and the counterexample to one
proposed critical-value budget. None closes unrestricted SFC or the
primitive monomial-linear sharp contact problem. Existing canonical
sources and status records are unchanged.

## 1. The precise target and its rational block-average consequence

Let `A(z)=z+O(z^2)` be a polynomial of degree `d>=2`, let `V=A^{-1}`,
and let `c` be the number of distinct finite **nonzero** critical values
of `A`. For integers `b>=1` and `N>=b`, the proposed bound is

\[
 \operatorname{ord}_0(V^b-Q)\le N+bc
 \qquad (\deg Q\le N).
 \tag{1}
\]

Substitution by `A` preserves order, so (1) is equivalent to the
corresponding bound for `z^b-Q(A(z))`. This equivalence is a reformulation,
not a proof. The retained nonzero-critical-value contact theorem proves
the case `b=1`.

Here is the exact factorial implication. Fix positive integers `a,b,v_j`
and distinct nonzero coefficient values `c_j`, with `j=1,...,h`. Use
`p_j=bv_j` ordinary variables in coefficient group `j`, and allocate a
total common exponent `R_j=av_j` among those variables in any nonnegative
integer way. Set

\[
 R(z)=\prod_{j=1}^h(1-c_jz)^{v_j},\qquad
 g_m=[z^m]R(z)^{-(am+b)},\qquad A_*(z)=zR(z^a).
 \tag{2}
\]

Direct multinomial expansion of the ordinary factorial functional gives
`g_m=Lambda(F^m)/(m! prod_i(mr_i)!)`. Lagrange inversion gives, for the
normalized inverse `V_*` of `A_*`,

\[
 [t^{am+b}]V_*(t)^b=\frac{b}{am+b}g_m.
 \tag{3}
\]

The remaining degree classes of `V_*^b` vanish. To count critical values,
write `C(u)=prod_j(1-c_j u)`. Then

\[
 A_*'(z)=\frac{R(z^a)}{C(z^a)}
 \left(C(z^a)+a z^a C(z^a)\frac{R'(z^a)}{R(z^a)}\right).
 \tag{4}
\]

The bracket has degree at most `h` in `z^a`. The other derivative roots
map to zero. Thus `#CV^times(A_*)<=ah`.

If `g_n,...,g_(n+bh-1)` vanished, truncate `V_*^b` at
`N=a(n-1)+b`. The next possible nonzero term would have degree at least
`N+a(bh+1)>N+b #CV^times(A_*)`, contradicting (1). Therefore (1) would
give a `bh`-moment window, no longer than the actual ordinary support
`b sum_j v_j`. This implication includes every starting order and any
allocation of the grouped exponents. It does **not** include arbitrary
unequal group-average exponents.

## 2. Taking a power does not lower the algebraic degree

**Proposition 1.** For every normalized polynomial `A` as above and every
integer `b>=1`,

\[
 \mathbb C(A(z),z^b)=\mathbb C(z).
 \tag{5}
\]

Consequently the algebraic function `V(t)^b` has degree exactly `d`
over `C(t)`. Its irreducible polynomial has degree exactly `b` in `t`
and exactly `d` in the function variable.

**Proof.** The extension `C(z)/C(z^b)` is Galois with automorphisms
`z -> omega z`, where `omega^b=1`. An automorphism fixing `A(z)` must
satisfy `A(omega z)=A(z)`. Comparing the linear terms gives `omega=1`.
The intermediate field in (5) therefore has trivial fixing subgroup.

In particular `[C(A,z^b):C(A)]=d`. Alternatively form

\[
 \mathcal P(t,y)=\prod_{\omega^b=1}(t-A(\omega z)),\qquad y=z^b.
 \tag{6}
\]

The product is a polynomial in `t,y`. Its degree in `t` is `b`, and its
degree in `y` is `d`: the leading product of the degree-`d` terms in
`A(omega z)` is a nonzero multiple of `y^d`. The same Galois orbit
argument makes (6) irreducible over `C(y)`, hence absolutely irreducible
as a polynomial. Finally `t=A(z),y=z^b` satisfies it. QED.

Thus taking an inverse power does not algebraically reduce the covering
rank to `b`. An order-`b` scalar differential system for the distinguished
branch would need a separate proof; algebraic degree is not the same
invariant as minimal differential order. No prohibition on all possible
differential or recurrence methods is asserted.

### What the ordinary trace actually says

Let `z_1(t),...,z_d(t)` denote all roots of `A(z)=t`. Newton identities
give

\[
 \sum_{i=1}^d z_i(t)^b\in\mathbb C[t],\qquad
 \deg_t\left(\sum_i z_i(t)^b\right)\le\lfloor b/d\rfloor.
 \tag{7}
\]

For indices below `d`, the Newton identities do not use the constant
coefficient `-t`; subsequently each use of `t` lowers the power index by
`d`. This proves the degree bound by induction. The trace is a sum over
all branches, whereas (3) uses the branch through zero.

For example, `A=z-z^3` has `sum_i z_i(t)^2=2`, while its distinguished
branch has

\[
 V(t)^2=t^2+2t^4+7t^6+30t^8+\cdots.
 \tag{8}
\]

The zero higher coefficients of the trace provide no zero-window
propagation for this individual branch. The first two resultant
identities are also explicit:

\[
 \begin{array}{ll}
 A=z-z^2,&t^2+2ty+y^2-y=0,\quad y=V^2;\\
 A=z-z^3,&t^2-y^3+2y^2-y=0,\quad y=V^2.
 \end{array}
 \tag{9}
\]

## 3. Exact criterion for descent to an ordinary inverse

**Proposition 2.** There exists a normalized polynomial `B` such that

\[
                     V(t)^b=B^{-1}(t^b)
 \tag{10}
\]

if and only if `A(z)=zR(z^b)` for a polynomial `R` with `R(0)=1`.
When it exists, it is uniquely

\[
                         B(u)=uR(u)^b.
 \tag{11}
\]

**Proof.** Substitution by `A` turns (10) into
`B(z^b)=A(z)^b`. The right side must therefore be invariant under
`z -> omega z`. The rational function `A(omega z)/A(z)` has constant
`b`th power one, so it is a constant root of unity. Its limit at zero is
`omega`. Hence `A(omega z)=omega A(z)` for every `b`th root of unity.
Precisely the degrees congruent to one modulo `b` can occur in `A`.
This gives the stated form and (11). Conversely (11) immediately gives
`B(z^b)=A(z)^b` and, by formal inversion, (10). QED.

**Contact corollary.** In this compatible cyclic class one has the
stronger all-order estimate

\[
                      \operatorname{ord}_0(V^b-Q)\le N+c.
 \tag{12}
\]

Indeed the nonzero critical values of `A` occur in free orbits under
multiplication by `b`th roots of unity, and raising them to the `b`th
power identifies precisely those orbits with the nonzero critical
values of `B`. To check this directly, away from `R=0` the two
derivatives vanish exactly when `R(u)+buR'(u)=0`, where `u=z^b`;
the discarded roots map to zero. Therefore `c=b #CV^times(B)`.
Apply the retained inverse gap theorem to `B^{-1}` and then substitute
`t^b`. The largest contact past a truncation at degree `N` is at most
`b floor(N/b)+b #CV^times(B)<=N+c`. If `Q` is not that truncation,
its contact is already at most `N`.

This is a valid reduction to the known theorem, not a proof of (1) in
general. For (2), after reducing `gcd(a,b)=1`, the criterion requires
`R(u)=S(u^b)`. Generic coefficient data do not have this extra cyclic
symmetry. For the smallest example `A=z-z^2,b=2`,
`A(z)^2=z^2-2z^3+z^4` is not a polynomial in `z^2`; (10) is unavailable
even though that example itself satisfies (1).

## 4. A separate proved case: one critical point

If `A` has exactly one distinct critical point, normalization forces

\[
 A(z)=\frac{(1+\lambda z)^d-1}{d\lambda},\quad\lambda\ne0,
 \qquad
 V(t)=\frac{(1+d\lambda t)^{1/d}-1}{\lambda}.
 \tag{13}
\]

For every `n>=1`, the coefficient of `t^n` in `V` is a strictly positive
rational multiple of `(-1)^(n-1) lambda^(n-1)`. Consequently every
coefficient of `V^b` at degree `n>=b` is a strictly positive rational
multiple of `(-1)^(n-b) lambda^(n-b)`: all summands in its convolution
have the same nonzero factor and positive remaining factors. Thus these
coefficients never vanish, for every `b>=1`. In particular (1) holds.

Having one **nonzero critical value** does not imply having one critical
point; the latter hypothesis in this paragraph is essential to the
proof. Polynomials with additional critical points over zero, or several
critical points sharing a nonzero value, are not covered by (13).

## 5. A genuine truncation defeats the divided-polynomial budget

A tempting construction is
`Psi(z)=Q(A(z))/z^(b-1)`. If `Q` starts with `t^b`, this is a normalized
polynomial and

\[
 \operatorname{ord}(\Psi-z)
   =\operatorname{ord}(V^b-Q)-(b-1).
 \tag{14}
\]

The ordinary critical-value theorem would give (1) if one could always
bound

\[
                 \#\operatorname{CV}^{\times}(\Psi)\le N+bc-b.
 \tag{15}
\]

But (15) is false even for an actual inverse-power truncation. Take

\[
 A=z-z^2,\quad b=2,\quad c=1,\quad N=4,\qquad
 Q(t)=t^2+2t^3+5t^4.
 \tag{16}
\]

Here `Q` is exactly the degree-four truncation of `V^2`, and

\[
 \Psi=z-14z^4+28z^5-20z^6+5z^7.
 \tag{17}
\]

Its exact critical-value resultant is

\[
 \operatorname{Res}_z(\Psi'(z),\Psi(z)-w)=25w\,C(w),
 \tag{18}
\]

where

\[
 C(w)=2573571875w^5-94035000w^4-54986400w^3
       -21376768w^2-225280w-29696.
 \tag{19}
\]

It has nonzero constant term and is squarefree. For example, its
reduction modulo `101` is
`-35w^5+40w^4+20w^3-17w^2-50w-2`, coprime to its derivative.
Thus `Psi` has **five** distinct nonzero critical values, whereas the
right side of (15) is `4`. Division by `z^(b-1)` has destroyed the
required critical-value budget.

This is not a counterexample to (1): the actual contact in (16) is five,
and (1) permits six. It refutes this proposed proof mechanism.

## 6. Primary-source comparison and what remains open

Alin Bostan's [2016 AEC course](https://mathexp.eu/bostan/AEC16/aec16-p3.pdf),
slide 22/35 (PDF page 27), explicitly labels as a conjecture a general
bound of `deg_t(P)(deg_y(P)-1)-1` consecutive zero coefficients for a
nonrational algebraic series in `Q[[t]]` with minimal polynomial
`P(t,y)`. The
actual two different degree subscripts were inspected; the search
snippet incorrectly duplicates one of them. By Proposition 1 its
degree-only specialization to `V^b` for rational-coefficient `A` would
allow `b(d-1)-1` consecutive zeros. Equation (1), with `c<=d-1`, is at
least as strong even on this rational subclass, and also asks for
arbitrary complex coefficients. This comparison does not assert the
present status of the general conjecture or import it as a theorem.

Seth Dutter, [Gaps in Taylor series of algebraic functions](https://doi.org/10.1017/S0004972715000040),
*Bull. Aust. Math. Soc.* **91** (2015), 412--418, Theorem 2.1 and its
proof on printed pages 414--417, were inspected directly. Its estimate
bounds the index of the `n`th nonzero term by a height plus `n-1` times
an explicit ramification/pole quantity. It controls density of nonzero
terms; it does not give (1) for each consecutive block. Its proof's
auxiliary derivative combination and pole count do not supply the
missing sharp constant here.

The exact missing statement for rational equal block averages remains
(1), or the weaker specialization only to `A_*=zR(z^a)` in (2). The
covering-degree and divided-polynomial obstructions above reject two
shortcuts. They do not refute either contact target. The proven cyclic
and unicritical subclasses do not cover arbitrary unequal primitive
exponents in the original monomial-linear observable.

The proved algebraic identities and subclasses hold over arbitrary
characteristic-zero fields. For the analytic inverse-contact input,
embed the field generated by the finitely many coefficients into
`C`; vanishing and nonvanishing of every specified finite coefficient
block are preserved. This descent does not promote the unproved
equation (1).

## Exact replay

Run from the repository root:

```text
python research/factorial/image-rigidity/verify_sfc_inverse_power_contact.py
```

The standard-library checker verifies the displayed polynomial
identities, two compatible cyclic examples, the incompatible quadratic
control, the unicritical coefficient formula on stated finite inputs,
and the complete resultant identity (18) by its degree bound and seven
exact determinant evaluations. It checks squarefreeness modulo `101`
and rejects altered and malformed critical-value certificates. The
unrestricted claims are proved above; finite controls are not offered
as evidence for (1).

Retained dependencies:
[POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md](POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md)
for the nonzero-critical-value inverse theorem, and
[WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md](WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md)
for comparison with the already closed integer-average family.
