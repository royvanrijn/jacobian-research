# Fixed-dimensional Gaussian-to-Jacobian implication

This note proves the precise implication needed by the 79-variable
cubic-homogeneous map:

\[
 \mathrm{GMC}(2r)
 \Longrightarrow \mathrm{SIC}(r)
 \Longrightarrow
 \text{every cubic-homogeneous Keller map in }r\text{ variables is invertible}.
                                                               \tag{0.1}
\]

The proof follows the architecture of Derksen--van den Essen--Zhao and Zhao,
but all steps needed for one fixed `r` are written here, including the
countable-union lemma and the formal inversion identity.  Applying the
contrapositive at `r=79` makes the repository's route to `not GMC(158)` locally
self-contained.  It remains nonconstructive at the Gaussian-witness step.

## 1. Two contraction maps

Let

\[
 A_r=\mathbb C[w_1,\ldots,w_r,z_1,\ldots,z_r].
\]

Define

\[
 \mathcal E_r(w^\alpha q(z))=\partial^\alpha q(z),
 \qquad
 \mathcal F_r(h)=\mathcal E_r(h)|_{z=0}.            \tag{1.1}
\]

On monomials,

\[
 \mathcal F_r(w^\alpha z^\beta)
 =\delta_{\alpha,\beta}\,\alpha!.                  \tag{1.2}
\]

A linear subspace `M` of an algebra is a Mathieu--Zhao space if

\[
 p^m\in M\text{ for every }m\geq1
 \quad\Longrightarrow\quad
 qp^m\in M\text{ for all sufficiently large }m
\]

for every `q`.

The Special Image Conjecture `SIC(r)` is exactly the assertion that
`ker E_r` is a Mathieu--Zhao space.

## 2. Gaussian moments give `ker F_r`

Let `X_1,...,X_r,Y_1,...,Y_r` be independent standard real Gaussians and put

\[
 W_j=\frac{X_j-iY_j}{\sqrt2},
 \qquad
 Z_j=\frac{X_j+iY_j}{\sqrt2}.
\]

The change between `(X,Y)` and `(W,Z)` is invertible over `C`.  Circular
Gaussian contraction gives

\[
 \mathbb E(W^\alpha Z^\beta)
 =\delta_{\alpha,\beta}\alpha!,                    \tag{2.1}
\]

so for every `h in A_r`,

\[
 \mathbb E(h(W,Z))=\mathcal F_r(h).                 \tag{2.2}
\]

Assume `GMC(2r)`.  If `F_r(p^m)=0` for every `m>=1`, then the Gaussian
polynomial `p(W,Z)` has every pure moment zero.  For every `q`, GMC gives

\[
 \mathcal F_r(qp^m)=\mathbb E(q(W,Z)p(W,Z)^m)=0
\]

for all sufficiently large `m`.  Thus

\[
 \mathrm{GMC}(2r)\Longrightarrow
 \ker\mathcal F_r\text{ is a Mathieu--Zhao space}.  \tag{2.3}
\]

## 3. From `ker F_r` to `ker E_r`

For `a in C^r`, let `tau_a` translate the `z` variables:

\[
 (\tau_a h)(w,z)=h(w,z+a).
\]

On a monomial, the binomial theorem gives

\[
 \mathcal F_r(\tau_a h)=\mathcal E_r(h)(a).         \tag{3.1}
\]

We need one elementary geometric lemma.

### Lemma 3.1

Affine space over an uncountable field is not a countable union of proper
Zariski-closed subsets.

#### Proof

It is enough to show that for any countable collection of nonzero polynomials
`f_n in C[x_1,...,x_r]`, there is a point where every `f_n` is nonzero.

For `r=1`, each zero set is finite, so their countable union cannot exhaust
the uncountable field.  Induct on `r`.  Regard each `f_n` as a polynomial in
`x_r` and choose one nonzero coefficient `c_n(x_1,...,x_{r-1})`.  By induction,
choose `a'` with every `c_n(a')` nonzero.  Then every `f_n(a',x_r)` is a
nonzero univariate polynomial and has finitely many roots.  Choose `a_r`
outside their countable union.  QED.

### Proposition 3.2

If `ker F_r` is a Mathieu--Zhao space, then `ker E_r` is a Mathieu--Zhao space.

#### Proof

Suppose

\[
 \mathcal E_r(p^m)=0\qquad(m\geq1).                 \tag{3.2}
\]

By (3.1), for every `a`, all powers of `tau_a p` lie in `ker F_r`.  Fix `q`.
The Mathieu--Zhao property gives, for each `a`, an integer `N(a)` such that

\[
 \mathcal E_r(qp^m)(a)
 =\mathcal F_r(\tau_a(qp^m))=0
 \qquad(m\geq N(a)).                                \tag{3.3}
\]

For `N>=1`, define

\[
 Z_N=\{a:\mathcal E_r(qp^m)(a)=0\text{ for every }m\geq N\}.
\]

Each `Z_N` is Zariski closed, being an intersection of polynomial zero sets,
and (3.3) says

\[
 \mathbb C^r=\bigcup_{N\geq1}Z_N.
\]

Lemma 3.1 implies that some `Z_N` is all of `C^r`.  Hence
`E_r(qp^m)=0` as a polynomial for every `m>=N`.  This is precisely the
Mathieu--Zhao property for `ker E_r`.  QED.

Combining (2.3) and Proposition 3.2 proves

\[
 \mathrm{GMC}(2r)\Longrightarrow\mathrm{SIC}(r).    \tag{3.4}
\]

The countable-union step explains why the contrapositive does not construct a
small explicit Gaussian witness.

## 4. The formal inversion identity

For the formal identity, let

\[
 F(z)=z-H(z),
\]

where every component of `H` has order at least two.  There is a unique formal
inverse `G`; its homogeneous pieces are obtained recursively because the
linear part of `F` is the identity.  Write `j(F)=det DF`.

### Lemma 4.1 — Abhyankar--Gurjar identity

For every formal power series `g`,

\[
 g(G(z))=
 \sum_{\alpha\in\mathbb N^r}\frac1{\alpha!}
 \partial^\alpha\!\left(H(z)^\alpha g(z)j(F)(z)\right).            \tag{4.1}
\]

#### Proof

We compare the coefficient of `z^beta`.  Formal residue substitution under
`z=F(u)` gives

\[
 \begin{aligned}
 [z^\beta]g(G(z))
 &=\operatorname{Res}_z g(G(z))z^{-\beta-\mathbf1}\,dz\\
 &=\operatorname{Res}_u
   g(u)j(F)(u)\prod_{i=1}^rF_i(u)^{-\beta_i-1}\,du.                \tag{4.2}
 \end{aligned}
\]

For completeness, the formal residue change rule follows from a one-parameter
argument.  Put `F_s(u)=u-sH(u)`.  It is a formal automorphism over `C[[s]]`.
For a top form `omega=phi(z) dz_1...dz_r`, differentiation of `F_s^* omega`
is the pullback of the Lie derivative along the formal velocity field of
`F_s`.  Since every top form is closed, Cartan's formula makes this Lie
derivative an exact form.  The formal residue of an exact top form is zero:
in a term `partial_i A`, the coefficient that could contribute exponent `-1`
in variable `i` is multiplied by zero.  Hence

\[
 \operatorname{Res}_u F_s^*\omega
\]

is independent of `s`.  Its values at `s=0` and `s=1` give precisely the
change rule used in (4.2).

Since `F_i(u)=u_i-H_i(u)`, expand

\[
 F_i(u)^{-\beta_i-1}
 =u_i^{-\beta_i-1}
  \sum_{\alpha_i\geq0}
  \binom{\beta_i+\alpha_i}{\alpha_i}
  \left(\frac{H_i(u)}{u_i}\right)^{\alpha_i}.       \tag{4.3}
\]

Substitution in (4.2) gives

\[
 [z^\beta]g(G(z))
 =\sum_\alpha
  \binom{\beta+\alpha}{\alpha}
  [u^{\beta+\alpha}]H(u)^\alpha g(u)j(F)(u).        \tag{4.4}
\]

On the other hand,

\[
 \frac1{\alpha!}[z^\beta]\partial^\alpha A(z)
 =\binom{\beta+\alpha}{\alpha}[z^{\beta+\alpha}]A(z).            \tag{4.5}
\]

Equations (4.4)--(4.5), with `A=H^alpha g j(F)`, prove (4.1) coefficient by
coefficient.  The order assumption on `H` makes every coefficient a finite
sum.  QED.

Thus the differential-operator inversion step is part of the local proof,
not a black-box invocation.

## 5. `SIC(r)` forces the fixed cubic map to be invertible

After an invertible linear target normalization, a Keller map with identity
linear part and cubic-homogeneous nonlinear part has the preceding form.  Now
suppose `H` is cubic homogeneous and

\[
 j(F)=1.
\]

In `A_r`, set

\[
 p(w,z)=\sum_{i=1}^r w_iH_i(z).                    \tag{5.1}
\]

The multinomial theorem and (1.1) give

\[
 \frac1{m!}\mathcal E_r(p^m)
 =\sum_{|\alpha|=m}\frac1{\alpha!}
   \partial^\alpha(H^\alpha).                      \tag{5.2}
\]

Apply (4.1) with `g=1`.  Since `j(F)=1`,

\[
 1=\sum_{m\geq0}\frac1{m!}\mathcal E_r(p^m).       \tag{5.3}
\]

For cubic `H`, the `m`-th summand is homogeneous of degree `2m`.  Distinct
values of `m` have distinct degrees, so (5.3) implies

\[
 \mathcal E_r(p^m)=0\qquad(m\geq1).                 \tag{5.4}
\]

Assume `SIC(r)`.  Since `ker E_r` is a Mathieu--Zhao space, for every
coordinate `z_i`,

\[
 \mathcal E_r(z_i p^m)=0
\]

for all sufficiently large `m`.  Apply (4.1) with `g=z_i`:

\[
 G_i(z)=\sum_{m\geq0}\frac1{m!}\mathcal E_r(z_i p^m).             \tag{5.5}
\]

The `m`-th term has degree `2m+1`, and all but finitely many terms vanish.
Therefore every `G_i` is a polynomial.  The formal inverse of `F` is a
polynomial inverse, so `F` is invertible.  This proves the second implication
in (0.1).

## 6. Application to the 79-variable map

The [BCW reproduction](LONG_SU2_AND_BCW_REPRODUCTIONS.md) constructs

\[
 V(Z)=Z+\mathcal H(Z),\qquad Z\in\mathbb A^{79},
\]

where `mathcal H` is cubic homogeneous, `det DV=1`, and three distinct rational
points have one common image.  Put `H=-mathcal H`; then `V=Z-H` has exactly the
form used above and is not invertible.

If `GMC(158)` were true, (3.4) would give `SIC(79)`, and Section 5 would make
`V` invertible, contradicting its explicit collision.  Hence

\[
 \boxed{\neg\mathrm{GMC}(158)}.
\]

This is now a complete local proof of the route-based failure.  It is still
nonexplicit as a Gaussian counterexample because Proposition 3.2 uses the
uncountable-field countable-union argument.  Long's direct three-variable
witness remains a separate and much stronger explicit external result.

## 7. Reproduction

Run

```bash
python3 scripts/verify_fixed_gmc_sic_bridge.py
.venv/bin/python scripts/verify_long_bcw_79_route.py
python3 scripts/audit_long_bcw_79_independent.py
```

The first checker verifies the exact coefficient identities underlying
Sections 1, 3, and 4 in bounded multi-index ranges.  The latter two generate
and independently replay the complete
[79-variable sparse artifact](../artifacts/generated-results/long_bcw_79_counterexample.json).

Primary provenance remains with
[Derksen--van den Essen--Zhao](https://arxiv.org/abs/1506.05192) and
[Zhao](https://arxiv.org/abs/0902.0210); local reproduction is evidence, not
authorship or external review.
