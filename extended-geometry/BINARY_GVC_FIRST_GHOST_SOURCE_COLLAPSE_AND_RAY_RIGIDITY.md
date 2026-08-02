# Binary GVC carry promotion: first-ghost source collapse and ray rigidity

## 1. Scope

This note audits the proposed final local step in binary GVC carry promotion.
It proves two complementary facts:

1. common translation of the first Cartier ghost loses the decomposition of a
   low state and remembers only its total source exponent; and
2. the all-scale factorial sequence of a pure repeated ray determines its
   positive multiplicity partition.

The first result disproves the attempted single-shell beta/Bessel exposure
lemma. The second shows that this loss is repaired for pure rays by genuinely
all-order data. Neither result proves unrestricted GVC(2). The residual target
is a mixed, nonfree return-semigroup packet theorem.

## 2. First-ghost source-total collapse

Let `E` be a finite set of low channels. Give channel `e` an independent
coefficient variable `X_e` and a source exponent

\[
 \gamma_e\in\mathbb N^s.
\]

Put

\[
 A(z)=\sum_{e\in E}X_e(1+z)^{\gamma_e},
 \qquad
 (1+z)^\gamma=\prod_{i=1}^s(1+z_i)^{\gamma_i}.
\tag{2.1}
\]

For a prime `p`, let

\[
 G_p(A)=\frac{A^p-\Phi_p(A)}p
\tag{2.2}
\]

be the first Cartier ghost. For a non-pure state `n=(n_e)` with
`|n|=p`, define

\[
 \Gamma(n)=\sum_en_e\gamma_e.
\]

> **Theorem 2.1 (source-total collapse).**
> For every such state,
> \[
> [X^n]G_p(A(z))
> =\frac{(p-1)!}{\prod_en_e!}(1+z)^{\Gamma(n)}.
> \tag{2.3}
> \]
> Hence for every multivariate Hasse order `kappa`,
> \[
> [X^nz^\kappa]G_p(A(z))
> =\frac{(p-1)!}{\prod_en_e!}
>   \binom{\Gamma(n)}\kappa.
> \tag{2.4}
> \]
> After dividing each state column by its nonzero divided-power unit,
> every common-translation row depends only on `Gamma(n)`.

### Proof

The coefficient of `X^n` in `A^p` is

\[
 \binom pn(1+z)^{\Gamma(n)}
 =p\frac{(p-1)!}{\prod_en_e!}(1+z)^{\Gamma(n)}.
\]

The Frobenius lift contributes only to states supported on one channel, so it
has no `X^n` term. Division by `p` gives (2.3), and Hasse extraction gives
(2.4). The identity is integral; reduction modulo `p` does not cause the
collapse. \(\square\)

## 3. Uniform invisible Bessel atom

Use one source coordinate with levels `0,1,2,3`. For every prime `p>=5`, set

\[
 n_+=2e_2+(p-2)e_3,
 \qquad
 n_-=e_1+(p-1)e_3.
\tag{3.1}
\]

Then

\[
 |n_+|=|n_-|=p,
 \qquad
 \Gamma(n_+)=\Gamma(n_-)=3p-2,
\]

while

\[
 n_+-n_-=2e_2-e_1-e_3,
\tag{3.2}
\]

which is the adjacent centered-triple/Bessel lattice atom.

> **Corollary 3.1.** The normalized common-translation signatures of
> `n_+` and `n_-` agree at every Hasse order. Therefore the implication
> \[
> \text{first ghost plus further common translations}
> \Longrightarrow
> \text{isolated beta/Bessel atom}
> \]
> is false.

This is a counterexample to the proposed local exposure lemma, not to GVC(2).

## 4. Why the isolated Bell calculation was misleading

For two source exponents `gamma,gamma'`, the isolated two-block Bell sum is

\[
 \sum_{\substack{\mu+\nu=\kappa\\\mu,\nu\ne0}}
 \binom\gamma\mu\binom{\gamma'}\nu.
\tag{4.1}
\]

It does depend on the pair. But the actual diagonal coproduct also contains
the one-block endpoint terms. Multivariate Vandermonde gives

\[
 \binom{\gamma+\gamma'}\kappa
 =\binom\gamma\kappa+\binom{\gamma'}\kappa
 +\sum_{\substack{\mu+\nu=\kappa\\\mu,\nu\ne0}}
  \binom\gamma\mu\binom{\gamma'}\nu.
\tag{4.2}
\]

The complete row therefore remembers only the sum. Moment-cumulant Möbius
inversion is nonlinear and does not make the two-block summand an independent
linear observation. This is the exact reason earlier finite-span and
filler-alignment arguments appeared to close the gap but did not.

## 5. Bounded tangent window

For

\[
 B(Y)=b_0(1+U(Y)),\qquad b_0\ne0,
\]

and total Taylor degree below `p`, the Frobenius term has no contribution.
Using

\[
 \frac1p\binom pk\equiv\frac{(-1)^{k-1}}k\pmod p
 \qquad(1\le k<p),
\]

one obtains

\[
 b_0^{-p}G_p(B)\equiv\log(B/b_0)\pmod p
\tag{5.1}
\]

in every bounded degree `<p`. Thus the first ghost retains the tangent
logarithm before source convolution, but Theorem 2.1 shows that the common
source action has already identified all decompositions with the same total.

## 6. Scaled-factorial partition rigidity

For a tuple of positive integers `a=(a_1,...,a_r)`, put

\[
 F_a(k)=\prod_{i=1}^r(ka_i)!.
\]

> **Theorem 6.1 (scaled-factorial rigidity).** If
> \[
> F_a(k)=F_b(k)\qquad(k\ge1),
> \tag{6.1}
> \]
> then the positive integer multisets `a` and `b` are equal.

### Proof

Successive ratios give equality of the polynomials

\[
 \frac{F_a(k+1)}{F_a(k)}
 =\prod_i\prod_{u=1}^{a_i}(ka_i+u)
\tag{6.2}
\]

and its `b` analogue. Their roots are `-u/a_i`. If `M=max a_i`, the root
closest to zero is `-1/M`. It cannot arise from a smaller part, and its
multiplicity is exactly the number of parts equal to `M`. Remove those factors
and repeat. Induction recovers the entire multiset. \(\square\)

If the parts have common total `d`, the same conclusion follows from the
all-scale multinomial sequence

\[
 W_a(k)=\frac{(kd)!}{\prod_i(ka_i)!}.
\tag{6.3}
\]

At one scale factorial collisions do occur, for example

\[
 4!1!1!1!=3!2!2!,
\]

but the scaled sequence separates them.

## 7. Consequence for a pure Hall return ray

For a primitive balanced return state `(p,q)` of order `d`, the two
multinomial factors are

\[
 \binom{kd}{kp}\binom{kd}{kq}
 =\frac{((kd)!)^2}
 {\prod_i(kp_i)!\prod_j(kq_j)!}.
\tag{7.1}
\]

After fixing the oriented radial profile, Theorem 6.1 determines the combined
positive multiplicity multiset

\[
 \{p_i:p_i>0\}\mathbin{\uplus}\{q_j:q_j>0\}.
\tag{7.2}
\]

It does not, by itself, distinguish the operator and polynomial sides: a
cross-side redistribution can preserve (7.1). When marked-side data separates
the two sides, the theorem applies to each side separately. Otherwise that
redistribution must remain explicitly in the finite packet.

Thus a pure-ray collision cannot be explained by two genuinely different
integer partitions. After the existing marked-side and minimal-circuit
reductions, any unresolved cancellation must involve mixed sums of primitive
rays rather than only their pure repeats.

## 8. Revised binary-GVC frontier

The following proposed route must be retired:

\[
 \text{common translation of one first ghost}
 \Longrightarrow
 \text{isolated two-block beta/Bessel row}.
\]

The first ghost collapses source decompositions exactly. The all-scale
factorial sequence recovers pure-ray multiplicity data, but it does not by
itself classify nonfree mixed sums.

The remaining valid target is therefore:

> **Mixed-semigroup packet theorem.** A Hall-reduced equal-profile affine
> packet whose primitive rays are individually terminal cannot sustain
> all-order cancellation through nonfree mixed sums, unless it is a complete
> profile, has a split-symbol separator, or loses support.

The canonical conditional theorem still requires scale-compatible promotion.
This note narrows that problem and prevents reuse of a false local lemma; it
does not claim unrestricted GVC(2).
