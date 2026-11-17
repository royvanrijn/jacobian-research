# Arithmetic indistinguishability without Keller equivalence

This note separates two attainable statements suggested by the universal
incidence cover.  The first is an exact degree-five construction: stably
inequivalent weighted Keller maps have identical fiber schemes after a
constructible relabeling of every finite-field target.  The second is the
genuinely Gassmann version: degree-seven Davenport polynomials give
nonisomorphic, arithmetically equivalent vertical fibers after seed
normalization and Hilbert specialization.

The later
[global Sunada theorem](GLOBAL_SUNADA_KELLER_COVERS.md) retains the Davenport
parameter `T`, computes the common branch surface, and simultaneously lifts
the two covers to stably inequivalent relative weighted Keller maps after a
finite tangent-mark base change.  The fixed-`T` treatment below remains the
number-field specialization of that result.

The distinction matters.  The degree-five incidence covers below are
isomorphic after an affine change of the pencil target, even though the
weighted Keller maps are not.  The degree-seven covers are not isomorphic;
their equality of arithmetic data comes from nonconjugate Gassmann-equivalent
subgroups.

## 1. A no-go theorem for the generic two-parameter pencil

For every degree-`N` seed `H`, the cover

\[
 H(W)-sW+t=0
\]

over `A^2_(s,t)` has geometric and arithmetic monodromy `S_N` in its natural
action.  Hence a nontrivial Gassmann pair cannot live over the full generic
two-parameter base.

Indeed, equality of two degree-`N` permutation characters in a common Galois
closure first gives equality of their kernels: an element is in the kernel
exactly when its character value is `N`.  Both actions therefore factor
faithfully through the same `S_N` quotient.  A connected degree-`N` natural
`S_N` cover corresponds to an `S_(N-1)` point stabilizer, and all such
stabilizers are conjugate.  The covers are then equivalent.

Thus a genuine Gassmann construction must force monodromy to drop on a proper
specialization locus.  A vertical polynomial slice is the natural place to
do this.

## 2. Exact degree-five arithmetic twins

Recall

\[
 H_\lambda(W)=
 \frac{W^2(W-1)(3W^2-(5\lambda+1)W+3\lambda)}{60},
 \qquad c_\lambda=\frac{\lambda-1}{30}.
\]

Put

\[
 \mu=\frac45-\lambda,\qquad \beta=\frac45,
\]

and

\[
 u_\lambda=\frac{2(5\lambda-2)}{1875},\qquad
 v_\lambda=-\frac{4(25\lambda+8)}{46875}.
\]

Direct expansion gives the key identity

\[
 \boxed{H_\mu(\beta-W)=-H_\lambda(W)+u_\lambda W+v_\lambda.} \tag{2.1}
\]

Consequently, with

\[
 s'=s-u_\lambda,
 \qquad
 t'=-t-v_\lambda+\beta(s-u_\lambda),                 \tag{2.2}
\]

one has

\[
 H_\mu(\beta-W)-s'(\beta-W)+t'
 =-(H_\lambda(W)-sW+t).                              \tag{2.3}
\]

This identifies the complete root schemes, including residue degrees and
multiplicities, over every field on which the displayed constants are
defined.  It is much stronger than equality of the leading Chebotarev law.

### Transfer to the open Keller target

For the weighted maps `F_lambda` and `F_mu`, write

\[
 s=BC,\qquad t=c_\lambda AC^2.
\]

On `C!=0`, equations (2.2) give the target relabeling

\[
 C'=C,
 \qquad B'=B-\frac{u_\lambda}{C},
\]

\[
 A'=\frac{-c_\lambda AC^2-v_\lambda+
                 \beta(BC-u_\lambda)}{c_\mu C^2}.    \tag{2.4}
\]

Simple roots reconstruct the affine source points and repeated roots are
reconstruction poles.  Therefore (2.3)--(2.4) identify every affine Keller
fiber over `C!=0`, as a reduced zero-dimensional scheme, after target
relabeling.

The pole in (2.4) is essential.  This pencil equivalence does not extend to a
polynomial target automorphism across `C=0`.

### The `C=0` fibers

Write `h_2(lambda)=-lambda/20` and

\[
 k_\lambda=\frac{h_2(\lambda)}{c_\lambda}
 =-\frac{3\lambda}{2(\lambda-1)}.
\]

The direct boundary calculation identifies the fiber over `(A,B,0)` with a
disjoint rational root-one sheet and

\[
 D_\lambda(A,B)x^2=1,
 \qquad
 D_\lambda(A,B)=\frac{B^2}{c_\lambda^2}-4k_\lambda A. \tag{2.5}
\]

The linear relabeling

\[
 (A,B,0)\longmapsto
 \left(\frac{k_\lambda}{k_\mu}A,
       \frac{c_\mu}{c_\lambda}B,0\right)             \tag{2.6}
\]

preserves `D` exactly.  Hence it identifies the boundary fiber schemes and
their zeta functions over all finite extensions.

For an odd residue field of size `q`, this also gives the exact boundary
histogram

\[
 B_1(q)=\frac{q(q+1)}2,\qquad
 B_3(q)=\frac{q(q-1)}2,\qquad B_j(q)=0\ (j\ne1,3),    \tag{2.7}
\]

with the nonsquare values of `D` refining the singleton bin by one closed
point of residue degree two.

Combining (2.4) on `C!=0` with (2.6) on `C=0` gives a constructible target
bijection.  At every good finite place, and over every finite residue-field
extension, corresponding fibers have the same rational-point count and the
same residue-degree data.  In particular, the complete target histograms of
`F_lambda` and `F_mu` agree exactly, including the lower-order and boundary
terms invisible to the universal `S_5` density.

### Why the Keller maps are nevertheless inequivalent

The unmarked Hessian triples are exchanged by `W -> 4/5-W`, but this affine
map moves the distinguished zero/root-one boundary marking.  On the
affine-marked clean open, exact seed recovery says that stable polynomial
left--right equivalence of `F_lambda` and `F_mu` forces `lambda=mu`.  Since

\[
 \lambda=\mu\quad\Longleftrightarrow\quad\lambda=\frac25,
\]

every parameter in the dense open

\[
 \{\lambda:\lambda,4/5-\lambda\in\Lambda^\circ,
                  \ \lambda\ne2/5\}
\]

gives stably inequivalent arithmetic twins.  This is already a Keller-map
analogue of arithmetic indistinguishability, but not yet an analogue of
nonisomorphic arithmetically equivalent fields: the open incidence fibers in
this construction are actually isomorphic.

## 3. The genuine Gassmann specialization in degree seven

Let `K=Q(a)`, where

\[
 a^2+a+2=0,
\]

so `K=Q(sqrt(-7))`, and let the bar denote the nontrivial automorphism of
`K`.  Cassou-Nogues and Couveignes give the one-parameter polynomial

\[
\begin{aligned}
 g_T(Y)={}&\frac17Y^7+(1+a)TY^5+(1+a)TY^4
 -(3-2a)T^2Y^3-2(1-2a)T^2Y^2\\
 &-\frac1{28}(5+3a)(28T-2-11a)T^2Y-(1+a)T^3,
\end{aligned}                                                   \tag{3.1}
\]

and `h_T=bar(g_T)`.  The cubic

\[
\begin{aligned}
 Q_T(Y,Z)={}&Y^3+aY^2Z+Y^2Z+aYZ^2-Z^3\\
 &+(5+3a)TY+(-2+3a)TZ+(2a+1)T
\end{aligned}                                                   \tag{3.2}
\]

divides `g_T(Y)-h_T(Z)`.

The common Galois closure has group

\[
 \Gamma=GL_3(F_2)\simeq PSL_2(7),
\]

of order `168`.  Its two degree-seven actions are on the points and the lines
of the Fano plane.  If `U` is a point stabilizer and `V` a line stabilizer,
then

\[
 U\not\sim_\Gamma V,
 \qquad
 \mathbf 1_U^\Gamma=\mathbf 1_V^\Gamma.             \tag{3.3}
\]

Thus the actions are nonisomorphic but have the same permutation character.
The three finite branch cycles can be taken as

\[
 (12)(36),\qquad(23)(45),\qquad(34)(67),             \tag{3.4}
\]

and their product is a seven-cycle at infinity.  This explains precisely how
the specialization evades the generic-pencil no-go theorem: three pairs of
simple critical points share three critical values, so local inertia is a
double transposition rather than a transposition.  The generated group is
`GL_3(F_2)`, not `S_7`.

Because character equality holds also on every power of an element, (3.3)
determines equal cycle types, not only equal fixed-point counts.  The two
covers therefore have identical unramified splitting statistics.

### Normalization into weighted seeds

The passage from (3.1) to the repository seed coordinates is elementary.
Over a finite constant extension, choose `r!=0` satisfying

\[
 g_T(r)-g_T(0)-r g_T'(0)=0                         \tag{3.5}
\]

and put, up to a nonzero scale,

\[
 H(W)=g_T(rW)-g_T(0)-r g_T'(0)W.                    \tag{3.6}
\]

Then

\[
 H(0)=H'(0)=H(1)=0,
\]

and `c_H=-H'(1)` supplies the weighted normalization.  Apply conjugation to
`r` and (3.6) to obtain the seed `G` from `h_T`, and work over a common finite
extension containing both tangent marks.  The generic Galois closure is
regular over the coefficient field, so this finite constant extension
preserves `Gamma`.  Removing the usual proper closed degeneracy loci gives
admissible, Hessian-clean, boundary-clean degree-seven seeds.

The vertical line

\[
 s_H=-r g_T'(0)
\]

in the `H` pencil is, after scaling and translating `t`, exactly the cover
`g_T(Y)-z=0`; similarly the corresponding vertical line for `G` is
`h_T(Z)-z=0`.  Taking `C=1`, these become explicit one-parameter target lines
in the two Keller maps.

For `T!=0`, the Hessian divisors of `g_T` and `h_T` are not affinely
equivalent.  Both second derivatives have leading coefficient `6` and zero
degree-four coefficient, so any affine equivalence has zero translation.  If
its scaling is `alpha`, comparison of the degree-three and degree-two
coefficients forces

\[
 \alpha^2=-\frac{1+a}{a},
 \qquad
 \alpha^3=-\frac{1+a}{a}.
\]

This would give `alpha=1` and `2a+1=0`, contradicting `a^2+a+2=0`.
Therefore the normalized seeds retain different Hessian divisors.  On the
clean locus, the Hessian obstruction already proves that the resulting
weighted Keller maps are stably polynomially left--right inequivalent.

### Hilbert specialization

Choose a nondegenerate algebraic value of `T` and compatible tangent marks,
then apply Hilbert irreducibility in `z` to the common regular
`Gamma`-closure of `g_T(Y)-z` and `h_T(Z)-z`.  Outside a thin set, the
specialized Galois group remains `Gamma`.  The two degree-seven fields are
then fixed fields of `U` and `V`; (3.3) gives equal Dedekind zeta functions,
while nonconjugacy gives nonisomorphic fields.

This yields the strong target:

> Two stably inequivalent degree-seven weighted Keller maps over a number
> field, with different Hessian divisors, contain target lines whose Hilbert
> specializations are nonisomorphic arithmetically equivalent fields.

The Davenport family starts over `Q(sqrt(-7))`, and seed normalization may
adjoin the finite tangent-mark field.  It does not start over `Q`: the
branch-cycle field of definition does not descend there.  A rational version
would require a different Gassmann realization or restriction of scalars,
and should be treated as a separate problem.

## 4. What is proved and what remains computational

The checker verifies:

1. the degree-five pencil identity and boundary discriminant matching;
2. the explicit divisor (3.2) of `g_T(Y)-h_T(Z)`;
3. non-equivalence of the two Hessian divisors for `T!=0`;
4. equality of the point/line permutation characters of `GL_3(F_2)`;
5. nonconjugacy of the two stabilizers; and
6. generation of the order-`168` group by (3.4).

Run

```bash
.venv/bin/python scripts/verify_arithmetic_indistinguishability.py
```

The remaining useful computational milestone is a small integral numerical
specialization.  Choose `T` and a tangent root `r`, construct `H,G` by (3.6),
certify all weighted clean-open conditions, and search a Hilbert parameter
`z` whose specialized group has order `168`.  Factorization at a few good
primes can certify the required Frobenius classes; a resolvent or direct
Galois-group computation then completes the certificate.

For the global cover, rather than one number-field specialization, see
[the relative Keller construction](GLOBAL_SUNADA_KELLER_COVERS.md).  It
proves equality of every good fiber zeta function over the common
two-dimensional target.  It deliberately does not claim an absolute
three-dimensional Cox-ledger suspension.

## References

- P. Cassou-Nogues and J.-M. Couveignes,
  [*Factorisations explicites de g(y)-h(z)*](https://matwbn.icm.edu.pl/ksiazki/aa/aa87/aa8741.pdf),
  Acta Arith. 87 (1999), 291--317, especially section 5.1.
- M. D. Fried,
  [*Variables separated equations: strikingly different roles for the Branch
  Cycle Lemma and the finite simple group classification*](https://arxiv.org/abs/1012.5297),
  Sci. China Math. 55 (2012), 1--72.
- W.-C. W. Li and Z. Rudnick,
  [*Pair arithmetical equivalence for quadratic fields*](https://arxiv.org/abs/2007.13147),
  for the classical fact that nonisomorphic arithmetically equivalent fields
  first occur in degree seven.
