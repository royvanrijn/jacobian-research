# Automorphisms and imprimitive slices of the incidence cover

Work over an algebraically closed field `k` of characteristic zero.  Let
`H in k[W]` have degree `N>=3`, put

\[
 E_H(W;s,t)=H(W)-sW+t,
\]

and let `pi_H:X_H -> A^2_(s,t)` be the finite degree-`N` incidence cover
cut out by `E_H=0`.  This note separates three notions which should not be
called by the same name:

1. deck transformations over the fixed two-dimensional target;
2. self-equivalences of the pencil which also move `(s,t)`; and
3. imprimitivity after restricting the pencil to a vertical line `s=s_0`.

The first locus is empty in every degree `N>=3`.  The second has a complete
elementary classification by the affine stabilizer of the Hessian divisor.
The third is the generally larger Ritt decomposition locus of `H-s_0W`.

## 1. There are no exceptional two-parameter deck groups

The [universal pencil theorem](../verified/WEIGHTED_SEED_THEOREM.md) gives
geometric and arithmetic monodromy `S_N`
for **every** degree-`N` polynomial `H`, not merely for a generic seed.  The
stabilizer of one sheet in the natural action is `S_(N-1)`, which is
self-normalizing for `N>=3`.  Therefore

\[
 \operatorname{Aut}_{\mathbb A^2_{s,t}}(X_H)=1.       \tag{1.1}
\]

The natural `S_N` action is also primitive.  Consequently there is no seed
for which the generic two-parameter cover is imprimitive.  Parity,
decomposability, and Chebyshev type do not create exceptions to (1.1).
They concern automorphisms which move the base or monodromy after a
one-dimensional specialization.

## 2. Pencil-compatible affine self-equivalences

Let `L=k[W]_(<=1)` and write `[H]` for the class of `H` in `k[W]/L`.  Define

\[
 \Gamma_{\rm pen}(H)=
 \{\phi\in\operatorname{Aff}_1:
       [H\circ\phi]=\gamma[H]\text{ for some }\gamma\in k^*\}. \tag{2.1}
\]

Thus `phi(W)=alpha W+beta` belongs to `Gamma_pen(H)` exactly when

\[
 H(\alpha W+\beta)=\gamma H(W)+uW+v                 \tag{2.2}
\]

for scalars `alpha,gamma!=0` and `beta,u,v`.  Equation (2.2) induces the
cover self-equivalence

\[
 W'=\alpha W+\beta,\qquad
 s'=\frac{\gamma s+u}{\alpha},\qquad
 t'=\gamma t-v+\beta s',                             \tag{2.3}
\]

because

\[
 E_H(W';s',t')=\gamma E_H(W;s,t).                    \tag{2.4}
\]

These are not deck transformations unless the transformation of `(s,t)` is
the identity.

### Hessian-divisor criterion

Twice differentiating (2.2) gives

\[
 \alpha^2H''(\alpha W+\beta)=\gamma H''(W).          \tag{2.5}
\]

Conversely, (2.5) integrates twice to (2.2).  Hence

\[
 \boxed{\Gamma_{\rm pen}(H)
 =\operatorname{Stab}_{\operatorname{Aff}_1}
       \bigl(\operatorname{div}(H'')\bigr),}          \tag{2.6}
\]

where the effective Hessian divisor retains all multiplicities.  This is
exactly the Fitting divisor already present in the decorated normalization.
It gives a direct scheme-theoretic detector for the symmetry locus.

The same group describes affine reparametrizations of the discriminant
normalization compatible with its plane embedding.  Indeed the repeated-root
curve is parametrized by

\[
 \nu_H(r)=\bigl(H'(r),rH'(r)-H(r)\bigr),             \tag{2.7}
\]

and (2.2) gives the corresponding affine target transformation.  Conversely,
a projective symmetry preserving the affine chart and the pole filtration at
the unique point over infinity induces an affine change of `r`; preservation
of the pulled-back Fitting divisor then gives (2.5).  Thus “extra projective
automorphism” should be qualified by these natural boundary conditions.

### Complete classification

Suppose first that `Gamma_pen(H)` is nontrivial.  A nonidentity translation
cannot stabilize the finite divisor of `H''`, so every nonidentity element
has a fixed point.  A finite subgroup of `Aff_1` in characteristic zero has a
common fixed point and is cyclic.  After writing `W=c+X` and subtracting the
tangent line of `H` at `c`, put

\[
 F_c(X)=H(c+X)-H(c)-H'(c)X=\sum_{j\ge2}a_jX^j.       \tag{2.8}
\]

There are precisely two cases.

* If `F_c(X)=aX^N`, then
  `Gamma_pen(H)=G_m`, acting by `X -> lambda X`.
* Otherwise put
  \[
   q=\gcd\{j-j':j\ne j',\ a_j a_{j'}\ne0\}.         \tag{2.9}
  \]
  The group is `mu_q`, acting by `X -> zeta X`.  It is nontrivial exactly
  when `q>1`.

Equivalently, the finite-symmetry normal form is

\[
 F_c(X)=X^eR(X^q),\qquad e\ge2,\quad q>1,            \tag{2.10}
\]

with at least two nonzero monomials.  The exact stabilizer is `mu_q` when
the differences of the exponents occurring in `R` have gcd one.  For a
primitive `q`-th root `zeta`, equation (2.2) has
`alpha=zeta` and `gamma=zeta^e`.

This proves that all pencil-compatible affine symmetry groups are cyclic,
apart from the centered-monomial `G_m` case.  It also gives coefficient
equations for the symmetry strata: translate by an unknown center `c`, remove
the constant and linear terms, and require all remaining exponents to lie in
one residue class modulo `q`.

This is the one-variable self-equivalence classification appearing, in a
closely related left-right form, as Lemma 3.17 of Zieve--Müller,
[*On Ritt's polynomial decomposition theorems*](https://arxiv.org/abs/0807.3578).
Here the quotient by linear polynomials in (2.1) is forced by the `-sW+t`
pencil.

## 3. Intersection with the normalized admissible seed space

Recall

\[
 \mathcal A_N=
 \{H:H(0)=H'(0)=H(1)=0,\ H'(1)=-1,\ H''(1)\ne-2\}. \tag{3.1}
\]

The general centered normal form can be imposed directly.  Start with

\[
 F(W)=(W-c)^eR((W-c)^q)
\]

or `F(W)=(W-c)^N` in the infinite case, and set

\[
 H(W)=A\{F(W)-F(0)-F'(0)W\}.                         \tag{3.2}
\]

The remaining endpoint conditions are

\[
 F(1)-F(0)-F'(0)=0,\qquad
 A=-\frac1{F'(1)-F'(0)},                             \tag{3.3}
\]

followed by the open condition `AF''(1)!=-2`.  Equations (3.2)--(3.3)
parametrize every symmetric admissible seed, including symmetries whose
center is not the boundary mark zero.

For the orbifold locus of the coarse decorated-normalization quotient, the
marks `0` and `infinity` force `c=0`.  If zero has exact multiplicity `e`,
the strata simplify to

\[
 \boxed{H(W)=W^eR(W^q),\quad
 N=e+mq,\quad R(0)\ne0,\quad R(1)=0,\quad qR'(1)=-1.} \tag{3.4}
\]

For fixed `(e,q,m)`, the locus with stabilizer containing `mu_q` has
dimension `m-1`; deleting the larger congruence subloci gives exact
stabilizer `mu_q`.  On the exact-double-zero locus, `e=2`, so possible
orders satisfy

\[
 q\mid N-2.                                          \tag{3.5}
\]

In particular the canonical seed

\[
 H_N(W)=\frac1{N-2}W^2(1-W^{N-2})                   \tag{3.6}
\]

has stabilizer `mu_(N-2)`.  Its nonzero primitive roots are precisely the
rerootings which collapse to the same normalized seed.  Thus (3.4) is an
explicit orbifold locus for the **coarse** decorated-normalization map.

The infinite-symmetry points are obtained from `F(W)=(W-c)^N`.  Their centers
satisfy the finite equation

\[
 (1-c)^N-(-c)^N-N(-c)^{N-1}=0,                      \tag{3.7}
\]

after which (3.3) fixes the scale.  For example, when `N=3`, `c=1/3`
gives the normalized cubic `H=W^2-W^3`.

The full affine-root-sheet decoration behaves differently from the coarse
quotient.  The
[decorated-normalization theorem](DECORATED_NORMALIZATION_INVARIANT.md)
shows that on the exact-double, boundary-clean locus, a nonidentity rerooting
sends the distinguished root-one affine sheet to an extra-root boundary
sheet.  It therefore does not preserve the regular-reconstruction open.
Consequently the cyclic stabilizers in (3.4) are orbifold stabilizers of the
coarse decoration, not counterexamples to the proved marked affine-sheet
faithfulness.  A failure of marked faithfulness could occur only where the
boundary-clean sheet distinction itself degenerates; that smaller
intersection must be tested using the reconstruction valuations, not merely
the functional equation.

## 4. Vertical specializations and Ritt imprimitivity

Fix `s_0` and consider the polynomial cover

\[
 f_{s_0}:\mathbb A^1_W\longrightarrow\mathbb A^1,
 \qquad f_{s_0}(W)=H(W)-s_0W.                        \tag{4.1}
\]

The monodromy of `f_(s_0)(W)-z` is imprimitive if and only if

\[
 f_{s_0}=A\circ B,qquad \deg A>1,\quad\deg B>1.     \tag{4.2}
\]

Indeed intermediate blocks in the monodromy action correspond to
intermediate fields between `k(f_(s_0)(W))` and `k(W)`; Lüroth's theorem and
total ramification at infinity turn such an intermediate field into a
polynomial decomposition.  This is the standard monodromy formulation of
Ritt decomposition; see Section 2, especially Lemma 2.2, of the
[Zieve--Müller paper](https://arxiv.org/abs/0807.3578).

Thus the imprimitive-specialization locus is

\[
 \mathcal D_{a,b}=
 \{H\in\mathcal A_N:\exists s_0,A,B,
   \ H-s_0W=A\circ B,\ \deg A=a,\deg B=b\},
 \qquad ab=N.                                        \tag{4.3}
\]

It is computable by normalizing the inner polynomial `B`, comparing
coefficients, and eliminating the coefficients of `A,B` and `s_0`.
Ritt moves describe intersections between the different `(a,b)` strata.
In particular:

* if `N` is prime, no vertical specialization is imprimitive;
* if `H` is decomposable, then `s_0=0` lies in (4.3);
* a Chebyshev specialization has dihedral monodromy, but is imprimitive only
  when its degree is composite; and
* the quartic case is ubiquitous: after translating a quartic to remove its
  cubic term, the free subtraction of `s_0W` removes its linear term, leaving
  a quadratic polynomial in `W^2`.

The functional equation (2.2) does **not** classify all of (4.3).  It
classifies the smaller locus where the decomposition has a cyclic affine
symmetry.  More precisely, the vertical cover has a nontrivial deck group
exactly when, after an affine change of `W`,

\[
 f_{s_0}(W)=A((W-c)^q),\qquad q>1.                   \tag{4.4}
\]

Equivalently there is a nonidentity solution of (2.2) with

\[
 \gamma=1,\qquad u=(\alpha-1)s_0,\qquad v=\beta s_0. \tag{4.5}
\]

Every case (4.4) is imprimitive, but a general composition $A\circ B$ has no
nontrivial deck transformation.  Therefore symmetry and imprimitivity should
be kept as two separate exceptional loci.

## 5. The examples in the question

The distinctions can be read off immediately.

* `H(-W)=H(W)` gives `gamma=1`.  At `s_0=0`, the vertical cover factors
  through `W^2` and has the deck involution `W -> -W`.
* `H(-W)=-H(W)` gives `gamma=-1`.  It acts on the pencil by
  `(W,s,t)->(-W,s,-t)`; it is not a target-fixed deck transformation and need
  not make an odd prime-degree vertical cover imprimitive.
* A decomposable seed gives an imprimitive `s_0=0` slice, usually without any
  affine symmetry.
* Chebyshev type controls the smaller dihedral monodromy of a vertical slice.
  Its visible affine self-equivalence is only the parity involution; dihedral
  monodromy itself is not a deck group of the non-Galois degree-`N` cover.
* Extra pencil-compatible projective symmetries of the discriminant are
  detected by the same full Hessian-divisor stabilizer (2.6).

## 6. Resulting classification programme

The automorphism half is complete: compute the affine stabilizer of the
effective divisor `(H'')`, or equivalently use the centered congruence normal
forms (2.10).  For normalized coarse moduli, equations (3.4)--(3.5) give the
orbifold strata explicitly.

The imprimitive half should be organized independently by the factor pairs
`ab=N` and the elimination loci (4.3).  A practical next step is to compute
their defining ideals and Ritt-intersection diagram degree by degree, then
intersect those loci with the boundary-degeneration conditions governing the
distinguished affine root sheet.  That intersection, rather than the whole
Hessian symmetry locus, is the candidate exceptional locus for marked
decorated-normalization faithfulness.
