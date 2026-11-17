# Positive-dimensional moduli of exact symplectic Keller lifts

This note combines three results already proved in the repository:

1. the decorated-normalization theorem, whose generically finite invariant
   has image dimension `N-3` on a nonempty ordinary boundary-clean weighted
   locus for every `N>=4`;
2. the explicit degree-five family `F_lambda:A^3->A^3`, for which full-cover
   faithfulness recovers `lambda` exactly while the Hessian-root divisor gives
   an independent coarse invariant; and
3. the exact cotangent lift of a Keller map, which is polynomially
   right-equivalent to the corresponding identity stabilization.

The combination gives an `(N-3)`-dimensional stable-moduli theorem for exact
symplectic etale endomorphisms of affine six-space in every generic degree
`N>=4`.  The cotangent construction also produces `(N-3)`-dimensional moduli
of injective non-surjective Weyl-algebra endomorphisms under
differential-order-preserving equivalence.  No corresponding claim is made
under arbitrary Weyl automorphisms.

## 1. The degree-five family

Let `Lambda` and

\[
F_\lambda:\mathbb A^3_\mathbb C\longrightarrow\mathbb A^3_\mathbb C,
\qquad \lambda\in\Lambda,
\]

be the explicit family in
[Degree-five stable moduli](DEGREE_FIVE_STABLE_MODULI.md).  Its weighted seed is

\[
H_\lambda(W)=
\frac{W^2(W-1)\bigl(3W^2-(5\lambda+1)W+3\lambda\bigr)}{60},
\]

with

\[
H_\lambda''(W)=
\frac{(W-\lambda)(10W^2-8W+1)}{10}.
\]

Put

\[
p=\frac{4-\sqrt6}{10},\qquad
q=\frac{4+\sqrt6}{10},\qquad
\chi(\lambda)=\frac{\lambda-p}{q-p}.
\]

The degree-five stable-moduli theorem proves that every `F_lambda` is a Keller
map of generic degree five with a five-point fiber, and that stable polynomial
left--right equivalence implies `lambda=mu`.  Independently, its coarse
Hessian/Fitting invariant implies

\[
\chi(\mu)\in
\left\{
\chi(\lambda),\ 1-\chi(\lambda),\ \frac1{\chi(\lambda)},
\ \frac1{1-\chi(\lambda)},\
\frac{\chi(\lambda)}{\chi(\lambda)-1},\
\frac{\chi(\lambda)-1}{\chi(\lambda)}
\right\}.
\tag{1.1}
\]

Thus (1.1) remains an explicit six-element coarse orbit, although full-cover
faithfulness separates every parameter exactly.


## 2. Exact cotangent lifts

Write

\[
J_\lambda=DF_\lambda,
\qquad
B_\lambda=J_\lambda^{-1}.
\]

Because `det J_lambda` is a nonzero constant, `B_lambda` has polynomial
entries.  Define

\[
\widehat F_\lambda:T^*\mathbb A^3\longrightarrow T^*\mathbb A^3,
\qquad
\widehat F_\lambda(x,p)
=\bigl(F_\lambda(x),B_\lambda(x)^Tp\bigr).
\tag{2.1}
\]

The general cotangent-lift theorem in
[Exact symplectic and Weyl lifts](SYMPLECTIC_WEYL_LIFT.md) gives

\[
\widehat F_\lambda^*(p^Tdx)=p^Tdx,
\qquad
\det D\widehat F_\lambda=1,
\tag{2.2}
\]

and the polynomial factorization

\[
\widehat F_\lambda
=(F_\lambda\times\operatorname{id}_{\mathbb A^3})\circ H_\lambda,
\qquad
H_\lambda(x,p)=(x,B_\lambda(x)^Tp),
\tag{2.3}
\]

where `H_lambda` is a polynomial automorphism with inverse

\[
H_\lambda^{-1}(x,q)=(x,J_\lambda(x)^Tq).
\tag{2.4}
\]

In particular, the cotangent lift is right-equivalent to the three-coordinate
identity stabilization of the original Keller map.

## 3. Main theorem

### Theorem 3.1 — uncountable exact symplectic stable moduli

For every `lambda in Lambda`, the map

\[
\widehat F_\lambda:\mathbb A^6_\mathbb C
\longrightarrow\mathbb A^6_\mathbb C
\]

is a noninjective, nonproper, exact symplectic etale polynomial endomorphism
with Jacobian one and generic degree five.  Every target momentum over a
five-point target fiber of `F_lambda` has five distinct preimages under
`widehat F_lambda`.

If `widehat F_lambda` and `widehat F_mu` are polynomially left--right
equivalent after adjoining any number of identity variables, then
`lambda=mu`.  Consequently

\[
\boxed{
\text{there are uncountably many pairwise stably polynomially
left--right inequivalent exact symplectic etale maps }
\mathbb A^6\to\mathbb A^6
\text{ of generic degree five.}
}
\tag{3.1}
\]

#### Proof

Exact symplecticity, etaleness, Jacobian one, generic degree, fiber schemes and
properness follow from the cotangent-lift theorem and (2.3).  More explicitly,
for every target `(y,q)` there is a canonical scheme isomorphism

\[
F_\lambda^{-1}(y)
\longrightarrow
\widehat F_\lambda^{-1}(y,q),
\qquad
x\longmapsto\bigl(x,J_\lambda(x)^Tq\bigr).
\tag{3.2}
\]

Hence a five-point fiber of `F_lambda` gives a five-point fiber of the lift over
every target momentum `q`, and noninjectivity and nonproperness transfer.

Suppose the two lifted maps become polynomially left--right equivalent after
adjoining `s` identity variables.  By (2.3), this is an equivalence between

\[
F_\lambda\times\operatorname{id}_{\mathbb A^{3+s}}
\quad\text{and}\quad
F_\mu\times\operatorname{id}_{\mathbb A^{3+s}}.
\]

Thus `F_lambda` and `F_mu` are stably polynomially left--right equivalent.
Generic affine-mark faithfulness then gives `lambda=mu`.  Since `Lambda` is
uncountable, this proves (3.1).  Equation (1.1) remains a weaker independent
consequence of the Hessian/Fitting divisor.

## 3A. Degreewise dimension theorem

### Theorem 3.2 — `(N-3)`-dimensional exact symplectic stable moduli

For every `N>=4`, there is a nonempty Zariski-open family of weighted
degree-`N` Keller maps whose cotangent lifts

\[
 \widehat F_H:\mathbb A^6\longrightarrow\mathbb A^6
\]

form an `(N-3)`-dimensional family of stable polynomial left--right classes
of noninjective, nonproper, exact symplectic etale maps with Jacobian one and
generic degree `N`.

#### Proof

The normalized admissible degree-`N` seed space has dimension `N-3`.  The
[decorated-normalization theorem](DECORATED_NORMALIZATION_INVARIANT.md) proves
that its ordinary boundary-clean locus is nonempty for every `N>=4`, and that
the decorated invariant map on this locus is generically finite with image
dimension `N-3`.  Because that invariant is preserved by stable polynomial
left--right equivalence, its image gives an `(N-3)`-dimensional family of
stable Keller classes.

For every seed `H`, the cotangent factorization (2.3) identifies
`widehat F_H` up to polynomial right equivalence with
`F_H x id_{A^3}`.  If two lifts became left--right equivalent after any
further stabilization, their base maps would be stably left--right
equivalent, and hence would have the same decorated invariant.  The lifted
family therefore has the same `N-3` lower bound on stable-class dimension.
Exact symplecticity, Jacobian one, generic degree, fiber schemes, and
nonproperness transfer by the general cotangent-lift theorem.

## 4. Filtered Weyl moduli

For each `lambda`, the inverse-Jacobian derivations define an endomorphism of
the third Weyl algebra:

\[
\Phi_\lambda(X_i)=F_{\lambda,i}(X),
\qquad
\Phi_\lambda(D_i)=
\sum_r(B_\lambda)_{ri}(X)D_r.
\tag{4.1}
\]

Every `Phi_lambda` is injective and non-surjective, and for the
differential-order filtration

\[
\operatorname{gr}(\Phi_\lambda)=\widehat F_\lambda^*.
\tag{4.2}
\]

Therefore the family supplies uncountably many explicit parameter values
producing injective non-surjective endomorphisms of `A_3(C)`, whose associated
graded maps range through uncountably many stable polynomial left--right
classes of exact symplectic maps.

Call two Weyl endomorphisms **stably filtered left--right equivalent** if,
after adjoining finitely many Weyl pairs on which both maps act identically,
they differ by pre- and post-composition with automorphisms preserving the
differential-order filtration in both directions.

### Theorem 4.1 — `(N-3)`-dimensional filtered Weyl moduli

For every `N>=4`, the inverse-Jacobian construction on the ordinary
boundary-clean degree-`N` seed open gives an `(N-3)`-dimensional family of
stable filtered left--right equivalence classes of injective non-surjective
endomorphisms of `A_3(C)`.  On the affine-marked cover of that open, filtered
equivalence compatible with the marked root sheet is parameter-faithful: it
forces equality of the normalized seeds.

To prove this, let `alpha` be a differential-order-preserving automorphism of
`A_n`.  Its restriction to the order-zero algebra `C[X_1,...,X_n]` is a
polynomial automorphism `a`.  Write the order-one symbol of `alpha(D_i)` as

\[
 \sum_r c_{ri}(x)\xi_r.
\]

The relations
`[alpha(D_i),alpha(X_j)]=delta_(ij)` say exactly that the matrix `C` is the
inverse-Jacobian matrix of `a`.  Hence

\[
 \operatorname{gr}(\alpha)=\widehat a^*,                       \tag{4.3}
\]

the polynomial cotangent lift of the base automorphism.  The same argument
applies after adjoining Weyl pairs.  Consequently a stable filtered
left--right equivalence between `Phi_H` and `Phi_K` induces a stable
polynomial left--right equivalence between their associated-graded cotangent
lifts `widehat F_H` and `widehat F_K`.  By the factorization (2.3), the base
maps `F_H` and `F_K` are then stably equivalent.  Decorated normalization
therefore gives the `N-3` lower bound, and generic affine-mark faithfulness
recovers `H` exactly on the marked open.

This conclusion is intentionally filtration-relative.  No pairwise
inequivalence statement is made under arbitrary pre- and post-composition by
Weyl automorphisms: such automorphisms need not preserve the
differential-order filtration or the polynomial function subalgebra, so the
associated-graded invariant need not descend to that larger equivalence
relation.

## 5. What is new in the combination

The inverse-Jacobian Weyl transfer is classical, and the cotangent lift is an
elementary exact symplectic construction.  The main new consequences recorded
here are the degreewise transfers

\[
\dim\overline{\operatorname{im}(\mathfrak h_N)}=N-3
\quad\Longrightarrow\quad
\begin{matrix}
\text{an }(N-3)\text{-dimensional family of stable exact symplectic classes
in dimension six},\\
\text{an }(N-3)\text{-dimensional family of stable filtered Weyl classes
in rank three}.
\end{matrix}
\]

This holds for every `N>=4`.  It is a moduli-dimension statement, distinct
from the finite arithmetic count
`1+(N-1)tau(N-1)-sigma(N-1)`, which additionally separates cancellation
branches from the weighted construction.

## 6. Reproduction and dependencies

No new symbolic identity is needed beyond the two source theorems.  Their
independent exact checks are:

```bash
.venv/bin/python scripts/verify_degree_five_stable_moduli.py
.venv/bin/python scripts/verify_symplectic_weyl_lift.py
```

The first checks the explicit degree-five seed, Keller map, boundary
exclusions, discriminant saturation and orbit invariant.  The last checks
the inverse Jacobian, Weyl commutators, canonical one-form identity,
six-variable Jacobian and symbolic-momentum collision for the foundational
map.

The load-bearing non-computational input for exact parameter recovery is
generic affine-mark faithfulness in `DECORATED_NORMALIZATION_INVARIANT.md`; the stable
Hessian-root invariant in `DEGREE_FIVE_STABLE_MODULI.md` supplies the
independent coarse orbit.  No external specialist review of either invariant
or of the combined stable-moduli consequence is currently recorded.  In
particular, Christopher D. Long's independent consequence papers do not
validate or review this theorem.
