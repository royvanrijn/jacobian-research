# Positive-dimensional moduli of exact symplectic Keller lifts

This note combines two independent results already proved in the repository:

1. the explicit degree-five family `F_lambda:A^3->A^3`, whose stable
   polynomial left--right class retains the affine-equivalence class of the
   Hessian-root divisor of its weighted seed; and
2. the exact cotangent lift of a Keller map, which is polynomially
   right-equivalent to the corresponding identity stabilization.

The combination gives a positive-dimensional stable-moduli theorem for exact
symplectic etale endomorphisms of affine six-space.  It also produces a
parameterized family of injective non-surjective Weyl-algebra endomorphisms,
while making no claim that the latter are pairwise inequivalent under arbitrary
Weyl automorphisms.

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
left--right equivalence implies

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

Thus every stable class meets the parameter line in at most six points.

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
equivalent after adjoining any number of identity variables, then (1.1)
holds.  Consequently every stable class meets this family in at most six
parameters, and

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
The degree-five stable-moduli theorem then gives (1.1).  Each orbit has at most
six parameters, while `Lambda` is uncountable, proving (3.1).

## 4. Weyl-algebra shadows

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

No pairwise inequivalence statement is made for the Weyl endomorphisms under
arbitrary pre- and post-composition by Weyl automorphisms.  Such automorphisms
need not preserve the differential-order filtration or the polynomial function
subalgebra, so the stable symplectic invariant does not automatically descend
to that larger equivalence relation.

## 5. What is new in the combination

The inverse-Jacobian Weyl transfer is classical, and the cotangent lift is an
elementary exact symplectic construction.  The new consequence recorded here
is the transfer of the repository's positive-dimensional stable-moduli theorem:

\[
\text{uncountably many stable Keller classes in dimension three}
\quad\Longrightarrow\quad
\text{uncountably many stable exact symplectic classes in dimension six}.
\]

This is stronger than transferring only the earlier divisor-count lower bound
`tau(N-1)`: it shows that fixed dimension, fixed generic degree, exact
symplecticity and Jacobian one still leave a positive-dimensional space of
stable inequivalent nonproper etale maps.

## 6. Reproduction and dependencies

No new symbolic identity is needed beyond the two source theorems.  Their
independent exact checks are:

```bash
.venv/bin/python scripts/verify_degree_five_stable_moduli.py
.venv/bin/python scripts/verify_symplectic_weyl_lift.py
```

The first checks the degree-five seed, Keller map, boundary exclusions,
discriminant saturation and orbit invariant.  The second checks the inverse
Jacobian, Weyl commutators, canonical one-form identity, six-variable
Jacobian and symbolic-momentum collision for the foundational map.

The load-bearing non-computational input is the stable Hessian-root invariant
proved in `DEGREE_FIVE_STABLE_MODULI.md`.  No external specialist review of
that invariant or of the combined stable-moduli consequence is currently
recorded.
