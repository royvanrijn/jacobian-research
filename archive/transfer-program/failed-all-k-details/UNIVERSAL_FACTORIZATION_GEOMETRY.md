# Universal factorization slices and affine complements

This note isolates the part of Factorization-Slice Obstruction that does not use the all-`k` transfer
claim Failed All-k Transfer Block.  Its subject is the following affine-geometry problem:

> For which factor degrees `(a,b)` and hyperplane contacts is the marked
> factorization complement `X_(a,b)(ell)` isomorphic to affine space?

There are four logically separate layers.

1. Sections 1--3 give unconditional universal factorization geometry.
2. Section 4 gives the exact cubic affine-space classification.
3. Section 5 gives an exact higher-dimensional negative example.
4. Sections 6--7 record related transfer and weighted constructions, but no
   assertion in those sections is used in Sections 1--5.

Thus a failure of Failed All-k Transfer Block does not affect the cubic theorem, the `(2,3)` point
count, or its affine-space obstruction.  The ground field is algebraically
closed of characteristic zero for the geometric classification, unless a
finite-field model is explicitly being counted.

Here "affine-space" means *isomorphic to* `A^(a+b)`, not merely affine as a
scheme.  These are different questions: the complements below have explicit
affine presentations in the consecutive-degree cases, while the main issue
is whether those affine varieties are actually affine spaces.

## 1. Universal factorization maps

For positive integers `a,b`, put `n=a+b` and consider

\[
 \mu_{a,b}:\operatorname{Sym}^a(\mathbb P^1)\times
 \operatorname{Sym}^b(\mathbb P^1)
 \longrightarrow\operatorname{Sym}^n(\mathbb P^1),\qquad
 (D,E)\longmapsto D+E.                                      \tag{1}
\]

After identifying `Sym^d(P^1)` with `P^d`, this is multiplication of binary
forms.  It is finite of generic degree `binom(n,a)`.  Its ramification divisor
is the resultant divisor

\[
 \mathcal R_{a,b}=\{([Q],[R]):\operatorname{Res}(Q,R)=0\};     \tag{2}
\]

away from it, a reduced divisor of degree `n` is split into two disjoint
subdivisors of degrees `a` and `b`.

Let `ell` be a nonzero linear functional on binary forms of degree `n`, and
write `H_ell` for its projective hyperplane.  The affine factorization
complement is

\[
 X_{a,b}(\ell)=
 (\mathbb P^a\times\mathbb P^b)\setminus
 (\mathcal R_{a,b}\cup\mu_{a,b}^{-1}(H_\ell)).                 \tag{3}
\]

The target complement is `P^n minus H_ell ~= A^n`, and (1) restricts to an
étale generically `binom(n,a)`-to-one map from (3) to `A^n`.  The affine-space
question is whether (3) is isomorphic to `A^n`.

## 2. The torus-normalized slice

Write

\[
 Q=\sum_{i=0}^a q_iU^{a-i}V^i,\qquad
 R=\sum_{j=0}^b r_jU^{b-j}V^j,
\]

and write

\[
 \ell\left(\sum_{k=0}^n m_kU^{n-k}V^k\right)
 =\sum_{k=0}^n\ell_km_k.
\]

The complement of the pulled-back hyperplane has the affine normalization

\[
 \widetilde X_{a,b}(\ell)=
 \left\{(Q,R):
   \sum_{i=0}^a\sum_{j=0}^b\ell_{i+j}q_ir_j=1,
   \ \operatorname{Res}(Q,R)\ne0\right\}.                     \tag{4}
\]

The torus acts freely by

\[
 \lambda\cdot(Q,R)=(\lambda Q,\lambda^{-1}R),                 \tag{5}
\]

and

\[
 X_{a,b}(\ell)=\widetilde X_{a,b}(\ell)/\mathbb G_m.          \tag{6}
\]

The resultant has weight `b-a` under (5).  If `d=b-a>0`, the section
`Res(Q,R)=1` meets every geometric orbit and leaves the residual action of
`mu_d`.  Consequently

\[
 X_{a,b}(\ell)\cong
 \left\{
 \sum_{i,j}\ell_{i+j}q_ir_j=1,
 \ \operatorname{Res}(Q,R)=1
 \right\}/\mu_d.                                               \tag{7}
\]

In particular, for consecutive degrees `b=a+1` there is no residual finite
quotient: (7) is an explicit affine complete intersection in
`A^(a+b+2)`.

There is also a quick divisor-class calculation.  In
`Pic(P^a x P^b)=Z^2`, the two boundary divisors in (3) have classes

\[
 [\mu^{-1}(H_\ell)]=(1,1),\qquad [\mathcal R_{a,b}]=(b,a).      \tag{8}
\]

When both are irreducible, the boundary sequence shows

\[
 \mathcal O(X_{a,b})^*/k^*=0,\qquad
 \operatorname{Pic}(X_{a,b})\cong\mathbb Z/|b-a|\mathbb Z.     \tag{9}
\]

For `b=a+1`, the units are constant and the Picard and divisor class groups
vanish.  Thus these elementary invariants do not distinguish the consecutive
factorization candidates from affine space.

## 3. Hyperplane contact types

Restrict `ell` to the small diagonal, the rational normal curve

\[
 \nu_n:[\alpha:\beta]\longmapsto
 [\alpha^n:\binom n1\alpha^{n-1}\beta:\cdots:\beta^n].        \tag{10}
\]

The resulting binary form

\[
 h_\ell(\alpha,\beta)
 =\sum_{k=0}^n\ell_k\binom nk\alpha^{n-k}\beta^k              \tag{11}
\]

classifies the contact divisor `H_ell cap nu_n(P^1)`.  Conversely, if
`h=sum h_k alpha^(n-k) beta^k`, then the corresponding normalized target and
factorization slices are

\[
 \sum_{k=0}^n{h_k\over\binom nk}m_k=1,
 \qquad
 \sum_{i=0}^a\sum_{j=0}^b
 {h_{i+j}\over\binom n{i+j}}q_ir_j=1.                         \tag{12}
\]

Thus contact types are partitions of `n`, together with the moduli of their
distinct support points modulo `PGL_2`.  A partition alone is a complete
orbit invariant only when the support has at most three points.  Formula
(12), rather than the partition by itself, is the general normalized slice.

## 4. Exact cubic affine-space classification

For `n=3`, `PGL_2` has exactly three hyperplane orbits, represented by

| Contact type | `h_ell(alpha,beta)` | normalized target slice |
|---|---|---|
| `(1,1,1)` | `alpha beta(alpha-beta)` | `m_1-m_2=1` |
| `(2,1)` | `alpha^2 beta` | `m_1=1` |
| `(3)` | `alpha^3` | `m_0=1` |

Take `(a,b)=(1,2)`.  Over a finite field with the displayed support points
rational, fix the marked linear factor `Q`.  Coprimality removes one line
from the residual `P^2`, and the hyperplane condition removes a second line.

- In type `(1,1,1)` the two lines are distinct for every marked point, so
  every one of the `q+1` fibers has `q^2-q` points.
- In type `(2,1)` the two lines coincide at the tangent point and are distinct
  at the other `q` points.
- In type `(3)` the hyperplane functional vanishes identically on the fiber
  over the osculation point, leaving an empty fiber; the other `q` fibers
  contain `q^2-q` points.

Therefore

\[
\begin{array}{c|c}
\text{contact type}&\#X_{1,2}(\ell)(\mathbb F_q)\\ \hline
(1,1,1)&q^3-q\\
(2,1)&q^3\\
(3)&q^3-q^2.
\end{array}                                                     \tag{13}
\]

The same fiber decomposition gives the corresponding Grothendieck classes
`L^3-L`, `L^3`, and `L^3-L^2`.  The middle slice is the tangent nonosculating
slice of Cubic Marked-Root Model, and the two-chart reconstruction proves that it is `A^3`.
The other two classes, or just their point counts at every good reduction,
exclude an isomorphism with `A^3`.  Hence:

> **Cubic uniqueness theorem.**  Over an algebraically closed field of
> characteristic zero, and up to automorphisms of `P^1` and rescaling the
> hyperplane equation, the tangent nonosculating contact type `(2,1)` is the
> unique hyperplane for which the unramified cubic factorization complement
> is isomorphic to affine three-space.

The positive assertion uses only the explicit two-chart theorem of Cubic Marked-Root Model.  The
two negative assertions use the displayed motivic classes (equivalently,
their good-reduction point counts).  Neither direction uses Failed All-k Transfer Block or any
transfer theorem.

## 5. Exact negative example: the `(2,3)` middle slice

Now take `(a,b)=(2,3)` and let `ell` extract the middle coefficient
`U^2V^3`; its contact divisor on the small diagonal has type `(2,3)`.  The
torus slice (7) is

\[
 \boxed{
 q_0r_3+q_1r_2+q_2r_1=1,
 \qquad \operatorname{Res}(Q,R)=1 }
 \ \subset\mathbb A^7.                                        \tag{14}
\]

It is a smooth affine fivefold.  Equations (8)--(9) show that it has only
constant units, trivial Picard and class groups, and trivial canonical class.
The coordinate-scaling torus has exactly one fixed coprime factorization in
the middle slice, so its Euler characteristic is one.  All of these agree
with `A^5`.

The decisive invariant comes from good reductions.  For a bilinear form of
rank `r` on vector spaces of dimensions `u,v`, the number of projective pairs
on which it is nonzero is

\[
 \left(\#\mathbb P^{u-1}-\#\mathbb P^{u-r-1}\right)q^{v-1}.    \tag{15}
\]

Apply homogeneous-polynomial Moebius inversion to the common divisor `G` of
`Q` and `R`.  Only `deg G=0,1,2` occur.

- For `deg G=0`, the middle catalecticant has rank three, contributing
  `(q^2+q+1)q^3`.
- For `deg G=1`, one of the `q+1` linear divisors gives rank one and the
  other `q` give rank two.  Their Moebius weight is `-1`, so the total
  contribution is `-q^3-q(q+1)q^2`.
- For `deg G=2`, every squarefree divisor gives a nonzero rank-one form.
  The sum of their Moebius weights is the coefficient `q` in
  `1/Z(P^1,t)=(1-t)(1-qt)`, giving contribution `q^2`.

Adding the three terms yields

\[
 \boxed{\#X_{2,3}(\ell)(\mathbb F_q)=q^5-q^3+q^2.}             \tag{16}
\]

This differs from `#A^5(F_q)=q^5` at every good reduction.  Any
characteristic-zero isomorphism is defined by finitely many coefficients and
relations; spreading them out and specializing gives an isomorphism over a
number field and then at almost every finite place, contradicting (16).
Thus the `(2,3)` candidate is not affine five-space.  The count is
polynomial; the standard polynomial-count comparison also gives compactly
supported Hodge--Deligne polynomial `L^5-L^3+L^2`.

## 6. Suspended transfer connection (not part of the independent theorem)

Let `Poly_d^mon` denote monic degree-`d` polynomials.  The strong transfer
scheme is the fiber product

\[
 \operatorname{Poly}_{3k}^{\rm mon}
 \mathop{\times}_{\operatorname{Poly}_{6k}^{\rm mon}}
 \operatorname{Poly}_{2k}^{\rm mon},                            \tag{17}
\]

where the two maps are `U mapsto U^2` and `V mapsto V^3`.  Its reduced
normalization locus is

\[
 S\longmapsto(U,V)=(S^3,S^2),\qquad\deg S=k.                    \tag{18}
\]

The formal completion of (17) along (18) is exactly `Z_k` by definition.  The
affine variant is the analogous fiber product after projecting degree-`6k`
polynomials modulo `span(1,Z)`.  The differential Wronskian argument proving
`Z_k^aff=Z_k` in characteristic zero is separate from the refuted ribbon
model.  The uniform finite-flat rank assertion formerly attached to this
statement is open; see
[ALL_K_CONDUCTOR_RIBBON_AUDIT.md](ALL_K_CONDUCTOR_RIBBON_AUDIT.md).  Thus the
transfer block is not merely analogous to a factorization singularity; it is
the local fiber of the square/cube universal factorization equalizer.  No
rank, flatness, or Boolean-structure conclusion from Failed All-k Transfer Block is imported into
the affine-complement results of Sections 1--5.

## 7. Separate weighted pullback connection

The universal marked-root factorization is `mu_(1,n-1)`.  For an admissible
weighted seed `H`, define the three-dimensional parameter map

\[
 \sigma_H:\mathbb A^3_{A,B,C}\longrightarrow
 \operatorname{Poly}_n,
 \qquad
 (A,B,C)\longmapsto E_{A,B,C}(W)
 =H(W)-BCW+cAC^2.                                      \tag{19}
\]

Pulling `mu_(1,n-1)` back along `sigma_H` gives the root incidence

\[
 \mathcal I_H=\{((A,B,C),W):E_{A,B,C}(W)=0\}.                   \tag{20}
\]

After normalization, retain the open on which the marked linear factor
reconstructs `(x,y,z)`.  The weighted marked-root theorem is precisely the
isomorphism

\[
 \mathbb A^3_{x,y,z}\xrightarrow{\sim}\mathcal R_H
 \subset\mathcal I_H^\nu,                                      \tag{21}
\]

and the weighted Keller map is the projection back to the parameter
threefold.  Although the polynomial in (19) depends on `(A,B,C)` through the
two pencil coordinates `BC` and `AC^2`, the third coordinate records the
reconstruction scale.  Thus the weighted family is a three-dimensional
parametrized pullback slice of the universal factorization geometry, not an
embedded three-plane in the full coefficient space.

This interpretation depends on the weighted marked-root theorem weighted marked-root theorem.  It is
context for the broader programme, not an input to the independent Factorization-Slice Obstruction
classification and obstruction results.

## 8. Resulting programme

The independent problem now has a clean sequence of targets.

1. Compute `[X_(a,b)(ell)]` in the Grothendieck ring, beginning with contact
   divisors supported at at most three points and then including support
   moduli.
2. Derive finite-field point-count polynomials (or detect non-polynomial
   counts) by homogeneous Moebius inversion.
3. Compute logarithmic invariants of the boundary
   `R_(a,b) union mu^(-1)(H_ell)`, including log Kodaira dimension and log
   Hodge data.
4. Combine these with units, class groups, canonical classes, topology, and
   cancellation invariants to characterize when a slice can be affine
   space.
5. Extend from two marked factors to ordered multi-factor multiplication
   maps and compare the resulting higher-dimensional complements.

The cubic `(2,1)` slice is the positive base case.  The `(2,3)` middle slice
is the first consecutive-degree warning: constant units, trivial class and
canonical groups, and Euler characteristic one do not characterize affine
space, while the finite-field/Grothendieck correction does.

## Verification

Run

```text
.venv/bin/python scripts/verify_universal_factorization_geometry.py
```

The script checks the contact normal forms, resultant weights, all three
cubic counts, and the `(2,3)` Moebius/rank formula by exact finite-field
enumeration.  Its final two checks concern the optional connections in
Sections 6--7 and are not dependencies of the affine-space results.
