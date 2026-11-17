# Bounded quantization from scalar Moyal anomalies

14 November 2026. **QCL1.** Every polynomial symplectic endomorphism
of total degree at most three, in every Weyl rank, has an exact
quantization with the same spatial degree bound and rational scalar
coefficients in hbar regular at zero. Thus every cubic classical point
belongs to the horizontal part of the degree-three relation scheme.

More generally, in any degree a spatially bounded tuple whose normalized
Moyal commutators are scalars can be normalized by a rational scalar
matrix. This is a positive bounded-realization construction, rather
than another order-by-order existence argument.

It supplies a concrete sufficient target for Long: a degree-three
representative under finite quantum-compatible coordinate changes, or
a finite scalar-anomaly tuple with the required classical reduction.
No such representative is constructed here.

## A scalar commutator matrix can be corrected without spatial growth

Use the standard Moyal product with
\[
 \Pi=\sum_i(\partial_{p_i}\otimes\partial_{x_i}
                    -\partial_{x_i}\otimes\partial_{p_i}),\qquad
 a\star b=m\exp(h\Pi/2)(a\otimes b).
\]
Order the coordinates as (x_1,p_1,...,x_n,p_n), and let J be the
block diagonal matrix with blocks
\[
                         \begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Let F be a column of 2n finite spatial polynomials, with scalar
coefficients rational in h and regular at h=0. Suppose
\[
       h^{-1}[F_i,F_j]_\star=B_{ij}(h),\qquad
                  B(0)=J,                             \tag{1}
\]
where every B_ij is scalar, with no spatial variables. Then there is
a rational matrix A(h), regular at zero with A(0)=1, such that
\[
                         A B A^t=J.                    \tag{2}
\]
Consequently
\[
                             \widehat F=A(h)F
\]
satisfies all Weyl relations exactly and retains the spatial degree
bound of F.

For completeness, construct A by symplectic Gram--Schmidt over K(h).
Start with the standard ordered basis. If a previously chosen pair
(u,v) has pairing u B v^t=-1, replace each remaining row w by
\[
                   w+(w B v^t)u-(w B u^t)v .
\]
It is orthogonal to that pair. For the next two residual rows (a,b),
keep a and replace b by -b/(a B b^t). The denominator is -1 at h=0.
Induction gives (2), and all divisions are by rational functions
nonzero at zero. Starting with the standard paired basis makes
A(0)=1. This proves the claim over K(h), with no field extension or
infinite spatial expansion.

One can also use a closed even formal expression when
B=J+h^2 C. Put L=1-h^2 C J, so B=LJ. Skew symmetry gives
LJ=JL^t. Then the binomial matrix series
\[
                         A=L^{-1/2}                    \tag{3}
\]
satisfies A J=J A^t and A^2 L=1, which imply ABA^t=J.
The rational Gram--Schmidt construction avoids requiring a square
root or termination of the series (3).

## Cubic symplectic maps automatically have scalar anomalies

Let f=(f_1,...,f_(2n)) be polynomial with total degree at most three
and
\[
                            \{f_i,f_j\}=J_{ij}.
\]
The odd Moyal expansion terminates at its third contraction:
\[
 \frac{[f_i,f_j]_\star}{h}
       =\{f_i,f_j\}+\frac{h^2}{24}
                    m\Pi^3(f_i\otimes f_j)
       =J_{ij}+h^2 C_{ij}.                              \tag{4}
\]
Every third derivative of a cubic is constant, so C is a scalar
skew matrix. Equations (1)--(2) give an exact full quantization
\[
                    \widehat f_i=\sum_j A_{ij}(h)f_j,
       \qquad \deg_{\rm spatial}\widehat f_i\le3,
       \qquad \widehat f_i(0)=f_i.                      \tag{5}
\]
This verifies all n(2n-1) independent commutator relations
simultaneously. Passing from symmetric symbols to normal PBW order
does not increase spatial degree.

The argument concerns a specified finite canonical tuple. It does
not define a Lie-algebra quantization preserving the Poisson bracket
of every cubic polynomial or every polynomial observable.

For the finite-type relation schemes of
[QFL2](BOUNDED_REALIZATION.md), (5) is a rational unramified curve
through every special cubic point, with h nonconstant. Thus every
such point is in H_3, the closure of h nonzero. This does not assert
that the entire relation scheme is flat or reduced.

By QFL2, a genuinely noninvertible cubic polynomial symplectic map
would yield a nonsurjective Weyl endomorphism of the same rank.
This statement neither constructs such a classical map nor proves
that all general cubic symplectic maps are invertible.

## A single homogeneous symplectic layer is invertible in every degree

Suppose
\[
                           f(z)=z+H(z)
\]
is symplectic and H is homogeneous of any degree d>=2. Write N=DH.
Comparing the distinct homogeneous degrees d-1 and 2d-2 in
\[
                    (1+N)J(1+N)^t=J
\]
gives
\[
                       NJ+JN^t=0,\qquad NJN^t=0.
\]
The first equation makes the second equal to -N^2J=0, so
\[
                               N^2=0.                  \tag{6}
\]
Euler's identity gives Nz=dH, hence NH=N^2z/d=0.

Let delta=sum_i H_i partial_(z_i). Then delta(H)=0 and delta(z)=H.
On any polynomial of degree e its (e+1)st iterate vanishes:
with H held constant by delta, repeated differentiation is the
finite directional Taylor expansion. Thus delta is locally
nilpotent and
\[
                    f=\exp(delta),\qquad
                    f^{-1}(z)=z-H(z).                  \tag{7}
\]
This excludes noninvertibility for every single homogeneous nonlinear
layer, not just for cubics. After removing a constant translation and
a symplectic linear part, it also applies to every map of degree at
most three with no quadratic homogeneous layer.

Accordingly, any noninvertible cubic symplectic map must have an
essential quadratic layer after tangent normalization. A reduction
to the form z+H_3 cannot preserve both symplecticity and
noninvertibility. This statement is valid in every rank.

## Why cubic Keller counterexamples do not automatically supply the input

Our [cubic Keller constructions](../../jacobian/cubic-reduction/README.md) satisfy
the determinant-one condition. That condition alone does not make
a map symplectic or put it under (4). A cotangent lift restores
symplecticity but can increase its degree.

This increase is especially explicit in homogeneous Keller form. Let
\[
                     F(z)=z+H_d(z),\qquad \det DF=1,
                     \qquad N=DH_d,\qquad d>=2.
\]
Homogeneity separates the terms of degrees d-1,2(d-1),... in
det(1+N)=1. Therefore the characteristic polynomial of N is a pure
power, N is nilpotent, and, if r is its nilpotency index,
\[
                       (DF)^{-1}=\sum_{j=0}^{r-1}(-N)^j.
 \tag{8}
\]
If N^2=0, the same Euler and locally nilpotent argument as above
makes F invertible. A noninvertible homogeneous Keller map of this
form must therefore have r>=3.

The nonzero terms N^j have distinct homogeneous degrees j(d-1).
Hence the cotangent lift
\[
              (z,p)\longmapsto
                  \bigl(F(z),(DF(z))^{-t}p\bigr)
\]
has exact total degree
\[
                    (d-1)(r-1)+1\ge2d-1               \tag{9}
\]
for a noninvertible homogeneous F. For d=3 this is 2r-1>=5.
There is no cancellation between those homogeneous degrees or between
distinct momentum variables. Thus homogenizing an existing Keller
collision to cubic Keller form and then taking its cotangent lift
cannot put it in the cubic symplectic class. The condition r>=3 is necessary for noninvertibility, not
sufficient; triangular automorphisms can already have r=3.

This is a structural separation between the two degree reductions,
not a new result about the minimum dimension of Keller collisions.

## The next bounded-lift gate

The scalar-anomaly criterion (1) applies in every spatial degree.
For Long, its displayed uncorrected Moyal tuple has nonconstant
defects, as recorded in the
[ordering audit](../dixmier-audit/NEXT_EXPERIMENTS.md), so that tuple
does not meet (1). It remains possible to seek a different bounded
tuple or finite quantum-compatible coordinate representative.

A degree-three representative would settle bounded realizability
immediately by (5), without calculating higher quantum corrections.
A homogeneous cubic representative without a quadratic layer is
excluded by (7). Finite source and target quantum automorphisms
with bounded-degree inverses preserve bounded realizability, as
proved in QFL2. Arbitrary formal gauges do not preserve a chosen
spatial bound.

Alternatively the
[right-symbol criterion](HAMILTONIAN_IDEALIZER_CRITERION.md) and
[compatible filtered-frame criterion](PROJECTIVE_FRAME_BOUNDARY.md)
remain ways to construct a bounded model from module data. No
right-symbol transversality or compatible frame for Long is supplied
by (4) alone.

[Exact controls](check_scalar_anomaly.py) check the cubic Moyal
formula, all six relations of a mixed rank-two cubic canonical
tuple, rational normalization of a nontrivial scalar skew matrix,
the homogeneous inverse identities and an index-three cotangent
degree-five example. The universal all-rank statements and the all-degree
homogeneous arguments are written proofs. No independent
theorem replay, formal verification, external review, novelty
priority, noninvertible cubic symplectic example, bounded Long lift
or DC2 solution is claimed.
