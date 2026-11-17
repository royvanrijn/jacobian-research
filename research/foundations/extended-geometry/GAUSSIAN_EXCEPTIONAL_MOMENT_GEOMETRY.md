# Exceptional seed geometry in Gaussian mixed-moment coordinates

This note combines the optimal-coordinate theorem in the
[weighted Gaussian bridge](WEIGHTED_GAUSSIAN_BRIDGE.md) with the affine
geometry of the nonsurjective seed locus.  The conclusion is formal but
strong: the exceptional partition-lattice geometry occurs unchanged inside
an explicit family of four-real-Gaussian `GMC` witnesses, and its coordinates
are finitely observable mixed moments.

## 1. The moment-coordinate isomorphism

Work over a characteristic-zero field, fix `N>=4`, and fix a nonzero bridge
scale `lambda`.  Let `S_N` be the ambient normalized affine seed slice

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1,qquad\deg H\le N.
\]

It is affine `(N-3)`-space.  For the four-real-Gaussian pair
`(P_(H,lambda),Q)` constructed by the bridge, put

\[
 M_m(H,\lambda)=\mathbb E(QP_{H,\lambda}^m).
\]

The optimal-coordinate theorem gives a polynomial isomorphism

\[
 \boxed{
 \Psi_{N,\lambda}:\mathcal S_N\xrightarrow{\sim}\mathbb A^{N-3},
 \qquad
 H\longmapsto(M_3,\ldots,M_{N-1}).}                 \tag{1}
\]

Indeed, after `mu_m=M_m/m!`, one has the triangular formulas

\[
 \mu_{j+1}=\lambda c_j+F_j(c_2,\ldots,c_{j-1}),
 \qquad 2\le j\le N-2,                              \tag{2}
\]

and the two endpoint equations recover `c_(N-1),c_N`.  The diagonal of (2)
is nonzero, so both (1) and its inverse are polynomial.

Retaining exact degree, weighted admissibility, or any of the ordinary and
boundary-clean conditions merely restricts (1) to the corresponding open
subsets.

## 2. The Gaussian invariant-transport theorem

> **Theorem 2.1 (Gaussian invariant transport).**  The isomorphism
> `Psi_(N,lambda)` transports every functorially defined seed-space locus to
> the optimal mixed-moment space, in the category in which that locus is
> defined.  In particular:
>
> 1. closed and locally closed subschemes, their intersections, and their
>    local and formal scheme structures transport scheme-theoretically;
> 2. constructible loci transport constructibly;
> 3. over `R`, semialgebraic loci and their relative Euclidean boundaries
>    transport semialgebraically; and
> 4. after spreading out the triangular isomorphism, arithmetic conditions
>    on the recovered seed or its associated covers transport over every
>    good residue field.

For a locally closed subscheme `X subset S_N`, define its moment model by

\[
 X^{\mathrm{mom}}=\Psi_{N,\lambda}(X).               \tag{3}
\]

Because (1) is an isomorphism, (3) preserves more than sets and dimensions:

\[
 \mathcal O_{X,x}\simeq
 \mathcal O_{X^{\mathrm{mom}},\Psi(x)},\qquad
 \widehat{\mathcal O}_{X,x}\simeq
 \widehat{\mathcal O}_{X^{\mathrm{mom}},\Psi(x)}.    \tag{4}
\]

It also preserves reduced structures, nilpotents, normalization morphisms,
conductors, scheme intersections, tangent cones, and transverse Artin
algebras.  Explicitly, if

\[
 \widetilde X\longrightarrow X
\]

is a normalization used in the seed geometry, then

\[
 \widetilde X\longrightarrow X
 \xrightarrow{\Psi_{N,\lambda}}X^{\mathrm{mom}}      \tag{5}
\]

is the normalization of the moment model.

This assertion is affine.  A nonlinear polynomial coordinate isomorphism
need not preserve the degree of a chosen projective closure, and no
projective seed-closure degree is inferred here.

For closed `X=V(I)`, the assertion is completely effective.  If

\[
 \Psi_{N,\lambda}^{\#}:
 k[M_3,\ldots,M_{N-1}]\xrightarrow{\sim}k[\mathcal S_N]
\]

is the triangular coordinate-ring isomorphism, then

\[
 I(X^{\rm mom})=(\Psi_{N,\lambda}^{\#})^{-1}(I).   \tag{5a}
\]

Thus equations are obtained by finite triangular substitution, without a
new elimination problem.  Constructible transport follows because an
isomorphism preserves finite Boolean combinations of closed subsets.  Over
`R`, `Psi` and its polynomial inverse are semialgebraic homeomorphisms, which
also proves the boundary assertion.

The same statement holds for loci with auxiliary parameters: for every
parameter scheme `T`, use `Psi_(N,lambda) times id_T`.  This is the form
needed for real root-count chambers in target space and for factorization or
Frobenius conditions on specialized covers.

For the arithmetic statement, choose a finitely generated coefficient ring
for `lambda` and invert `lambda` and the finitely many denominators in the
factorial-normalized formulas.  Equations (2) then spread `Psi` and its
inverse to an isomorphism over that ring.  Reduction at every remaining good
prime preserves rational points, residue-field extensions, and every
splitting condition computed from the recovered seed.  This is the precise
sense in which finite-field conditions transfer; it is not a claim at the
excluded primes.

The theorem also applies to stabilizer and orbifold strata by transporting
their groupoid or inertia-locus presentation through `Psi`.  It does not by
itself prove that an informally specified complexity class is algebraic or
constructible.  An `LR`-complexity locus transfers once such a locus has
first been defined in one of the categories above.

> **Corollary 2.2.**  Hessian symmetry and orbifold loci, vertical Ritt
> decomposition loci, real-sheet chambers and their walls, and good-prime
> splitting conditions all have canonical optimal-moment models.  Any
> explicit seed-coordinate equations or tests for these loci therefore give
> explicit equations or tests on the corresponding four-real-Gaussian
> witness family.

## 3. The nonsurjective locus in moment space

Let `N_N subset S_N` be the normalized nonsurjective seed locus and set

\[
 \mathcal N_N^{\mathrm{mom}}
 =\Psi_{N,\lambda}(\mathcal N_N).                   \tag{6}
\]

Its irreducible components are

\[
 \mathcal C_{a,b}^{\mathrm{mom}}
 =\Psi_{N,\lambda}(\mathcal C_{a,b}),
 \qquad 2a+3b=N.                                    \tag{7}
\]

Therefore the following seed-space results hold verbatim in the coordinates
`(M_3,...,M_(N-1))`:

* the period-six component count;
* the consecutive dimension and codimension staircases;
* the explicit smooth component normalizations;
* the collision-stratum and common-coarsening intersection posets;
* the binomial normalization-branch profiles; and
* all generic minimal-intersection transverse algebras.

For example, two distinct components `C_(a,b)^mom` and `C_(c,d)^mom` have
unique minimal common coarsening

\[
 \nu_{\min}=6^t3^{\min(b,d)}2^{\min(a,c)},
 \qquad
 t=\frac{|a-c|}{3}=\frac{|b-d|}{2},                 \tag{8}
\]

and, at the geometric generic point of that stratum, transverse algebra

\[
 \boxed{
 k[\epsilon_1,\ldots,\epsilon_t]/
 (\epsilon_1^2,\ldots,\epsilon_t^2),
 \qquad\operatorname{length}=2^t.}                  \tag{9}
\]

Thus the partition-lattice geometry in moment space is not an analogy: it is
the same affine scheme geometry under the explicit isomorphism (1).

## 4. The degree-eighteen three-sheet face

For `N=18`, the component `C_(3,4)` has three normalization branches over
the exact collision stratum `E_(6,6,6)`.  The proved ordered-triple transverse
algebra transports to

\[
 \boxed{
 D_{123}^{\mathrm{tr}}
 \simeq k[\epsilon_1,\epsilon_2,\epsilon_3]/
 (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2),
 \qquad\operatorname{length}=8.}                    \tag{10}
\]

In particular its cubic socle class

\[
 \epsilon_1\epsilon_2\epsilon_3\ne0
\]

is visible in the moment-space local algebra and disappears from every
pairwise slice.

The scope of (10) is important.  It is the ordered-triple transverse
intersection algebra on the proved nonempty open.  Reconstruction of the
entire completed three-branch component ring and its total conductor
equalizer remains open.  Moment coordinates transport that open problem but
do not solve it.

## 5. Realization inside Gaussian coefficient space

For fixed `N` and `lambda`, all pairs `(P_(H,lambda),Q)` lie in a
finite-dimensional affine coefficient space `V_(N,lambda)`.  The bridge
construction gives a polynomial morphism

\[
 \iota_{N,\lambda}:\mathcal S_N\longrightarrow
 \mathcal V_{N,\lambda}.                             \tag{11}
\]

Wick contraction makes each mixed moment a polynomial function on `V`.  By
applying the triangular inverse of (1) to those functions, one obtains a
polynomial retraction

\[
 R_{N,\lambda}:\mathcal V_{N,\lambda}\longrightarrow\mathcal S_N,
 \qquad R_{N,\lambda}\circ\iota_{N,\lambda}=\mathrm{id}.        \tag{12}
\]

Hence `iota_(N,lambda)` is a closed immersion.  The image of (6) is therefore
an explicit stratified subfamily of four-real-Gaussian `GMC` counterexamples,
not merely an abstract reparameterization.  Its seed parameters and all the
affine exceptional strata above are recovered from the finite observable
vector `(M_3,...,M_(N-1))`.

Here closed immersion follows because (11) is a section of the morphism
(12), and sections of morphisms between affine schemes are closed
immersions.

This does not classify the image modulo arbitrary transformations of the
Gaussian variables or a broader equivalence of polynomial pairs.  It proves
parameter faithfulness and exact affine scheme geometry before taking any
such quotient.

## 6. The first explicit moment-space partition complex: degree six

The first degree with two exceptional components is `N=6`:

\[
 \mathcal C_{3,0}^{\mathrm{mom}}
 \quad\text{and}\quad
 \mathcal C_{0,2}^{\mathrm{mom}}.
\]

Set the bridge scale to `lambda=1` and use factorial-normalized moment
coordinates

\[
 u=\mu_3=\frac{M_3}{3!},\qquad
 v=\mu_4=\frac{M_4}{4!},\qquad
 \omega=\mu_5=\frac{M_5}{5!}.                       \tag{13}
\]

If `H=sum_(i=2)^6 c_iW^i`, then triangular recovery is especially simple:

\[
 u=c_2,\qquad v=c_3,\qquad\omega=c_4+2c_2^2.         \tag{14}
\]

The all-double component is the irreducible sextic surface

\[
 \mathcal C_{3,0}^{\mathrm{mom}}=V(\mathscr F_6),    \tag{15}
\]

where

\[
\begin{aligned}
\mathscr F_6={}&
16\omega^2-32\omega^2uv-48\omega^2u^2\\
&+192\omega u^4+128\omega u^3v-320\omega u^3
 -384\omega u^2v+128\omega u^2\\
&-104\omega uv^2+112\omega uv+24\omega u
 +8\omega v^3+40\omega v-16\omega\\
&-192u^6-128u^5v+640u^5+768u^4v-768u^4\\
&+208u^3v^2-992u^3v+400u^3
 -16u^2v^3-384u^2v^2+400u^2v-64u^2\\
&-16uv^3+80uv^2+16uv-16u
 +21v^4-28v^3+46v^2-28v+5.
\end{aligned}                                       \tag{16}
\]

The degree rises from four in the seed coordinates to six in the moment
coordinates.  This is a concrete reminder that a nonlinear affine
automorphism preserves local rings and intersection algebras, but not the
ordinary degree of a coordinate presentation.

The all-triple component has the following rational normalization chart.  Put

\[
 A=15q^3-15q^2+6q-1,\qquad B=5q^2-5q+1.
\]

Then

\[
\begin{aligned}
u(q)&={A\over3qB},\\
v(q)&=-{(3q^2-3q+1)(45q^4-15q^2+6q-1)
       \over27q^4B},\\
\omega(q)&={A(30q^5-30q^4-3q^3+18q^2-8q+1)
       \over9q^4B^2}.
\end{aligned}                                       \tag{17}
\]

The two components meet at the four all-six collision points.  On the chart
(17), their reduced support is

\[
 \Delta_6(q)=45q^4-30q^3+15q^2-6q+1=0.             \tag{18}
\]

The polynomial `Delta_6` is squarefree, its roots avoid the denominators in
(17), and the four corresponding all-six seeds are weighted-admissible.  The
scheme-theoretic intersection is visible directly from

\[
 \boxed{
 \mathscr F_6(u(q),v(q),\omega(q))
 =-{(3q-1)^8\Delta_6(q)^2
    \over3^{11}q^{16}B^4}.}                         \tag{19}
\]

Thus each of the four intersection points has local length two.  Equation
(19) is the first completely explicit Gaussian-moment realization of the
primitive dual-number weight in the exceptional partition complex.

## 7. Vertical Ritt loci among degree-six Gaussian witnesses

The degree-six solution of
[vertical OP-RITT](HESSIAN_RITT_INCIDENCE.md) now gives the first application
of Theorem 2.1 beyond the nonsurjective partition geometry.  Continue with
`lambda=1` and the moment coordinates `(u,v,omega)` from (13).  In the
normalized seed chart

\[
 H=\sum_{j=3}^6h_j(W^j-W^2),
\]

the triangular inverse is

\[
\begin{aligned}
h_3&=v,\\
h_4&=\omega-2u^2,\\
h_5&=1-4u-3v-2\omega+4u^2,\\
h_6&=\omega-2u^2+3u+2v-1.
\end{aligned}                                      \tag{20}
\]

Let `R_(a,b)^mom` be the locus where some vertical specialization
`H-s_0W` factors as `A o B` with `(deg A,deg B)=(a,b)`.  The two composite
factor orders are hypersurfaces

\[
 \boxed{
 \mathcal R_{2,3}^{\rm mom}=V(\mathscr F_6),\qquad
 \mathcal R_{3,2}^{\rm mom}=V(\mathscr G_6),}       \tag{21}
\]

where `F_6` is exactly the polynomial in (16), and

\[
\begin{aligned}
\mathscr G_6={}&
-4\omega^3+24\omega^2u^2-60\omega^2u-27\omega^2v+6\omega^2\\
&-48\omega u^4+240\omega u^3+108\omega u^2v-288\omega u^2
 -252\omega uv+114\omega u\\
&-54\omega v^2+36\omega v-12\omega
 +32u^6-240u^5-108u^4v+552u^4\\
&+504u^3v-548u^3+108u^2v^2-549u^2v+264u^2
 -216uv^2+198uv-60u\\
&-27v^3+27v^2-18v+5.
\end{aligned}                                      \tag{22}
\]

This produces a sharper coincidence than transport alone predicts:

\[
 \boxed{\mathcal R_{2,3}^{\rm mom}
 =\mathcal C_{3,0}^{\rm mom}}.                     \tag{23}
\]

Thus, for normalized sextic seeds, the entire all-double nonsurjective
component is exactly the locus whose pencil has a quadratic-after-cubic
imprimitive vertical specialization.  The other factor order gives the new
sextic surface `V(G_6)`.

Their exact-degree scheme intersection is

\[
 V(\mathscr F_6,\mathscr G_6).                     \tag{24}
\]

On the marked clean open it is the one-dimensional monomial/Chebyshev Ritt
collision family.  An explicit incidence chart is obtained from

\[
 g_{c,q}(W)=((W-c)^3+q(W-c))^2
\]

subject to

\[
15c^4-20c^3+12c^2q+15c^2-8cq-6c+q^2+2q+1=0,       \tag{25}
\]

by normalizing

\[
 H_{c,q}(W)=
 -\frac{g_{c,q}(W)-g_{c,q}(0)-g'_{c,q}(0)W}
 {g'_{c,q}(1)-g'_{c,q}(0)}                         \tag{26}
\]

and then applying `(u,v,omega)=(c_2,c_3,c_4+2c_2^2)`.  The rational clean
witness `(c,q)=(2/5,-1/5)` maps to

\[
 (u,v,\omega)=\left(-\frac14,4,-\frac{99}{8}\right),             \tag{27}
\]

which lies on both hypersurfaces in (21).  Consequently low-degree OP-RITT
does not merely promise future Gaussian equations: degree six already gives
explicit imprimitive-specialization equations among four-real-Gaussian
`GMC` witnesses.

## 8. Status and reproduction

The transport theorem adds no independent geometric hypothesis: it is the
functoriality of schemes, normalization, and completed local rings under the
polynomial isomorphism (1).  Its computational inputs are checked by

```bash
.venv/bin/python scripts/verify_gaussian_moment_fingerprint.py
.venv/bin/python scripts/verify_nonsurjective_enumerative_geometry.py
.venv/bin/python scripts/verify_exceptional_partition_lattice.py
.venv/bin/python scripts/verify_degree18_triple_intersection.py
.venv/bin/python scripts/verify_degree_six_gaussian_moment_geometry.py
.venv/bin/python scripts/verify_hessian_ritt_degree_six.py
```

The first checker verifies the triangular moment isomorphism; the remaining
checkers verify the geometric packages transported through it, the explicit
degree-six partition equations (13)--(19), and the Ritt equations (20)--(27).
