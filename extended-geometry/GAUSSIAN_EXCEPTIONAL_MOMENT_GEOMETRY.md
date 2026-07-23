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

## 2. Scheme-theoretic transport

For every locally closed subscheme `X subset S_N`, define its moment model by

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

## 6. Status and reproduction

The transport theorem adds no independent geometric hypothesis: it is the
functoriality of schemes, normalization, and completed local rings under the
polynomial isomorphism (1).  Its computational inputs are checked by

```bash
.venv/bin/python scripts/verify_gaussian_moment_fingerprint.py
.venv/bin/python scripts/verify_nonsurjective_enumerative_geometry.py
.venv/bin/python scripts/verify_exceptional_partition_lattice.py
.venv/bin/python scripts/verify_degree18_triple_intersection.py
```

The first checker verifies the triangular moment isomorphism; the remaining
three verify the geometric packages transported through it.
