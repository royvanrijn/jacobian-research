# Geometry of the `HC4` exceptional Schur locus

## Status

This note records a first component-directed calculation for the
two-parameter sextic family

\[
h_{\mu,\nu}=\frac{x^6+y^6+z^6}{30}
 +\mu x^2y^2z^2
 +\nu\sum_{i\ne j}x_i^4x_j^2.                 \tag{0.1}
\]

The current result is a **modular reconstruction and exact special-fiber
analysis**, not yet a classification theorem for the full 15-coefficient
quartic.  In the six-dimensional even-quartic subspace, three good primes
reconstruct only the Fermat and radial parameter points.  A quartic lying
in one nontrivial sign-character block is excluded exactly.  Quartics
mixing several sign-character blocks and the exact generic denominator
remain open.

The machine-readable transcript is
[`hc4_exceptional_schur_locus_modular.json`](artifacts/generated-results/hc4_exceptional_schur_locus_modular.json).
Replay its rational reconstruction and exact special-fiber identities with

```bash
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py
```

## 1. Generic system and denominator attempt

On \(\nu\ne0\), the six quotient pivots of `HC4QSG2` eliminate a generic
quadratic quotient.  Clearing \(2\nu^2\) leaves 114 homogeneous quadratic
equations in the fifteen coefficients of a quartic \(s_4\).  Over
\(\mathbb Q(\mu,\nu)\), their 117-element Gröbner basis contains every
coefficient cube.

The transformation-aware command

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --extract-denominators --basis-profile
```

did not finish its `liftstd` calculation in 900 seconds.  This is a
diagnostic only.  In particular, no polynomial printed by an interrupted
calculation is used below.  A new `--extract-basis-denominators` mode
separates the smaller basis-only calculation from the lift-certificate
denominators.

The basis-only calculation finishes and gives

\[
D_{\mathrm{basis}}=2.                              \tag{1.1}
\]

Thus the reduced 117-element basis itself has no parameter-dependent
coefficient denominator.  This does not make the exceptional locus empty:
the denominators can occur in the transformation expressing that basis in
the original 114 generators.  The relevant discriminant is therefore the
lift-certificate denominator, not the coefficient denominator (1.1).

Two narrower exact attempts locate the next computational boundary.  The
primitive degree-two coefficient matrix has rank \(99\) at
\((1,1),(1,2),(2,1),(2,3),(0,1)\).  Its selected \(99\)-minor did not
finish in a one-hour fraction-free determinant run.  At degree three, the
trivial sign-character block has 191 cubic monomials and 441
equation-times-linear-form columns.  A direct lift of \(s_0^3\) over
\(\mathbb Q(\mu,\nu)\) reached the 900-second Singular bound without
returning a certificate.  Neither interrupted calculation supplies a
discriminant factor.

## 2. Signed-permutation decomposition

The signed coordinate group splits the quartic representation into the
six-dimensional even block

\[
\langle x^4,y^4,z^4,x^2y^2,x^2z^2,y^2z^2\rangle
\]

and three three-dimensional nontrivial character blocks, represented by

\[
xy\langle x^2,y^2,z^2\rangle.                 \tag{2.1}
\]

For a quartic contained in (2.1), its Schur quotient is even.  Direct
coefficient comparison gives 36 equations.  Each of the projective charts
\(a=1,b=1,c=1\) for

\[
s_4=xy(ax^2+by^2+cz^2)
\]

has unit ideal over \(\mathbb Q[\mu,\nu]\).  The other two character
blocks follow by permutation.

This does **not** yet show that an arbitrary solution is even: polar cross
terms can couple different character blocks.  That mixed-character
possibility must be checked on the exact generic denominator components.

## 3. Modular even-chart elimination

Write

\[
s_4=ax^4+by^4+cz^4+dx^2y^2+ex^2z^2+fy^2z^2.
\]

The quotient is necessarily
\(q_xx^2+q_yy^2+q_zz^2\), and coefficient comparison again gives 36
equations.  Up to permutation, the projective space is covered by the
two chart types \(a=1\) and \(d=1\).

Elimination modulo \(47,101,103\) has stable Gröbner-basis sizes \(88\)
and \(56\), respectively.  Chinese remaindering and rational
reconstruction on \(a=1\) gives the parameter ideal

\[
\begin{aligned}
&\mu\nu^2-2\nu^3-\frac{\mu^2}{10}
 +\frac{3\mu\nu}{10}-\frac{\nu^2}{5},\\
&\mu^2\nu-4\nu^3-\frac{\mu^2}{2}
 +\frac{8\mu\nu}{5}-\frac{6\nu^2}{5},\\
&\mu^3-8\nu^3-\frac{9\mu^2}{5}
 +6\mu\nu-\frac{24\nu^2}{5},\\
&\nu^4-\frac{\nu^3}{5}-\frac{\mu^2}{100}
 +\frac{\mu\nu}{25}-\frac{3\nu^2}{100}.
                                                        \tag{3.1}
\end{aligned}
\]

Its quotient has length seven and its exact radical is

\[
(\mu-2\nu,\;10\nu^2-\nu),                       \tag{3.2}
\]

with candidate support

\[
(\mu,\nu)=(0,0),\qquad (1/5,1/10).              \tag{3.3}
\]

On \(d=1\), reconstruction gives

\[
\begin{aligned}
&\mu\nu-3\nu^2-\frac{\mu}{10}
 +\frac{2\nu}{5}-\frac1{100},\\
&\mu^2-9\nu^2-\frac{2\mu}{5}
 +\frac{9\nu}{5}-\frac1{20},\\
&\left(\nu-\frac1{10}\right)^3.                  \tag{3.4}
\end{aligned}
\]

Thus this chart is supported only at the radial point.  The powers in
(3.1) and (3.4) explain why a specialization calculation sees a thick
exceptional fiber even though the reconstructed reduced support is
zero-dimensional.

Equations (3.1)--(3.4) are reconstructed from three good characteristics.
They are candidate rational components until exact characteristic-zero
ideal-containment certificates are recorded.

A direct characteristic-zero replay is available as

```bash
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py \
  --exact-pure-chart --singular-timeout 900
```

The \(a=1\) calculation reached that 900-second bound while computing the
36-equation rational standard basis, before returning an elimination
ideal.  It therefore supplies no containment certificate.

## 4. Complex component tracking

For the reconstructed reduced ideal (3.2), a generic complex affine slice
has two constant branches,

\[
(\mu,\nu)=(0,0),\qquad (0.2,0.1).
\]

There is consequently no numerical evidence for an exceptional curve in
the even block.  Reduction at all three good primes preserves the two
rational support points; the second carries the nonreduced transverse
factor \((\nu-1/10)^3\) on the mixed chart.  Since both branches are
rational and stationary, rational reconstruction is unambiguous.

This tracking concerns the reconstructed even-chart parameter ideal, not
the still-uncomputed generic denominator divisor.

## 5. Hessian discriminants and automorphisms

At the Fermat point,

\[
h_{0,0}=\frac{x^6+y^6+z^6}{30},\qquad
\det\operatorname{Hess}(h_{0,0})=x^4y^4z^4.       \tag{5.1}
\]

The Hessian discriminant is the three coordinate lines, each with
multiplicity four.  The projective linear stabilizer is

\[
(\mu_6^2)\rtimes S_3,
\]

of order \(216\).  The polynomial Schur quartics are precisely
\(ax^4+by^4+cz^4\), with quotient
\(16(a^2x^2+b^2y^2+c^2z^2)\).

At the radial point, put \(R=x^2+y^2+z^2\).  Then

\[
h_{1/5,1/10}=\frac{R^3}{30},\qquad
\det\operatorname{Hess}(h_{1/5,1/10})=\frac{R^6}{25}. \tag{5.2}
\]

The discriminant is the smooth conic \(R=0\) with multiplicity six, and
the projective stabilizer is the infinite group
\(\operatorname{PO}_3(\mathbb C)\simeq\operatorname{PGL}_2(\mathbb C)\).

The radial Schur quartic is unique up to scale.  Indeed the inverse
Hessian metric gives, for a homogeneous quartic \(s\),

\[
\nabla s^{\mathsf T}\operatorname{Hess}(h)^{-1}\nabla s
=\frac{5R|\nabla s|^2-64s^2}{R^3}.                \tag{5.3}
\]

Polynomiality first forces \(R\mid s\), say \(s=Rq\).  Equation (5.3)
then becomes

\[
\frac{-4q^2+5R|\nabla q|^2}{R},
\]

so \(R\mid q\), and \(s\) is a scalar multiple of \(R^2\).

## 6. Componentwise lower-face prolongation

Only the two reconstructed reduced strata need lower-face tests:

* The Fermat stratum, including its full projective plane of diagonal
  Schur quartics, is excluded by `HC4QF1`.
* The radial stratum and its unique quartic \(R^2\) are excluded by
  `HC4QSE2`; the incompatible coefficients are
  \(2(3\delta-32)/25\) and
  \(64(99\delta-1040)/25\).

No full-family saturation is used in these prolongation steps.

## 7. Remaining classification gates

The following calculations are still required before promoting (3.3) to
the exceptional-locus classification:

1. compute the transformation/lift denominator, preferably from
   sign-character blocks of the fifteen cube certificates (the basis
   denominator itself is the constant \(2\));
2. compare their reductions in several good characteristics and track
   every complex divisor component;
3. run the full mixed-character projective quartic system only over those
   components;
4. certify (3.1) and (3.4) by characteristic-zero ideal containment.

Until then, the proved statements remain `HC4QF1`, `HC4QSG2`,
`HC4QSE1`, and `HC4QSE2`; this note is a reproducible research frontier.

## 8. Later cube-torsion experiment

The follow-up
[`HC4_FITTING_DENOMINATOR_EXTRACTION.md`](HC4_FITTING_DENOMINATOR_EXTRACTION.md)
constructs the canonical degree-three module behind the fifteen cube
certificates.  The raw 114 quadrics map to the 120-dimensional quadratic
coefficient space; after multiplication by the fifteen coefficient
variables, the relevant presentation is \(A^{1710}\to A^{680}\), with
the cubes defining a separate \(A^{15}\to A^{680}\) map.

On the valid chart \(\nu\ne0\), finite-field scans at four primes
reconstruct the radial point and an additional rational point
\((-5/3,-1/6)\).  Exact specialization shows that the latter is only a
jump from cube to fourth-power nilpotence: all coefficient fourth powers
vanish and the reduced fiber is still the origin.  Thus cube-certificate
torsion does not equal the reduced exceptional Schur locus.  The exact
integral zeroth Fitting ideal and associated-prime equality remain open
after 900-second Singular timeouts.

## 9. Full reduced-fiber fourth-power scan

The cube calculation is refined by testing the fourth powers of all
fifteen quartic coefficients.  The relevant zero-character degree-four
block has 819 monomials and 3474 equation-times-quadratic-monomial
columns.  If every coefficient fourth power belongs to the specialized
homogeneous Schur ideal, then its radical is the coefficient maximal
ideal; equivalently, the reduced projective Schur fiber is empty.

Exhaustive scans on \(D(\nu)\) give:

| field | parameter points | certified empty | remaining point |
|---|---:|---:|---|
| \(\mathbb F_7\) | 42 | 41 | \((3,5)\) |
| \(\mathbb F_{11}\) | 110 | 109 | \((9,10)\) |
| \(\mathbb F_{13}\) | 156 | 155 | \((8,4)\) |

In every row the remaining point is the reduction of
\((1/5,1/10)\).  At that point the degree-four multiplication map has rank
756 and all fifteen tested fourth powers lie outside its image.  The
reductions of the nilpotence-jump point \((-5/3,-1/6)\) are certified
reduced-empty at all three primes.

This is the first full 15-coefficient finite-field test of the reduced
fiber, rather than an even-quartic slice or a cube-torsion proxy.  It does
not yet exclude parameter points defined only over proper extensions of
the finite fields, nor does it compute a characteristic-zero support
ideal.  The transcript is
[`hc4_fourth_power_support.json`](artifacts/generated-results/hc4_fourth_power_support.json).

The direct symbolic annihilator attempt over
\(\mathbb F_7[\mu,\nu]\) reached its 900-second Singular bound before
returning a standard basis.  Thus it supplies no radical or component
factorization and does not close the finite-extension limitation.
