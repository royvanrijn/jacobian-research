# Geometry of the `HC4` exceptional Schur locus

## Status

This note classifies the reduced projective Schur incidence for the
two-parameter sextic family

\[
h_{\mu,\nu}=\frac{x^6+y^6+z^6}{30}
 +\mu x^2y^2z^2
 +\nu\sum_{i\ne j}x_i^4x_j^2.                 \tag{0.1}
\]

For a completely general fifteen-coefficient quartic, the reduced
incidence scheme is

\[
\begin{aligned}
\mathcal I_{\mathrm{red}}
={}&\{(\mu,\nu)=(0,0)\}\times
 \mathbb P\langle x^4,y^4,z^4\rangle\\
&{}\sqcup
\left\{\left(
\frac15,\frac1{10},
[ (x^2+y^2+z^2)^2]
\right)\right\}.                                      \tag{0.2}
\end{aligned}
\]

Thus there is no exceptional parameter curve and no mixed-character
component.  Both components die before the antipodal collision equations
can contribute: the Fermat plane is excluded by `HC4QF1`, and the radial
point by `HC4QSE2`.  This proves `HC4QSE4`.  Formula (0.2) is a statement
about the reduction of the 120-equation projective coefficient scheme;
the nilpotent thickness of that scheme at its two support points is not
classified.

Replay the exact discriminant and incidence atlas with

```bash
.venv/bin/python scripts/verify_hc4_exceptional_schur_atlas.py
```

The machine-readable transcript is
[`hc4_exceptional_schur_locus_modular.json`](artifacts/generated-results/hc4_exceptional_schur_locus_modular.json).
It records the earlier modular route.  Replay that historical
reconstruction and the exact special-fiber identities with

```bash
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py
```

## Exact squarefree-discriminant closure

Put \(X=x^2,Y=y^2,Z=z^2\) and

\[
R=X+Y+Z,\qquad P_2=XY+XZ+YZ,\qquad P_3=XYZ.
\]

The Hessian determinant is

\[
D=c_0R^6+c_1R^4P_2+c_2R^3P_3+c_3R^2P_2^2
  +c_4RP_2P_3+c_5P_2^3+c_6P_3^2,                    \tag{0.3}
\]

where

\[
\begin{aligned}
c_0&=4\nu^2,\\
c_1&=4\nu(\mu-20\nu^2),\\
c_2&=-4(3\mu^2-130\mu\nu^2+13\mu\nu+240\nu^3
             -46\nu^2+\nu),\\
c_3&=-2\nu(10\nu-1)(4\mu-18\nu+1),\\
c_4&=-2(20\mu^2\nu-18\mu^2+780\mu\nu^2-44\mu\nu
             +3\mu-2040\nu^3+240\nu^2-10\nu),\\
c_5&=2(\mu-2\nu)(10\nu-1)^2,\\
c_6&=(2\mu+6\nu-1)^2(10\mu-30\nu+1).
\end{aligned}                                         \tag{0.4}
\]

The key lemma is elementary.  If \(D\) is squarefree and
\[
\nabla s^{\mathsf T}\operatorname{adj}(H)\nabla s\in(D),
\qquad H=\operatorname{Hess}(h_6),
\]
then \(s=0\).  Modulo each irreducible factor \(f\mid D\), squarefreeness
makes \(H\) have generic rank two.  Its symmetric adjugate has rank one,
so the scalar Schur identity forces
\(\operatorname{adj}(H)\nabla s=0\bmod f\).  Hence \(D\) divides this
degree-eleven vector.  The vector is zero, and invertibility of the
adjugate over the ambient fraction field gives \(\nabla s=0\).

Regard (0.3) as \(aP_3^2+bP_3+c\) over
\(\mathbb Q[R,P_2]\).  The radical of the four coefficients of
\(b^2-4ac\) is

\[
(\mu,\nu)\cap(5\mu-1,10\nu-1)
\cap(7\mu-5,14\nu+1).                                \tag{0.5}
\]

The third point is harmless:

\[
D_{5/7,-1/14}
=\frac{(12P_2-R^2)(144P_2^2-24P_2R^2-7R^4)}{343},    \tag{0.6}
\]

and the two factors are distinct and coprime.  On the two lines
\(c_6=0\), the coefficient \(b\) is squarefree away from the points in
(0.5); this excludes a repeated factor when the quadratic drops to a
linear polynomial in \(P_3\).

Finally, the maps
\((x,y,z)\mapsto(X,Y,Z)\mapsto(R,P_2,P_3)\) are ramified over
\(P_3=0\) and the cubic-root discriminant.  Exact coefficient comparison
shows that \(P_3\mid D\) only at Fermat and that \(D\) is never a scalar
multiple of the cubic-root discriminant.  Thus pullback introduces no
additional repeated component.  The nonsquarefree locus is exactly
Fermat plus radial, and the exact fiber classifications give (0.2).

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
They remain modular rather than characteristic-zero ideal-containment
certificates, but their reduced support is now independently proved by
the exact squarefree-discriminant argument above.

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

This tracking concerns the reconstructed even-chart parameter ideal.
The exact full-coefficient classification no longer depends on the
uncomputed generic denominator divisor.

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

Only the two reduced strata in (0.2) need lower-face tests:

| component | first continued face | second continued face | result |
|---|---|---|---|
| Fermat, \(s_4=ax^4+by^4+cz^4\) | \(\lambda^{13}t x^4y^4z^4\) forces \(\delta=32(a^3+b^3+c^3)/3\) | three \(\lambda^{11}t^3\) coefficients become \(1024a^5,1024b^5,1024c^5\) | projectively empty |
| radial, \(s_4=R^2\) | \(\lambda^{13}t x^{12}\) is \(2(3\delta-32)/25\) | \(\lambda^{11}t^3x^8\) is then \(1024/25\) | empty |

No full-family saturation is used in these prolongation steps.
Both calculations retain arbitrary lower homogeneous coefficients inside
the collision-normalized potential.  Their determinant-face ideals are
already empty on the projective Schur components, so adjoining the
antipodal collision equations cannot create a solution.

## 7. Closed classification gates

The squarefree-discriminant argument above bypasses the four expensive
gates proposed in the original component-directed calculation:

1. lift-certificate denominators are unnecessary for reduced support;
2. proper finite-field extension points are covered in characteristic
   zero by (0.5);
3. mixed sign-character quartics are excluded whenever \(D\) is
   squarefree, without projective chart elimination;
4. the exact Fermat and radial fiber calculations replace the modular
   containments (3.1) and (3.4).

This promotes the reduced atlas to theorem `HC4QSE4`.  The remaining
scheme-theoretic question is only the nilpotent/embedded thickness of the
120-equation incidence scheme at the two support points; it cannot create
a new reduced Schur pair or an `HC(4)` candidate.

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

This was the first full 15-coefficient finite-field test of the reduced
fiber, rather than an even-quartic slice or a cube-torsion proxy.  By
itself it did not exclude parameter points defined only over proper
extensions of the finite fields.  The exact argument (0.3)--(0.6) now
closes that limitation.  The transcript is
[`hc4_fourth_power_support.json`](artifacts/generated-results/hc4_fourth_power_support.json).

The historical direct symbolic annihilator attempt over
\(\mathbb F_7[\mu,\nu]\) reached its 900-second Singular bound before
returning a standard basis.  It supplies no certificate and is not used
by `HC4QSE4`.
