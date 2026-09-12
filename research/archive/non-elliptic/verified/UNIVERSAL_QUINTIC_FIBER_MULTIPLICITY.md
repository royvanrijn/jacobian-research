# Universal quintic fiber multiplicity

This note proves the rank-five case of universal stable multiplicity for
Keller fibers.  Unlike the quartic proof, it does not use weighted tangent
chords or a rational-point theorem: ordinary translation of one suitably
chosen primitive generator moves nontrivially in the one-dimensional stable
moduli of quintic quadratic-gauge maps.

## The theorem

> **Universal quintic multiplicity theorem.**  
> For every characteristic-zero field `K` and every rank-five finite etale
> `K`-algebra `A`,
> \[
>  \boxed{|\mathcal R_K(A)|=\infty.}
> \]
> The infinitely many classes may all be represented by determinant-one,
> geometric-degree-five quadratic-gauge Keller maps
> `A^3_K -> A^3_K`.

Together with the
[quartic multiplicity theorem](UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md),
this proves the stronger universal claim in ranks four and five, with
infinite multiplicity rather than the requested lower bound three.

## 1. Choosing a nonexceptional primitive generator

Let

\[
 A_0=\ker(\operatorname{Tr}_{A/K}).
\]

The trace pairing on `A` is nondegenerate, and
`\operatorname{Tr}(1)=5`, so its restriction to `A_0=1^\perp` is
nondegenerate.  The primitive elements form a nonempty Zariski open in
`A_0`: translate any primitive element by one fifth of its trace.

Consequently there is a primitive `eta in A_0` satisfying

\[
 \operatorname{Tr}(\eta^2)\ne0.                        \tag{1.1}
\]

Indeed the nonprimitive locus and the trace-isotropic quadric are proper
closed subsets of the four-dimensional affine space `A_0`.  Their
complement has a `K`-point because every characteristic-zero field is
infinite.

Write the characteristic polynomial of `eta` as

\[
 P(T)=T^5+aT^3+bT^2+cT+d.
                                                               \tag{1.2}
\]

Newton's identity gives

\[
 \boxed{a=-\frac12\operatorname{Tr}(\eta^2)\ne0.}       \tag{1.3}
\]

The nonzero coefficient `a` is the only generator choice needed below.

## 2. Translation and the quadratic-gauge stable invariant

For `s in K`, use the translated primitive generator `T=s+S` and put

\[
 G_s(S)=P(s+S)-P(s)
       =g_1(s)S+\cdots+g_5(s)S^5.                      \tag{2.1}
\]

The relevant derivative jets are

\[
 g_1=P'(s),\qquad
 g_3=\frac{P'''(s)}6=10s^2+a,\qquad
 g_4=\frac{P''''(s)}{24}=5s,\qquad
 g_5=1.                                                \tag{2.2}
\]

The quadratic coefficient `g_2` is removable by a polynomial target shear.
After dividing by `g_1`, the coefficient-torus coordinates are

\[
 a_3=\frac{g_3}{g_1},\qquad
 a_4=\frac{g_4}{g_1},\qquad
 a_5=\frac{g_5}{g_1}.
\]

For degree five, the
[quadratic-gauge stable-moduli theorem](QUADRATIC_GAUGE_STABLE_MODULI.md)
has weights

\[
 (-2,-1),\qquad(-3,-4),\qquad(-4,-5).
\]

Their primitive integer relation is `(-1,-6,5)`.  Hence the stable quotient
is detected by

\[
 \boxed{
 I(s)=\frac{a_5(s)^5}{a_3(s)a_4(s)^6}
     =\frac{g_5(s)^5g_1(s)^2}
            {g_3(s)g_4(s)^6}.
 }                                                        \tag{2.3}
\]

## 3. The invariant is nonconstant

At `s=0`, equation (2.2) gives

\[
 g_3(0)=a\ne0,\qquad g_4(s)=5s.
\]

Also

\[
 g_1(s)=5s^4+3as^2+2bs+c.
\]

Because `a!=0`, the order of `g_1` at zero is at most two.  Therefore the
order of (2.3) at zero is at most

\[
 2\cdot2-6=-2.
\]

Thus `I(s)` has a pole at zero and is nonconstant.

This also identifies the elementary exceptional presentation.  If
`P(T)=T^5+d`, then

\[
 I(s)=\frac1{6250}
\]

is constant.  Condition (1.3) excludes precisely this easy centered
pure-power obstruction; a non-affine change of primitive generator is not
needed after the initial choice of `eta`.

A nonconstant rational function over the infinite field `K` takes infinitely
many values on `K` away from its finite pole set.  Remove the finitely many
`s` for which

\[
 g_1(s)g_3(s)g_4(s)g_5(s)=0.
\]

The remaining translations lie in the clean coefficient torus, and
`I(s)` still takes infinitely many `K`-values.

## 4. The common abstract fiber

Let `F_s` be the quadratic-gauge map attached to `G_s`.  At

\[
 y_s=
 \left(1,0,-\frac{2P(s)}{g_1(s)}\right)                \tag{4.1}
\]

its inverse polynomial is

\[
 G_s(S)-\frac{g_1(s)}2
 \left(-\frac{2P(s)}{g_1(s)}\right)
 =P(s+S).                                              \tag{4.2}
\]

Since `eta` is primitive and `A` is finite etale, `P` is squarefree of
degree five.  The root-engineered reconstruction theorem therefore gives
the complete fiber

\[
 F_s^{-1}(y_s)
 \simeq
 \operatorname{Spec}K[S]/(P(s+S))
 \simeq\operatorname{Spec}A.                           \tag{4.3}
\]

The maps have determinant `-2`; a fixed target scaling makes every
determinant one without changing the stable class or the fiber.

If `I(s_1)!=I(s_2)`, the stable-moduli theorem proves that `F_{s_1}` and
`F_{s_2}` are not stably polynomially left--right equivalent.  The
infinitely many values from Section 3 therefore prove

\[
 |\mathcal R_K(A)|=\infty.
\]

## 5. Scope and the next degree

The argument is uniform over every characteristic-zero field, including
products and the split algebra.  It uses only:

1. a primitive generator with nonzero second trace moment;
2. translation of that generator;
3. the exact stable-orbit classification of the quadratic-gauge
   coefficient torus.

For degree four the quadratic-gauge coarse stable quotient is a point, which
is why the separate weighted trace-quadric proof is necessary.  For every
degree `N>=6`, the top-three-coefficient invariant gives the uniform
continuation proved in the
[all-rank multiplicity note](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md).

## 6. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_quintic_fiber_multiplicity.py
```

The checker verifies the translated derivative jets, the primitive
weight-lattice relation, the stable invariant (2.3), its pole forced by
`Tr(eta^2)!=0`, and the constant pure-power exception.
