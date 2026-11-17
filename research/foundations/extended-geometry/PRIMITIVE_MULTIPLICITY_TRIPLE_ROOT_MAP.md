# Primitive-multiplicity triple-root Cox map

The affine-source triple-root construction has a simpler target than the
five-factor square orientation.  The primitive source character `u=ac`
occurs with multiplicity two in the determinant ledger, leading to the
four-branch hypersurface

\[
 u^2v=x\,y(x+1)(x-y+1).
\]

The resulting map has constant hypersurface-residue Jacobian `-1`, is
birational, and has three distinct dicritical target components.  Within
the natural class of monomial factor-splitting targets `U^mV=f(x,y)`, this
choice is unique.

## 1. Minimal coordinates

On the triple-root affine completion `A^3_(a,b,c)`, put

\[
 x=bc,\qquad y=a^2c,\qquad
 \tau=x-y+1.                                        \tag{1}
\]

The simplified five-branch target from the
[square-oriented construction](AFFINE_SOURCE_TRIPLE_ROOT_COX_MAP.md) is

\[
 S^2=C\,x\,y(x+1)\tau,
\]

with

\[
 C=b(x+1)\tau,\qquad
 S=ax(x+1)\tau.                                     \tag{2}
\]

The common primitive character hidden in (2) is

\[
 u=ac.
\]

Define

\[
 v=b(x+1)\tau.                                      \tag{3}
\]

Then

\[
 u^2v
 =a^2c^2b(x+1)\tau
 =x\,y(x+1)\tau.                                    \tag{4}
\]

## 2. The primitive-multiplicity target

Let

\[
 T_{\rm pm}
 =
 \{u^2v=x\,y(x+1)(x-y+1)\}
 \subset\mathbb A^4_{x,y,u,v}.                      \tag{5}
\]

Equations (1) and (3) give a polynomial morphism

\[
\boxed{
\Phi_{\rm pm}:\mathbb A^3_{a,b,c}\longrightarrow T_{\rm pm},
}
                                                            \tag{6}
\]

\[
(a,b,c)\longmapsto
\left(
bc,\ a^2c,\ ac,\ b(bc+1)(bc-a^2c+1)
\right).
\]

This target is an affine modification of `A^3_(x,y,u)` along the
coefficient `u^2`.  Its singular locus consists of three affine `v`-lines
over the intersection points

\[
(x,y)=(0,0),\quad(0,1),\quad(-1,0)                  \tag{7}
\]

of the four-line arrangement.  The target is smooth at the generic point
of each of its four boundary components.  A principal open can exclude
(7) while retaining all four generic points.

## 3. Constant residue Jacobian

For

\[
 G=u^2v-x\,y(x+1)(x-y+1),
\]

one has

\[
\frac{\partial G}{\partial v}=u^2.
\]

Use the hypersurface residue form

\[
\Omega_{\rm pm}
=\frac{dx\wedge dy\wedge du}{u^2}.                  \tag{8}
\]

Direct differentiation gives

\[
\det\frac{\partial(x,y,u)}{\partial(a,b,c)}
=-a^2c^2=-u^2.                                      \tag{9}
\]

Consequently

\[
\boxed{
\Phi_{\rm pm}^*\Omega_{\rm pm}
=-da\wedge db\wedge dc.
}                                                     \tag{10}
\]

Thus (6) has constant residue Jacobian `-1` on every smooth principal
target chart.

## 4. Birational inverse and dicritical components

On `uy!=0`, the inverse is

\[
\boxed{
a=\frac yu,\qquad
c=\frac{u^2}{y},\qquad
b=\frac{xy}{u^2}.
}                                                     \tag{11}
\]

Hence (6) is birational.

Over `u=0`, equation (5) has four irreducible components, corresponding to

\[
y=0,\qquad x=0,\qquad x+1=0,\qquad x-y+1=0.         \tag{12}
\]

At a generic component, `u` is a uniformizer and the selected line factor
has order two.  Formula (11) yields:

\[
\begin{array}{c|ccc|c}
\text{component}&v(a)&v(b)&v(c)&\text{behavior}\\ \hline
y=0&1&0&0&\text{finite: }a=0,\\
x=0&-1&0&2&\text{dicritical},\\
x+1=0&-1&-2&2&\text{dicritical},\\
x-y+1=0&-1&-2&2&\text{dicritical}.
\end{array}                                          \tag{13}
\]

The three dicritical components have distinct target images.

## 5. Uniqueness in the monomial factor-splitting class

Let

\[
 f=x\,y(x+1)(x-y+1)
 =a^2bc^2(bc+1)(bc-a^2c+1).
\]

Choose a monomial divisor `U` of the displayed five source factors and put
`V=f/U`.  A target

\[
 U^mV=f(x,y)
\]

has residue denominator `U^m`.  Exhausting all divisor exponents shows:

> The ratio
> \[
> \frac{1}{U}
> \det\frac{\partial(x,y,U)}{\partial(a,b,c)}
> \]
> is always an integer multiple of `ac`, whenever it is nonzero.

The only choice for which the determinant is a constant times an integral
power `U^m` is

\[
 U=ac,\qquad m=2,                                   \tag{14}
\]

up to multiplication by a nonzero constant and exchanging presentation
conventions.  Equation (5) is therefore the unique constant-residue target
within this natural factor-splitting class.

## 6. Finite fields

For every finite field `F_q`, the four-line arrangement has

\[
4(q-1)
\]

rational points.  Counting (5) by `u` gives

\[
\boxed{
\#T_{\rm pm}(\mathbb F_q)
=q^2(q-1)+4q(q-1)
=q(q-1)(q+4).
}                                                     \tag{15}
\]

This differs from the source count `q^3` by `q(3q-4)`.  Therefore (6)
cannot be a permutation of full rational point sets over any finite field.

The simplified target improves the boundary ledger and singular locus, but
does not solve the arithmetic gap.

Its full resolution is also explicit.  The
[affine-modification chain](TRIPLE_ROOT_AFFINE_MODIFICATION_CHAIN.md)
successively adjoins `a=y/u`, `c=u/a`, and `b=x/c`.  Every step preserves
the residue form; the target becomes `A^3` only at the final step, when all
dicritical boundary information has been filled.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_primitive_multiplicity_triple_root_map.py
```

The checker verifies the target equation, residue Jacobian, inverse,
singular locus, dicritical valuations, finite-field count, and uniqueness
among all monomial divisors of `f`.
