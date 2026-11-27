# Pure quartic and binary triple-line Keller components

**Theorem QuinticKellerPureQuarticComponentCompletion.** A
three-variable Keller map of degree at most five in characteristic
zero is invertible if a component has fourth-degree part equal
to a nonzero fourth power of a linear form.

**Proposition QuinticKellerBinaryTripleLineScalar.** Suppose a
three-variable Keller map of degree at most five is noninvertible
and has a quartic component with leading form \(x^3 y\) in independent
linear coordinates. After affine source changes and affine rescaling
of that component the scalar has the form
\[
                             R=x+x^2(z+xy).
\]
In particular, in the primitive degree-one quintic pencil, where
all three degree-four and degree-five forms are binary, this is the
only remaining possibility for the cancelling quartic.

## The cubic on the repeated line

Work over \(\mathbb C\) by embedding a finitely generated coefficient
field. Suppose the map is noninvertible. For leading quartic \(x^4\), the
[pure-quartic fibre argument](KELLER_PURE_QUARTIC_NONZERO_BINARY_CUBIC.md)
forces the cubic to be divisible by \(x\). The same holds for \(x^3y\)
by the [binary quartic cubic restriction](QUINTIC_KELLER_BINARY_QUARTIC_REPEATED_CUBIC_VANISHING.md),
since \(x=0\) is its unique repeated line.

Let \(U=(R_3/x)|_{x=0}\). The
[exceptional section argument for a triple line](QUINTIC_KELLER_TRIPLE_FIXED_LINE_PENCIL_EXCLUSION.md#a-square-restriction-from-the-line-blowup)
applies to this single component, also when the leading quartic is
\(x^4\): on blowing up \(x=w=0\) in the chart \(w=xW\), the
leading quartic is of order at least three and the exceptional
section \(W=0\) has residue proportional to \(dy/U(y,1)\).
If \(U\ne0\), the tube periods therefore force \(U\) to be a square.

Assume \(U\ne0\). Choose coordinates \(u,v\) complementary to \(x\)
so \(U\) is a nonzero scalar times \(u^2\).
If \(R_4=x^3y\) and \(u\) is independent of \(y\), take \(v=y\).
If \(u\) is proportional to \(y\), take \(u=y\). In the case
\(R_4=x^4\) any complementary \(v\) works.
The component is quadratic in \(v\) with constant quadratic coefficient,
coming only from \(R_2\). A nonzero coefficient gives inversion by the
[quadratic-source theorem](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
In the remaining affine-linear case the coefficient of \(v\) has
cubic part either zero or a nonzero multiple of \(x^3\).
The [cubic coefficient classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
therefore makes it zero or gives the repeated-parallel residual
\(R-\mu=\lambda s+s^2H\). In the latter case the linear direction
of \(s\) is the repeated line \(x=0\), since its square divides
\(R_4\) and that line is unique. But the coefficient of \(u^2\)
in \(R\) is \(kx+d\) with \(k\ne0\), and this cannot be divisible
by the square of the affine coordinate \(s\). Hence the
coefficient of \(v\) must be zero.

In the resulting planar case \(R_4\) is independent of \(v\).
Thus \(R=f(x)+A(x)u+(kx+d)u^2\), where \(k\ne0\),
\(\deg f\le4,\deg A\le3\), and either \(\deg A=3\) or
\(\deg f=4\) with \(\deg A\le2\).
Translate \(x\) to eliminate \(d\). If \(A(0)=0\), then
\(R_u=0\) on \(x=0\) and the quadratic \(R_x(0,u)\)
has a root. If \(A(0)\ne0\), solve \(R_u=0\) at \(x\ne0\)
by \(u=-A(x)/(2kx)\).
The stationary value \(f(x)-A(x)^2/(4kx)\) has a simple
pole at zero and nonconstant growth at infinity: its
second term has degree-five growth if \(\deg A=3\),
and otherwise the quartic term of \(f\) dominates.
The derivative multiplied by \(x^2\) is a nonconstant
polynomial nonzero at zero, so it has a nonzero root.
This is again a critical point. Therefore \(U=0\)
and \(x^2\mid R_3\).

## The remaining scalar forms

Return to a frame where \(R_4=x^4\) or \(R_4=x^3 y\).
A nonzero \(z^2\) coefficient in \(R_2\) gives quadratic-source
completion, since \(x^2\mid R_3\) and \(R_4\) is independent
of \(z\). Otherwise the component is affine linear in \(z\),
with coefficient of degree at most two. If the coefficient
is nonzero, the [linear-source classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
forces the same repeated-parallel form with square in the
\(x\)-direction. By an affine change and rescaling we can put
it in one of the forms
\[
     R=x+x^2(z+a y+f(x))\quad(R_4=x^4),\qquad
     R=x+x^2(z+xy)\quad(R_4=x^3y).
 \tag{1}
\]
For the second form the quadratic part of \(H\) is \(xy\);
its lower affine part can be absorbed into \(z\).
For the first, \(f\) has degree two, since \(R_4=x^4\).

If the coefficient of \(z\) is zero, the component is planar.
In that case any \(y^2\) term in \(R_2\) again gives quadratic
completion. Thus \(R=f(x)+A(x)y\). The polynomial \(A\) is
nonzero because \(R\) is a quartic submersion.
A nonzero constant gives a bounded graph completion, and a
simple root would give a critical point by solving \(R_x=0\)
there. In the pure fourth-power case \(\deg A\le2\), so
\(A\) is a repeated square. In the triple-line case \(\deg A=3\),
so \(A\) is a cube. Translate its root to zero, use
\(f'(0)\ne0\) there, and take an affine change in \(y\) and target
scaling. This gives \(R=x+x^2(y+f_2(x))\) with \(\deg f_2=2\)
in the pure fourth-power case, or \(R=x+d x^2+x^3y\) in
the triple-line case.

## Degree-four exceptional scalars

The first case of (1) becomes a planar scalar by an affine change
in the \(y,z\)-plane. All the planar scalars just obtained have
the form \(r=g(x)+x^j V\), of ordinary degree four, with \(j=2\)
or \(3\) and \(g'(0)\ne0\). They admit no quintic Keller mates.
Indeed \(D=-x^j\partial_x+(\partial_x r)\partial_V\) has rational
constants \(\mathbb C(r)\). A polynomial primitive of a constant
would differ from a multiple of \(x^{1-j}/(j-1)\) by a rational
function in \(r\). Its denominator after substitution would
divide \(x^{j-1}\), which has degree smaller than \(\deg r=4\).
The [exceptional bracket restrictions](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-through-degree-five)
apply to the mates \(P,Q\) in the free variable, with degrees
\(m,n\) and leading coefficients \(p_m,q_n\).
Equal degrees give a constant ratio and can be canceled by
constant subtraction when \(m,n\ge2\), since neither coefficient
has degree four. An affine pair is already forbidden.

For unequal degrees \(m>n\) with \(m\ge3\), the leading coefficient
\(p_m\) has degree at most two. If \(q_n\) has a full fibre factor,
it is a scalar multiple of that factor by its bound at most four.
Then the leading relation \(p_m^n/q_n^m\in\mathbb C(r)\) would
make \(p_m\) constant. Otherwise the ratio is constant because no
coefficient contains a full fibre factor. The degree relation
\(n\deg p_m=m\deg q_n\) is impossible: \(m=3\) requires
\(\deg p_m\ge3\), and \(m=4\) requires \(\deg p_m\ge2\)
though the bound is one. Degree \(m\ge5\) would have a constant
leading coefficient. For \(m=2,n=1\), a nonconstant ratio
requires a full fibre factor in \(q_1\), since \(\deg p_2\le3\).
Reducing the bracket modulo that fibre forces \(p_2\) to contain
it as well. Writing the free variable as \(y\), the reduced
bracket equation is \(-\overline{P_y}\,\overline{D(Q|_{y=0})}=\kappa\).
The second factor is independent of \(y\), so the coefficient
of \(y\) in \(\overline{P_y}\) vanishes on each component. This is impossible in degree three for
a degree-four fibre. A constant ratio permits subtraction of a
multiple of \(Q^2\), giving an affine or independent mate.
Thus the planar scalars are excluded.

The pure fourth power case is now exhausted. In the binary
primitive degree-one quintic pencil the cancelling quartic has
a repeated line by the
[double-conic closure](QUINTIC_KELLER_RANK_TWO_PURE_FIFTH_TRANSVERSE_DOUBLE_CONIC.md#exclusion-of-the-double-conic-boundary).
A repeated line with zero cubic restriction exists by the
binary quartic result. Multiplicity two is excluded by the
[planar](QUINTIC_KELLER_BINARY_QUARTIC_REPEATED_CUBIC_VANISHING.md#excluding-the-planar-multiplicity-two-case)
and [affine scalar](QUINTIC_KELLER_BINARY_REPEATED_AFFINE_SCALAR_EXCLUSION.md)
arguments. Multiplicity four is the pure quartic just
excluded. Thus the cancelling quartic has a triple line,
and its only remaining scalar in this reduction is the second form in (1).
The [binary triple-line theorem](QUINTIC_KELLER_BINARY_TRIPLE_LINE_SCALAR_EXCLUSION.md)
excludes that scalar with the independent binary quintic leading forms.
