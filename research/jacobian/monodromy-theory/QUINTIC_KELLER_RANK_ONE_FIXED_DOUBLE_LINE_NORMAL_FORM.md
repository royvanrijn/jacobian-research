# A normal form for the fixed double line in top rank one

**Theorem QuinticRankOneFixedDoubleLineNormalForm.** Let \(F\) be a
noninvertible characteristic-zero three-variable Keller map of degree
five and top rank one. Suppose its cancelling leading quartic pencil
has primitive degree one and the repeated line in its fixed factor has
multiplicity exactly two. Over the algebraic closure there are affine
source and target changes for which
\[
 Q=x+x^2(z+xy),\qquad
 R_4=x^2 B_2(x,y),\qquad B_2(0,y)=c y^2,\quad c\ne0,
 \qquad x\mid R_3,\qquad P_5\in k[x,y].
\]
Here \(k\) denotes the field in the chosen frame. In particular the
cancelling quartic with leading term divisible by \(x^3\) is unique
up to scalar and its leading term is \(x^3 y\) in this frame. The
normal form is a necessary condition.

## The exceptional divisor over the double line

Work over \(\mathbb C\). The
[primitive pencil and valuation restriction](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil)
makes \(P_5,Q_4,R_4\) binary in the two pencil directions. Take \(x=0\)
as the repeated line and use \(x,y\) as the binary coordinates. The two
cancelling quartics are independent. By the exact multiplicity two
assumption we can choose \(R_4=x^2 B_2(x,y)\) with
\(r=B_2(0,y)=c y^2\ne0\), and choose \(Q_4=x^3(a x+b y)\).
The [cubic restriction](QUINTIC_KELLER_RANK_ONE_LINE_PENCIL_CUBIC_RESTRICTION.md)
gives \(Q_3=xU,\ R_3=xV\) for homogeneous quadratics \(U,V\).

Restrict quadratics to the line \(x=0\) by a subscript \(0\), and
write \(D=(Q_2)_0,\ W=(R_2)_0\). In the projective fibre of
\(R+\lambda Q\), blow up the line \(x=w=0\), with \(w\) the
homogenizing coordinate. On the chart \(x=wX\) the exceptional
curve over the line has quadratic equation
\[
             r X^2+(V_0+\lambda U_0)X+W+\lambda D=0.
\]
Its discriminant is the binary quartic
\(\Delta_\lambda=(V_0+\lambda U_0)^2-4r(W+\lambda D)\).
The [exceptional divisor residue argument for a fixed double line](QUINTIC_KELLER_RANK_ONE_DOUBLE_LINE_CONIC_PENCIL.md#the-exceptional-discriminant)
depends only on these quadratic, cubic and quartic restrictions.
It applies with the present \(r\): if the discriminant is nonzero,
the residue on the double cover is a nonzero constant times
\(dy/\sqrt{\Delta_\lambda}\) in an affine coordinate of the base
line. Four simple roots give a holomorphic differential with a
nonzero period, and a double root with no higher root gives a
simple pole with nonzero residue on the normalization. Both
contradict exactness on the original affine fibre by a normal tube.
The same argument works for a split double cover. Thus the
discriminant vanishes identically in \(\lambda\), or has a triple
root for general \(\lambda\). In the latter case the triple root
is fixed. Indeed, writing the coefficient map as the product of
a triple linear factor and a linear factor gives pullback degree
\(3s_1+s_2\), where \(s_1,s_2\) are the degrees of the rational
root maps on the parameter line. The unique triple root and its
remaining factor are rational, and the degree of the coefficient
map is at most two. Hence \(s_1=0\).

When \(\Delta_\lambda=0\) identically, its quadratic and linear
coefficients in \(\lambda\) give \(U_0=D=0\). In the other case let
\(L\) be the fixed triple factor. First \(U_0=k L^2\). If \(L\)
does not divide \(r\), the linear coefficient gives \(L^2\mid D\).
If \(L\) divides \(r\), the constant coefficient first gives
\(L\mid V_0\), and the linear coefficient then gives
\(L^3\mid rD\), or \(L\mid D\) since \(r\) is proportional to \(L^2\).

## The distinguished quartic component

Take \(u=L\) if a triple factor is used, or \(u=y\) when \(U_0=D=0\).
Complete \(u\) to a binary frame \(u,v\) on the line \(x=0\). If
\(u\) is not proportional to \(y\), take \(v=y\), and otherwise choose
\(v\) independent of \(u=y\). In these affine coordinates the conditions
on \(Q_4,Q_3,Q_2\) give
\[
 Q=f(x)+A(x)u+(kx+d)u^2+(a_*(x)+e u)v,\qquad
 \deg f\le4,\quad \deg A,\deg a_*\le3.
 \tag{1}
\]
Here \(e\ne0\) is possible only for \(u=y\), in which case
\(\deg a_*\le2\). Its coefficient \(a_*+e u\) is squarefree
because its \(u\)-derivative is a nonzero constant. The
[reduced quadratic completion](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient)
excludes this case.

Assume \(e=0\). If \(a_*\ne0\), the
[parallel coefficient classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#coordinates-on-the-parallel-support-fibres)
gives a repeated root \(x=\xi\) of \(a_*\). At that root
\(Q|_{x=\xi}\) is constant and \(Q_x|_{x=\xi}\) is a nonzero
constant. It follows that \(kx+d\) is divisible by
\((x-\xi)^2\), so \(k=d=0\).

If \(a_*=0\), first suppose \(k\ne0\). Translate \(x\) to
make \(d=0\). If \(A(0)=0\), then \(Q_u=0\) on \(x=0\)
and its \(x\)-derivative is a quadratic in \(u\) with nonzero
leading coefficient there, giving a critical point. If \(A(0)\ne0\),
set \(u=-A(x)/(2kx)\) to solve \(Q_u=0\) off \(x=0\).
The equation for \(Q_x=0\) is the derivative of
\[
                      f(x)-\frac{A(x)^2}{4kx}
\]
equal to zero. This rational function has a simple pole at zero and
grows polynomially at infinity. Indeed, if \(\deg A=3\) its second
term has degree-five growth, and otherwise \(\deg f=4\) since \(Q\)
is quartic. Hence its derivative after multiplication by \(x^2\)
is a nonconstant polynomial with nonzero constant term. A nonzero
root produces a critical point. If instead \(k=0\) and \(d\ne0\), the
[constant quadratic coefficient completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component)
applies using the \(u\)-direction. Thus in the planar case
\(k=d=0\) also. The coefficient \(A\) is nonzero, since a
univariate quartic has a critical point, and the same parallel
coefficient classification applies with \(u\) as the linear
source variable.

Translate the repeated root \(\xi\) to zero and normalize the
constant and linear terms of \(Q\) there. The two conclusions of
the classification give
\[
        Q=x+x^2\bigl(f_2(x)+(B_0+x B_1)\cdot(y,z)\bigr),\qquad
        \deg f_2\le2 ,
 \tag{2}
\]
in the original binary directions \(x,y\), with constant covectors
\(B_0,B_1\) in the \(y,z\)-plane. The repeated-coordinate
translation preserves \(x\mid R_3\), because every first derivative
of \(R_4=x^2 B_2\) is divisible by \(x\). The covectors cannot both
vanish by submersion, and \(B_1\) is proportional to the covector
\(y\) by the leading form of \(Q\).

## Parallel covectors

Suppose \(B_0,B_1\) are dependent. If \(B_1\ne0\), write their
common direction as \(Y\). Its coefficient in \(Q\) is
\(x^2(a_1x+a_0)\) with \(a_1\ne0\). It cannot have a nonzero
simple root: at such a root \(Q_Y=0\) and one can choose \(Y\)
so that \(Q_x=0\). Therefore \(a_0=0\). If \(B_1=0\),
the nonzero coefficient is a scalar times \(x^2\).
In either case choose binary coordinates \(Y,Z\) with \(Y\)
the common direction and scale to write \(Q=g(x)+x^hY\) with
\(h=2\) or \(3\). This affine frame retains degree four for \(Q\).

The other components satisfy
\[
   D(P)R_Z-P_Z D(R)=\pm\kappa,\qquad
   D=-x^h\partial_x+(g'(x)+h x^{h-1}Y)\partial_Y .
\]
The rational constants of \(D\) are \(\mathbb C(Q)\); the
field \(\mathbb C(x,Y)=\mathbb C(Q,x)\) has \(D(x)=-x^h\).
There is no polynomial \(D\)-primitive of a nonzero constant.
A rational primitive of \(1\) is \(x^{1-h}/(h-1)\);
adding a rational function of \(Q\) cannot remove its pole.
In fact a reduced denominator in \(Q\) would have to divide
\(x^{h-1}\) by a univariate Bézout identity, impossible in
degree at least four. The residue and leading-coefficient
arguments of the
[exceptional scalar proof](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-through-degree-five)
therefore apply: both \(Z\)-degrees are positive, neither leading
coefficient is constant, and an affine pair in \(Z\) is impossible.
These three consequences only use the rational constants and the
absence of a polynomial primitive of a nonzero constant.

Let \(m=\deg_Z P,\ n=\deg_Z R\), with leading coefficients \(p_m,r_n\).
Then \(n\le2\) since \(R_4=x^2 B_2(x,y)\) and \(x\mid R_3\).
If \(n=1\), then \(\deg r_1\ge2\), for a nonzero linear or constant
coefficient would give reduced quadratic completion. The top
determinant coefficient gives \(p_m^n/r_n^m\in\mathbb C(Q)\).
For \(m\ge2\) both coefficients have degrees at most three, so the
ratio is constant: any finite zero or pole of a nonconstant rational
function of \(Q\) would force a full fibre polynomial \(Q-c\) to
divide one coefficient. These fibres are squarefree by submersion
and have degree four. In the case \(m=1,n=2\), the only possible
full fibre factor is in \(p_1\). If it occurs, reduce the bracket
modulo \(Q-c\), using \(D(Q)=0\). The reduced equation is
\(D(p_0)(2r_2 Z+r_1)=\pm\kappa\). On each irreducible component
of that fibre it forces \(r_2=0\), impossible since
\(\deg r_2\le2\). Thus the ratio is constant here too.

For \(n=1,\ m\ge2\), proportionality forces
\(2m\le m\deg r_1=\deg p_m\le5-m\), impossible. For \(n=2,\ m\ge3\)
it forces \(2\deg p_m=m\deg r_2\); the degree bound excludes this
when \(m\ge4\), and at \(m=3\) an odd right side or a right
side at least six excludes it. At \(m=n=2\), subtract a constant
multiple of \(R\) from \(P\) to reduce \(m\) to one, preserving
the degree bounds. The ratio for \(m=1,n=2\) makes \(r_2\)
proportional to \(p_1^2\), so subtract the corresponding multiple
of \(P^2\) from \(R\) to obtain an affine pair or a component
independent of \(Z\). The three degree-independent restrictions
exclude both outcomes. This removes all dependent covectors.

## The surviving affine frame

In the independent case choose \(y=B_1\cdot(y,z)\) and
\(z=B_0\cdot(y,z)\) in (2), using the original symbols on the
right. Absorb the \(x^2\) term of \(f_2\) into \(xy\) by adding
a multiple of \(x\) to the new \(y\), and absorb its affine
part into the new \(z\). These are affine changes and preserve
the binary leading forms and the factor \(x\) in \(R_3\).
This is exactly the asserted frame. All reductions over
\(\mathbb C\) apply in characteristic zero by embedding the
finitely generated coefficient field, with geometric changes
and polynomial inversion descending as usual.

## Exclusion of the surviving frame

**Theorem QuinticRankOneFixedDoubleLineExclusion.** A noninvertible
degree-five three-variable Keller map of top rank one and primitive
pencil degree one cannot have a fixed repeated line of multiplicity
exactly two.

Use the normal form just proved. Put \(Y=z+xy\), an automorphism
with determinant one, and regard \(y\) as the free variable for the
component \(Q=x+x^2Y\). On \(Q=t\) with \(x\ne0\), set
\(Y=T=t/x^2-1/x\). There the determinant equation reads, up to
changing the sign of \(\kappa\),
\[
                   (P_t)_x(R_t)_y-(P_t)_y(R_t)_x
                   =\kappa/x^2.                                      \tag{3}
\]
Before restricting \(Q=t\), the two other components satisfy the
exceptional cubic bracket in \(x,Y,y\). Thus their degrees in \(y\)
are positive and their leading coefficients in \(y\) are nonconstant,
by the degree-independent restrictions cited above. The total
degrees in \(y,z\) of terms of \(R\) are at most two by the form
of \(R_4\) and \(x\mid R_3\). Hence \(\deg_y R_t\le2\).
A term of \(P\) with \(y,z\)-degree five can only be a constant
multiple of \(y^5\), since \(P_5\) is binary in \(x,y\). It would
give a constant leading coefficient in \(y\) after the shear.
Hence \(\deg_y P_t\le4\). The passage from \(\mathbb C[x,Y]\) to
\(\mathbb C(t)[x,x^{-1}]\) is injective, so positive degrees and
nonzero coefficients persist after restriction to a generic fibre.

Expand
\[
\begin{aligned}
R_4&=a x^4+b x^3y+c x^2y^2, &&c\ne0,\\
R_3&=d x^3+e x^2y+f x^2z+g xy^2+h xyz+i xz^2,\\
R_2&=j x^2+k xy+l xz+m y^2+n yz+o z^2,\qquad
R_{\le1}=p x+q y+r z+s_0 .
\end{aligned}
\]
Writing \(R_t=A y^2+B y+C\), substitution \(z=T-xy\) gives
\[
\begin{aligned}
 A&=i x^3+(c-h+o)x^2+(g-n)x+m,\\
 B&=(b-f)x^3+(e-l)x^2+(k-r)x+q+((h-2o)x-2i x^2+n)T,\\
 C&=a x^4+d x^3+j x^2+p x+s_0+(f x^2+l x+r)T+(i x+o)T^2.
\end{aligned}                                                         \tag{4}
\]

### Affine case

Suppose \(A=0\). Here \(i=0\), \(c-h+o=0\), \(g=n\), \(m=0\).
The double-line discriminant of \(R\) in the argument above is
\[
          y^2\bigl((g y+h z)^2-4c(m y^2+n yz+o z^2)\bigr).
\]
It is unchanged by taking a pencil combination with \(Q\), since
\((Q_3/x)_0=(Q_2)_0=0\). It vanishes or has a triple root, forcing
the parenthesis to be divisible by \(y\). Thus \(4co=h^2\), and
the relation \(c-h+o=0\) yields \(h=2c,\ o=c\).
On \(R_t=\xi\), solve \(y=(\xi-C)/B\). Equation (3) says the
rational differential \(\kappa\,dx/(x^2 B)\) has a rational
primitive in \(x\). If \(n\ne0\), the polynomial
\[
 x^2B=(b-f)x^5+(e-l)x^4+(k-r)x^3+q x^2-nx+nt
\]
has generically simple nonzero roots. Indeed its derivative is
independent of \(t\) and has nonzero constant coefficient, and
its constant term is \(nt\). The resulting simple poles contradict
exactness. Hence \(n=0\). At \(x=0\), \(dQ=dx\) and
\(R|_{x=0}=c z^2+r z+q y+s_0\). If \(q=0\), choose \(z\) with
\(R_z=0\), violating the Keller condition. Thus \(q\ne0\), and
\(B\) is a polynomial of degree at most three with \(B(0)=q\).
A simple root of \(B\) gives a nonzero residue. If \(B\) is
nonconstant and all its roots are repeated, it is proportional
to \((x-\alpha)^d\) with \(d=2\) or \(3\), \(\alpha\ne0\).
Then \(B'(0)\ne0\) gives a nonzero residue at zero.
Finally a constant \(B\) would be a constant leading coefficient
of \(R\) in the free variable before restricting \(Q=t\), already
excluded. This rules out the affine case.

### Quadratic case

Suppose \(A\ne0\), which is nonconstant by the leading-coefficient
restriction. Put \(K=C-B^2/(4A)\) and consider \(R_t=\xi\) over
\(\mathbb C(t,\xi)\). With \(U=2Ay+B\), the fibre has
\(U^2=4A(\xi-K)\). It is quadratic over \(\mathbb C(t,\xi,x)\),
since the linear factor \(\xi-K\) has odd exponent.
The derivative of \(P_t\) at fixed \(\xi\) is
\(\kappa/(x^2 U)\) by (3). Its primitive's odd part under
\(U\mapsto-U\) is \(U H(x,\xi)\), with \(H\) polynomial in
\(\xi\) over \(\mathbb C(t,x)\). This follows by substituting
\(y=(U-B)/(2A)\) into \(P_t\) and reducing powers of \(U^2\).
The equation
\[
  2\bigl(4A(\xi-K)\bigr)H_x+\bigl(4A(\xi-K)\bigr)_x H
                           =2\kappa/x^2
\]
forces \(2A(H_{\mathrm{top}})_x+A' H_{\mathrm{top}}=0\) at the
highest degree in \(\xi\). Here \(H\ne0\) by the nonzero right
side. Therefore \(A H_{\mathrm{top}}^2\in\mathbb C(t)^\times\).
Every root of \(A\) has even multiplicity, so
\(A=a_2(x-\alpha)^2\) by \(\deg A\le3\) and nonconstancy.
Expand the identity \(d(UH)=\kappa\,dx/(x^2 U)\) formally
at \(\xi=\infty\) in half-integral powers. The coefficient of
\(\xi^{-1/2}\) on the right is a nonzero multiple of
\(dx/(x^2(x-\alpha))\), while on the left it is an exact
rational differential in \(x\). The residue at \(\alpha\ne0\)
would be nonzero. Hence \(\alpha=0\), and \(A=a_2 x^2\).

Take \(\lambda^2=a_2\), \(\lambda\ne0\), and set
\(v=\lambda xy+B/(2\lambda x)\), so \(R_t=v^2+K\).
Equation (3) becomes
\[
                   2v(P_t)_x-K_x(P_t)_v=\kappa/(\lambda x^3)
 \tag{5}
\]
in \(x,v\). The odd part of \(P_t\) in \(v\) has degree one or
three, since \(\deg_y P_t\le4\) and (5) has nonzero even
right side. If the odd part is \(L(x,t)v\), its leading
equation gives \(L_x=0\), and (5) gives
\(K=\beta(t)+\alpha_2(t)x^{-2}\) with \(\alpha_2\ne0\).
If the odd part is \(M(x,t)v^3+L(x,t)v\), (5) gives
\(M_x=0,\ L_x=3M K_x/2\). Its constant equation then makes
a nonconstant quadratic polynomial in \(K\) equal to an
affine polynomial in \(x^{-2}\) with nonzero slope. As
\(K\in\mathbb C(t)[x,x^{-1}]\), its only pole is then of
order one at \(x=0\), so \(K=\beta(t)+\alpha_1(t)x^{-1}\)
with \(\alpha_1\ne0\).

Compare these two possible pole profiles with (4).
The identity \(A=a_2x^2\) gives
\(i=m=0,\ g=n,\ a_2=c-h+o\).
The coefficient of \(x^{-6}\) in \(K\) is \(-n^2t^2/(4a_2)\),
so \(n=g=0\). Put \(J=h-2o\). Then the coefficients of
\(t^2x^{-4}\) and \(tx^{-3}\) are
\[
                   o-\frac{J^2}{4a_2},\qquad
                  -2o-\frac{J(q-J)}{2a_2}.
\]
Both vanish, giving \(4a_2 o=J^2\) and \(Jq=0\).
If \(J\ne0\) then \(q=0\), \(o\ne0\), and at \(x=0\)
the restriction \(R=o z^2+r z+s_0\) has a critical point
along that plane. Since \(dQ=dx\), this is impossible.
If \(J=0\) the first equality forces \(o=0\), and hence
\(h=0\). Thus \(a_2=c\), \(h=o=0\).

Now the coefficients of \(x^{-2}, x^{-1}\) in \(K\) are
\[
             r t-\frac{q^2}{4c},\qquad
             l t-r-\frac{q(k-r)}{2c}.
\]
Pole order one would make \(r=q=0\), so \(R|_{x=0}\)
would be constant, impossible. Pole order two therefore holds,
with \(l=0,\ r=-q(k-r)/(2c)\) and \(q\ne0\).
The odd part must be \(L(t)v\). It has
\(L(t)\in\mathbb C[t]\), since the substitution into \(v\)
expresses \(P_t\) with coefficients in
\(\mathbb C[t,x,x^{-1}]\). By (5),
\(L(t)(rt-q^2/(4c))\) is a nonzero constant, hence
\(r=0,\ k=0\), and \(L\) is a nonzero constant.

The even part of \(P_t\) is a polynomial in \(R_t\) with coefficients
in \(\mathbb C[t]\). To see this, take its top coefficient in \(v\)
in the homogeneous equation (5) with zero right side. It is
independent of \(x\), hence in \(\mathbb C[t]\), and one can
subtract that multiple of the corresponding power of \(R_t\).
Descend on its even degree. After replacing \(t\) by \(Q\), this
even part is polynomial in the original variables. But the odd
part \(L v\) has a nonzero simple pole at \(x=0\), because
\(B=(b-f)x^3+e x^2+q\), \(\lambda\ne0\) and \(q\ne0\).
It cannot occur in the original polynomial \(P\). This contradiction
finishes the quadratic case and the proof.

## Consequence for quintic inversion

**Corollary QuinticKellerTopRankTwoNonreducedBoundary.** A noninvertible
characteristic-zero three-variable Keller map of degree five has
degree-five parts of rank exactly two. The unique quartic component
up to scalar obtained by cancelling degree five has geometrically
nonreduced leading form.

The [preceding quartic-boundary reduction](QUINTIC_KELLER_RANK_ONE_CONIC_PENCIL.md#consequence-for-the-quintic-inversion-problem)
leaves top rank one or two, and in rank one reduces the primitive
pencil to degree one with a fixed repeated line of multiplicity two
or three. The double-line exclusion above and the
[triple-line exclusion](QUINTIC_KELLER_TRIPLE_FIXED_LINE_PENCIL_EXCLUSION.md)
eliminate both. The preceding quartic-boundary reduction gives
the asserted necessary condition in top rank two.
