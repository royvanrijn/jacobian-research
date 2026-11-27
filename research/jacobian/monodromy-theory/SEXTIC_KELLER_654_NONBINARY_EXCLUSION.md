# Exclusion of the nonbinary reduced \((6,5,4)\) profile

**Theorem SexticKeller654NonbinaryReducedPairExclusion.** Let
\(F=(F_6,F_5,F_4)\) be a characteristic-zero three-variable Keller map
with the indicated component degrees. Suppose the highest forms of
\(F_5,F_4\) are coprime and squarefree over the algebraic closure. Then
the three highest forms belong to a polynomial ring in the same two
linear source forms.

Work over an algebraically closed field \(k\). By the
[reduced lower pencil and first correction theorems](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md),
it suffices to exclude a Keller completion of the leading tuple
\[
 C=x^4\big((y/x)^j(z/x)+b(y/x)\big),\qquad
 B=x(C+\lambda x^4),\qquad A=x^2C+\tfrac53\lambda x^6,
\]
where \(j=2\) or \(3\), \(b(0)=-5\lambda\), and \(\lambda\ne0\). Write
\[
 h=C+5\lambda x^4=y^m H,\qquad 1\le m\le j.
\]
Here \(H\) is irreducible and coprime to \(y\). Indeed, the affine level
equation after division by \(y^m\) is linear in \(z\) with coprime
coefficient and constant term, by the maximality of \(m\) and \(m\le j\).
Its homogenization has degree \(4-m\): for \(j=3\) this follows from
the \(z\) term, and for \(j=2\) the two-place form has a nonzero
\(y^4\) coefficient. In particular \(H\) is also coprime to \(x\).

## The first error and two charts

Start with reference polynomials \(G_A(x,F_4),G_B(x,F_4)\) with the
given leading terms, for instance \(x^2F_4+\frac53\lambda x^6\) and
\(xF_4+\lambda x^5\). Give their arguments weights \(1,4\). Let \(r\)
be the first degree loss in the errors and \(E_A,E_B\) the errors
in degrees \(6-r,5-r\). A missing part is zero. The determinant and
critical-level arguments in the cited descent give
\(E_A-2xE_B\in k[x,C]\), and \(E_B(1,y,z)\) is constant on the open
part of \(H(1,y,z)=0\) where \(y\ne0\). Subtract these constants in
the following way. Absorb \(E_A-2xE_B\) into the sextic reference
by evaluating the homogeneous polynomial expression at \(x,F_4\).
If the value of \(E_B(1,y,z)\) on \(H\) is \(c\), add \(c x^{5-r}\)
to the quintic reference and \(2c x^{6-r}\) to the sextic
reference. The resulting error parts satisfy
\[
 E_A=2xE_B,\qquad H\mid E_B.
 \tag{1}
\]
If \(r=6\), only an additive constant in the sextic component is
left and it can be absorbed immediately. For \(r\le5\), the
preceding adjustments either increase \(r\) or leave a nonzero
\(E_B\) with \(5-r\ge4-m\). Thus it suffices to consider
\(r\le m+1\) and to prove \(y^m\mid E_B\). Then \(h\mid E_B\)
and \(\deg E_B\le4\), so both errors in (1) can be absorbed
into the reference. Iteration would factor the whole map through
\(x,F_4\), contradicting the Keller condition.

Rescale the components by \(F_{d,\epsilon}=\epsilon^d
F_d(x/\epsilon,y/\epsilon,z/\epsilon)\). Their determinant is
a nonzero constant times \(\epsilon^{12}\). Fix the rescaled
quartic value to an indeterminate \(S\), and let \(a,b\) denote
the rescaled references at that value. At \(\epsilon=0\),
\(b=xS+\lambda x^5\), \(a=x^2S+\tfrac53\lambda x^6\).
Take a root \(x_0\) of \(S=-5\lambda x_0^4\), and work over a
finite extension of \(k(S)\) containing \(x_0\) and a square
root of the nonzero quadratic coefficient at the critical point.
The simple zero of \(b_x\) lifts to a critical point
\(\xi(\epsilon)\). A formal square root of
\((b(x)-b(\xi))/(x-\xi)^2\) gives a parameter \(s\) with
\(b(x)=b(\xi)+s^2\) and \(x_s(0,0)\ne0\).

Put \(U=F_{5,\epsilon}-b\). Subtract from \(F_{6,\epsilon}\) a
multiple \(\beta(\epsilon)F_{5,\epsilon}\) and a constant in
\(s\) so that the same operations on the reference remove its
constant and \(s^2\) coefficients. Call the resulting functions
\(Q,Q_{\mathrm{ref}}\). Then
\[
 Q_{\mathrm{ref}}=\alpha(\epsilon)s+\gamma(\epsilon)s^3+O(s^4),
 \qquad \alpha(0)=0,\quad \beta(0)=2x_0,\quad
 \gamma(0)=2d/3,\quad d=2x_s(0,0)\ne0.
\]
The last identities follow from \(a_x(x,0)=2x b_x(x,0)\).
Both \(U\) and \(Q-Q_{\mathrm{ref}}\) have order at least \(r\)
in \(\epsilon\). Denote coefficients in \(\epsilon\) by
subscripts. At order \(r\), (1) gives
\[
 (Q-Q_{\mathrm{ref}})_r=2(x(s,0)-x_0)U_r.
 \tag{2}
\]

There is a series \(R(s,\epsilon)\) with scalar Laurent series
coefficients in \(s\), satisfying
\[
 Q=\sum_{n\ge0}\frac{U^n}{n!}\,\delta^n R
 \pmod{\epsilon^{12}},\qquad \delta=\frac{1}{2s}\partial_s.
 \tag{3}
\]
To see this, use \(x,y,S\) as generic coordinates for the
rescaled source. Solving for \(z\) in the quartic equation has
rational coefficients at each order of \(\epsilon\), since
the leading quartic is linear in \(z\). Replace \(x\) by a
formal coordinate \(X\) such that \(F_{5,\epsilon}=b(X)\).
The inversion in \(\epsilon\) uses \(b_x(x,0)\ne0\) in the
generic field. The Keller determinant shows that at fixed
\(X,S\) the \(y\)-derivative of the sextic component vanishes
modulo \(\epsilon^{12}\), since both coordinate changes have
Jacobian units in that field. Its coefficients thus lie in
\(k(X,S)\). Making the subtractions defining \(Q\) and expanding
at the reference critical point gives \(R\).
In the original coordinates the new parameter has square
\(s^2+U\), which gives (3) by formal Taylor expansion. All
comparisons below use orders in \(\epsilon\) strictly less
than \(12\). By (2) and \(\delta R_0=2(x(s,0)-x_0)\) we have
\(R_i=(Q_{\mathrm{ref}})_i\) for \(i\le r\).

Use (3) at the two generic points of \(h=0\) on the leading
level \(C=S\) at \(x=x_0\). On \(H\) take \(y\) as the free
parameter and solve the quartic for \(z\), with derivative
a unit there. The series \(U,Q\) are regular in \(s,\epsilon\),
and \(U_r\) is divisible by \(s\). On \(y=0\) take \(z\)
as the free parameter and use formal Weierstrass preparation
in the other source coordinate. The resulting ring is free
over the regular series ring in \(s,\epsilon\) with basis
\(1,Y,\ldots,Y^{m-1}\), and relation
\[
 Y^m=s+\epsilon\sum_{i=0}^{m-1}b_i(s,\epsilon)Y^i.
 \tag{4}
\]
In fact at \(\epsilon=0\) the derivative of \(C-S\) in \(s\)
at \(y=0\) is nonzero, and \(s\) has order exactly \(m\)
in \(y\), since \(H(x_0,0,z)\ne0\) generically. A unit change
of \(y\) makes \(s=Y^m\) on this leading level, after adjoining
a root of the leading unit in the coefficient field if needed.
Lift that change to the formal chart and apply preparation
over the series ring in \(s,\epsilon\). The functions \(U,Q\)
are regular in this free module, whereas (3) holds after
inverting \(s\) and truncating in \(\epsilon\). The generic
inversions used to obtain (3) remain valid there, since
modulo \(\epsilon\) both \(y\) and \(b_x\) are units after
inverting \(s\). Write \(u=U_r\) in this chart. The
order of \(u\) in \(Y\) equals the multiplicity of \(y\)
in \(E_B\) if this multiplicity is less than \(m\):
after factoring out that power of \(y\), the constant term
at \(x=x_0,y=0\) is nonzero for generic \(z,x_0\).
Thus the desired conclusion is \(u\in(s)\) modulo \(\epsilon\).

## The first pole comparisons

If the order \(p\) of \(\alpha\) is less than \(r\), inspect
order \(p+r\) in (3). The only term that can have a pole besides
\(R_{p+r}\) is \(\alpha_p U_r/(2s)\). It is regular on \(H\),
so \(R_{p+r}\) is regular. On the other chart this implies
\(u\in(s)\).

Suppose \(p\ge r\), including \(\alpha=0\). On \(H\) all
\(R_i\) with \(i\le2r\) are regular. Indeed derivatives of
\(R_i\) for \(i<r\) use a reference with zero linear
coefficient at those orders, while at order \(2r\) the possible
poles from \(\delta R_r\) and \(\delta^2 R_0\) are simple
and have a factor \(U_r\). At order \(2r\) in the other chart,
the simple-pole numerator gives
\[
 u(u+2\alpha_r/d)\in(s).
 \tag{5}
\]
If \(r=4\), necessarily \(m=j=3\) and \(E_B\) is a nonzero
constant times \(H\). In that case \(u\) at \(Y=0\) varies
linearly with the free parameter \(z\), contrary to (5).
We can now assume \(r\le3\).

On \(H\), \(R_{3r}\) has at most a simple pole. For the terms
with \(n=1\) in (3) use regularity through order \(2r\).
For \(n=2\) a pole of order greater than one requires the
linear coefficient of \(R_r\) and the factor \(U_r^2\).
For \(n=3\) the term uses \(R_0\) and has a factor \(U_r^3\).
The factors compensate the higher poles on \(H\). Thus in
the other chart, at order \(3r\) the only terms with
denominators beyond \(s\) are the cubic-pole parts
\[
 -\frac{\alpha_r u^2}{8s^3}-\frac{d u^3}{24s^3}.
\]
The first comes from the linear term of \(R_r\) and the
second from \((2d/3)s^3\) in \(R_0\). In both terms the
products are at their lowest possible \(\epsilon\)-order,
so they are taken in the leading ring \(Y^m=s\).
All other terms have at most a simple pole even after reduction
in the free basis (4). It follows that
\[
 u^2(3\alpha_r+d u)\in(s^2).
 \tag{6}
\]
Equations (5) and (6) at \(Y=0\) force \(u(0)=0\).
If \(\alpha_r\ne0\), (5) then gives \(u\in(s)\).
If \(\alpha_r=0\), (6) gives \(\operatorname{ord}_Y u\ge
\lceil 2m/3\rceil\). This suffices for \(m=1,2\).
The remaining possibility has \(m=3\), \(\alpha_r=0\) and
\(\operatorname{ord}_Y u=2\). Then \(y^2H\mid E_B\) and
\(5-r\ge3\), so \(r=1\) or \(2\).

## The triple line at loss one

Suppose \(m=3\), \(r=1\). On \(H\), \(R_3\) is regular
and \(R_4\) has at most a simple pole. For the former, the
possible pole from \(\delta R_2\) has a factor \(U_1\), and
the \(n=2,3\) terms either have no pole or have enough
factors of \(U_1\) to cancel it. More explicitly, \(\delta^2
R_i\) for \(i=0,1\) has at worst a simple pole, and
\(\delta^3 R_0\) has at worst a cubic pole. At order four the
new higher poles are \(\delta^2 R_2\) times \(U_1^2\),
\(\delta^3 R_0\) times \(U_1^2 U_2\), \(\delta^3 R_1\) times
\(U_1^3\), and \(\delta^4 R_0\) times \(U_1^4\). Their pole
orders before multiplication are at most \(3,3,3,5\),
respectively.

In the free module (4), using a correction term of
\(\epsilon\)-degree \(k\ge1\) instead of \(s\) in reducing \(Y^3\)
can lower the \(Y\)-order by at most \(3k\), with \(s\) of
order three. At order four in (3) the term
\((U^4/4!)\delta^4 R_0\) has \(Y\)-order \(-7\) from
\(\gamma(0)\binom{3/2}{4}u^4/s^5\).
The term \(R_4\) and the \(n=1\) terms have order at least
\(-3\). For \(n=2\) only the linear part of \(R_2\) can
give more than a simple pole, with order at least \(4-9=-5\).
For \(n=3\) the term using \(R_1\) has order at least
\(6-9=-3\). The terms using \(R_0\) have order at least
\(4-9=-5\) from \(u^2U_2\) or \(6-3-9=-6\) from a
first-order correction in reducing \(u^3\). Thus the term
of order \(-7\) cannot cancel. This contradicts regularity
of \(Q\).

## The triple line at loss two

Now \(m=3\), \(r=2\), \(\alpha_2=0\). In (4) write
\(u=a_2Y^2+O(s)\) in the free basis, with \(a_2\ne0\).
Set \(L_3=[s]R_3\). At order five on \(H\), \(R_5\) is
regular, as every potential pole from a Taylor term has
a factor \(U_2\) and is at most simple. On \(y=0\), compare
the coefficient of \(Y\) in the simple-pole numerator
modulo \(s\) at that order. Only the first-order reduction
of \(u^2\) in the \(\delta^2 R_0\) term can contribute:
the \(\delta R_3\) term and the mixed product \(uU_3\)
are multiples of \(Y^2\), and \(u^2\) itself vanishes
modulo \(s\) at order zero in the multiplication.
Reducing \(a_2^2Y^4\) by (4) gives the
coefficient \(d a_2^2 b_0(0,0)/4\) of \(Y\).
Hence \(b_0(0,0)=0\).

At order six on \(H\), the residue of the simple pole
of \(R_6\) is
\[
 -L_3 v_H/2-d v_H^2/4,\qquad
 v_H=U_3(0)\ \text{on }H.
\]
The corresponding Taylor terms are \(\delta R_3\,U_3\)
and \(\frac12\delta^2 R_0\,U_3^2\).
All other possible poles have enough factors of \(U_2\)
to cancel them, since \(\alpha_i=0\) for \(i\le2\).
The residue of \(R_6\) is independent of the free parameter
\(y\) on \(H\). Therefore the rational function \(v_H\)
in \(y\) is constant over the scalar coefficient field.
For this local comparison only, add \(\epsilon^3v_H\) to the
reference \(b\) to make \(U_3(0)=0\) on \(H\).
This shift changes neither \(s\) nor (4). With this shifted
reference let \(v=U_3(0,0)\) on \(y=0\), where the
arguments are \(s,Y\), and let \(L=[s]R_3\).

Give \(\epsilon,s,Y\) weights \(2,3,1\), respectively, and
work modulo \(\epsilon^9\). In the chart on \(H\) the
weights of \(U\) are at least seven: \(U_2\) vanishes at
\(s=0\), \(U_3(0)=0\) by the shift, and later terms have
weight at least eight. Thus each nonidentity Taylor term
in (3) strictly raises weight, as \(\delta\) lowers it by
six. The only \(s\)-dependent terms of \(Q\) of weight at
most nine there are \((2d/3)s^3\) and its order-three-in-
\(\epsilon\) linear term. Indeed the reference has no constant
or quadratic term and \(\alpha_i=0\) for \(i\le2\);
the order-two error (2) vanishes to order two in \(s\).
Comparison in increasing weight shows that \(R\) has the
same terms through weight nine up to scalar constants in \(s\):
constants have zero Taylor derivatives, and any lower-weight
nonconstant term of \(R\) would survive at its lowest weight.
In particular the linear term is \(L\epsilon^3s\).

In the chart on \(y=0\), the initial relation is
\(Y^3=s+\eta\epsilon Y\), where \(\eta=b_1(0,0)\).
Indeed \(b_0(0,0)=0\), so all omitted terms in (4) have
weight above three. The initial part of \(U\) has weight six
and is \(w=a_2\epsilon^2Y^2+v\epsilon^3\).
Thus (3) in weight nine, apart from scalar constants in
\(s\), is the reduction of
\[
 \frac{2d}{3}s^3\sum_{n=0}^{4}\binom{3/2}{n}(w/s^2)^n
 + L\epsilon^3 s\sum_{n=0}^{2}\binom{1/2}{n}(w/s^2)^n
 \pmod{\epsilon^9}
 \tag{7}
\]
by that initial relation. The weight comparisons are made
on the coefficients in the fixed free basis. Reduction
preserves or increases weight, and never decreases
\(\epsilon\)-order, so the higher-weight terms of the
relation and the discarded orders cannot affect (7).

In (7) the coefficients of \(\epsilon^5Y^2/s\),
\(\epsilon^6/s\), \(\epsilon^7Y/s^2\), and
\(\epsilon^8Y^2/s^3\), respectively, are
\[
 \begin{aligned}
 &a_2(2L+a_2\eta d+2dv)/4,\\
 &(12Lv-a_2^3d+6dv^2)/24,\\
 &-a_2^2(3L+2a_2\eta d+3dv)/24,\\
 &-a_2(24La_2\eta+48Lv-3a_2^3d+8a_2^2\eta^2 d
                         +24a_2\eta dv+24dv^2)/192 .
 \end{aligned}
\]
These are obtained by expanding the two finite sums and
replacing \(Y^3\) by \(s+\eta\epsilon Y\).
All four vanish by regularity of \(Q\) in the free basis.
Since \(a_2d\ne0\), the first and third force
\(\eta=0\) and \(L=-dv\). The second and fourth now give
\(a_2^3=-6v^2\) and \(a_2^3=-8v^2\), respectively.
This contradicts \(a_2\ne0\).

All cases yield \(y^m\mid E_B\) or a contradiction. The
reference descent therefore excludes the nonbinary tuple.
