# Aligned powers in mixed sextic degree profiles

**Theorem SexticKellerAlignedLinearPowerFirstWeight.**
Let \((P,Q,T)\) be a characteristic-zero three-variable Keller map
with ordinary highest forms \(x^e,x^6,H\), where \(e\in\{4,5\}\),
\(H\) has degree five, and \(h=H(0,y,z)\ne0\). Assume that
\(H/x^5\) is a closed parameter in the projective plane function
field over the algebraic closure. One may subtract a constant
multiple of \(T\) from \(Q\) so that the degree-five part of \(Q\)
vanishes at \(x=0\). In this frame decrease the weight of \(x\)
from one, keeping the weights of \(y,z\) equal to one. There is a
first weight \(w\in(0,1)\) where another monomial ties \(x^e\) or
\(x^6\) in \(P\) or \(Q\), respectively. Then the following hold
over the algebraic closure.

- If \(e=5\), \(h\) is a fifth power of a linear form \(y_0\),
  and the highest forms at weights \((w,1,1)\) depend only on
  \(x,y_0\).
- If \(e=4\), \(w=1/2\) and the
  highest forms at that weight are \((W^2,W^3,h)\), where
  \(W=x^2+\ell(y,z)\) for a nonzero linear form \(\ell\).

The first case applies to the nonbinary
[reduced profile \((6,5,5)\)](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md).
The second applies to the linear-power alternative of the
irreducible-quintic reduction in profile \((6,5,4)\).

**Theorem SexticKeller655AlignedPowersExclusion.**
There is no characteristic-zero three-variable Keller map with
ordinary highest forms \(x^5,x^6,H\), where \(H\) has degree five,
\(H(0,y,z)\ne0\), and \(H/x^5\) is a closed projective parameter.
Consequently a Keller map of ordinary component degrees
\((6,5,5)\) whose quintic highest forms are coprime and squarefree
has all three highest forms in the same two linear source forms
over the algebraic closure.

**Theorem KellerPureBoundaryQuinticExactArea.**
Suppose two components of a three-variable Keller map have ordinary
highest forms \(x^e\) and \(H\), where \(e\ge1\), \(H\) has degree five,
and \(H(0,y,z)\) is a nonzero fifth power of a linear form. In a linear
source and scalar frame with \(H(0,y,z)=y^5\), the plane polynomial
\(K(y,z)=H(1,y,z)\) has one of the following forms.
\[
\begin{array}{ll}
\text{(i)} & K\in k[y];\\[2pt]
\text{(ii)} & K=\mu(y-\rho)^j z+b(y),\qquad
         \mu\ne0,\quad j\in\{0,2,3\},\quad b(y)=y^5+O(y^4);\\[2pt]
\text{(iii)} &
 K=a(y)(z+q(y))^2+b_0(z+q(y))+c_0,\qquad
 \deg a=1,\quad \deg q=2 .
\end{array}
\]
Here the coefficients are taken over the algebraic closure, and
the coefficient of \(y^5\) in (iii) is one. The closedness of
\(H/x^5\) excludes (i).

**Theorem SexticKeller654AlignedNonpureQuinticBoundary.**
In the case \(e=4\) of the first-weight theorem, write
\(W=x^2+\ell(y,z)\) for the quadratic parameter of the
first weighted pair. If \(h\) is not a fifth power, then
\(h\) is proportional to \(\ell^3 L^2\) for a linear
form \(L\) independent of \(\ell\).

## The first weight at which the pure powers change

Work over an algebraic closure \(k\). Write \(P_{e-1}\) and \(Q_5\)
for the indicated homogeneous parts. The determinant equation
one degree below the highest implies
\[
 \det J(x,H,Q_5-(6/e)x^{6-e}P_{e-1})=0.
\]
The rational kernel of \(\det J(x,H,\cdot)\) is \(k(x,H)\):
closedness of \(H/x^5\) means its function field is relatively
algebraically closed in the projective plane field. The extension
is regular in characteristic zero, so adjoining the independent
radial coordinate \(x\) preserves this relative closedness.
The pair \((x,H)\) is onto \(\mathbb A^2\), since \(h\) is a
nonzero degree-five form and \(H(a,y,z)\) has highest part \(h\)
for every \(a\in k\). A reduced fraction in \(k(x,H)\) with
nonconstant denominator has a denominator zero off its numerator
in the target plane, and a preimage there precludes a polynomial
pullback. Thus the polynomial kernel is \(k[x,H]\). Its degree-five
part is \(\langle x^5,H\rangle\), which gives the stated target
subtraction.

If both \(P,Q\) were polynomials in \(x\), the full determinant
would vanish. Every monomial \(x^i y^j z^l\) in \(P\) with
\(j+l>0\) ties \(x^e\) at weight \(w=(j+l)/(e-i)\in(0,1)\).
Use \(6\) in place of \(e\) for \(Q\). Take the largest such
weight. The highest forms \(P_*,Q_*\) have weights \(ew,6w\)
and are monic of degrees \(e,6\) in \(x\). The highest form of
\(T\) is \(h\) because \(w<1\). The top determinant vanishes:
its weight \((e+5)w+3\) is positive. Thus \(P_*,Q_*\) are
algebraically dependent over \(k(h)\).

Here we use the following polynomial parametrization observation.
Let \(L\) be an algebraic closure of \(k(y,z)\) and \(K_0\) the
algebraic closure of \(k(h)\) inside \(L\). The irreducible affine
curve with ring \(K_0[P_*,Q_*]\) has normalization an affine line
over \(K_0\). Indeed the coefficients of \(P_*,Q_*\) generate a
finite type \(K_0\)-subalgebra of \(L\) and can be specialized to
\(K_0\) at a closed point. Monicity in \(x\) preserves the degrees,
so the map from the coordinate ring of the curve to \(K_0[x]\)
remains injective: a nonzero prime kernel would make the image
zero-dimensional. The resulting morphism from the affine line
is finite since \(x\) is integral over the image of \(K_0[P_*]\).
It lifts to the normalization, since every integral element of
the curve function field maps into the integrally closed
\(K_0[x]\). The corresponding projective curve is rational by
Lüroth's theorem and has exactly one point omitted by the affine
normalization, since any omitted point can be covered only by
the point at infinity of the projective source line.

Let \(U\) be the parameter of this normalization, viewed in
\(L[x]\) by integrality in the original realization. Then
\(P_*,Q_*\in K_0[U]\) and \(\deg_x U=m\) divides \(e\) and \(6\).
The leading coefficient of \(U\) is in \(K_0^\times\), since
\(P_*\) is monic, and we can make \(U\) monic by rescaling.
Translate \(U\) by an element of \(K_0\) to remove the coefficient
of \(U^{e/m-1}\) in \(P_*\). Now \(U\) is the polynomial part
in \(x\) of the formal monic root \(P_*^{m/e}\), since the
terms after \(U^{e/m}\) decrease the \(x\)-degree by at least
\(2m\). Expanding that root at \(x=\infty\) shows that \(U\)
is in \(k[y,z][x]\) and is homogeneous of weight \(mw\).
The coefficients of \(P_*,Q_*\) as polynomials in \(U\) lie in
\(K_0\cap k[y,z]\), by monic division, and each is homogeneous
of the residual weight.

A homogeneous \(f\in K_0\cap k[y,z]\) of degree \(s\), with
\(0<s<5\), forces \(h\) to be a fifth power up to scalar.
Indeed \(df\wedge dh=0\) and the Euler derivation give
\(5h\,df=s f\,dh\), so \(f^5/h^s\) is a nonzero constant.
The multiplicities of all factors of \(h\) are divisible by
five, as \(\gcd(5,s)=1\). If \(h=y_0^5\), any such polynomial
coefficient, of any degree, belongs to \(k[y_0]\), since
\(k(y_0)\) is relatively algebraically closed in \(k(y,z)\).

If \(e=5\), then \(m=1\) and \(U=x\), as \(0<w<1\) leaves
no term of weight \(w\) in \(y,z\). If \(h\) is not a fifth
power, every nonleading coefficient of \(P_*,Q_*\) vanishes:
the residual weights are positive and less than five, except
possibly for a degree-five coefficient at \(x^0\) in \(Q_*\).
That coefficient vanishes by the target normalization.
This contradicts the choice of \(w\). For \(h=y_0^5\) the same
parametrization gives \(P_*,Q_*\in k[x,y_0]\).

For \(e=4\) with \(h\) not a fifth power, \(m=1\) again
gives a contradiction. For \(m=2\), the coefficient of \(x\)
in \(U\) must vanish by its weight \(w\).
If \(U=x^2\), the same residual-weight argument gives only pure
powers, which is impossible. Hence \(2w=1\) and
\(U=x^2+\ell(y,z)\) with \(\ell\ne0\) linear.
Now the weights of \(P_*,Q_*\) are two and three, so each
nonleading coefficient in \(U\) vanishes. This proves the
first-weight statement in that case. For \(e=4\) with \(h=y^5\),
we use the formal expansion below. If \(w\ne1/2\), then even
when \(m=2\) the parameter \(U\) is \(x^m\) by its weight, so
\(P_*,Q_*\in k[x,y]\).

## Formal expansion at the first weight

Use the normalized frame from the first-weight theorem, with
\(h=y^5\), and write \(w=p/q\) in lowest terms. By construction
\(q\le6\), and the only option with \(p=5\) is \(w=5/6\).
This weight could only come from a degree-five monomial of \(Q\)
independent of \(x\), which was removed. Thus \(1\le p\le4\).

Use ordinary scaling \(F_\epsilon=\epsilon^d
F(x/\epsilon,y/\epsilon,z/\epsilon)\) on a component of degree
\(d\). Its determinant for the tuple is a nonzero constant
times \(\epsilon^N\), where \(N=e+8\). The formal monic root
\(X=P_\epsilon^{1/e}=x+O(\epsilon)\) belongs to
\(k[x,x^{-1},y,z][[\epsilon]]\). In this ring,
\[
 Q_\epsilon \equiv \sum_{i=0}^{N-1}\epsilon^i D_i(X,T_\epsilon)
          \pmod{\epsilon^N},\qquad
 D_0(X,T_\epsilon)=X^6,\quad D_i\in k[X,X^{-1},T_\epsilon].
 \tag{1}
\]
In fact, after subtracting any sum of the displayed kind, its
bracket with \(P_\epsilon,T_\epsilon\) is unchanged since
\(X^e=P_\epsilon\). At the first unmatched order \(i<N\),
the residual coefficient is therefore in the rational kernel
of \(\det J(x,H,\cdot)\), hence in \(k(x,H)\). It is a
Laurent polynomial localized at \(x\), so multiplying by a
large enough power of \(x\) and using polynomial contraction
puts it in \(k[x,x^{-1},H]\). Taking this coefficient for
\(D_i(x,H)\) and lifting to \(X,T_\epsilon\) raises the
unmatched order. With \(X,T_\epsilon\) of ordinary scaling
degrees one and five, \(D_i\) is homogeneous of degree \(6-i\).
Its monomials thus have form \(c_{ij} X^{6-i-5j}T_\epsilon^j\)
with \(j\ge0\).
This is a two-first-integral version of formal root subtraction
for plane Jacobian pairs.
[Guccione–Guccione–Valqui, Proposition 1.7](https://arxiv.org/pdf/1406.0886v3)

Give \(x,\epsilon,y,z\) the weights \(a=1-w,1,0,0\), and
now take *minimum* weight. All weights of \(P_\epsilon,Q_\epsilon\)
are at least \(ea,6a\), respectively. Their minimum-weight parts
are the earlier highest \((w,1,1)\)-parts lifted to the rescaling.
In the cases just shown to have binary highest forms they are
\[
             x^e f(t),\qquad x^6 g(t),\qquad
             t=\epsilon^{q-p}y^p/x^q,
\]
where \(f(0)=g(0)=1\), \(\deg f\le\lfloor e/q\rfloor\) and
\(\deg g\le\lfloor6/q\rfloor\). At least one of \(f,g\) is
nonconstant. Also the minimum-weight part of \(T_\epsilon\) is
\(y^5\). The expansion \(X^v=x^v(1+(P_\epsilon/x^e-1))^{v/e}\)
for any integer \(v\) has weights at least \(va\); terms of weight
\(va\) are obtained from \(x^v f(t)^{v/e}\).

Every nonzero monomial \(\epsilon^i c_{ij}X^{6-i-5j}T_\epsilon^j\)
in (1) must have minimum weight \(i+a(6-i-5j)\ge6a\).
To see this, induct on \(i\) in (1). The coefficient of
\(\epsilon^i\), after the earlier terms are removed, has weight
at least \(6a\) including its \(\epsilon^i\) factor.
The distinct \(j\) in \(D_i\) have different weights at that
order, so no term of smaller weight can cancel. At equality
\(i p=5(q-p)j\). Thus \(i=n(q-p)\), \(j=np/5\),
and \(n\) is a multiple of five since \(p\le4\). The identity
in weight \(6a\) is
\[
 g(t)\equiv f(t)^{6/e}
       +\sum_{\substack{n\ge5,\;5\mid n\\n(q-p)<N}}
                  c_n\, t^n f(t)^{(6-nq)/e}
       \pmod{t^{\lceil N/(q-p)\rceil}}.
 \tag{2}
\]
This comparison only uses powers of \(\epsilon\) below \(N\).
All binomial expansions have constant term one in the base.

When \(q=6\), \(f=1\) and (2) at order \(t\) gives \(g=1\).
When \(q=4,5\), both \(f,g\) have degree at most one.
Vanishing of the coefficient at \(t^2\) on the right gives
\(\binom{6/e}{2}([t]f)^2=0\). Then \(f=g=1\).
When \(q=3\), \(f\) has degree at most one and \(g\) has degree
at most two. The coefficient at \(t^3\) gives the same
conclusion since \(\binom{6/e}{3}\ne0\).

For \(e=5\) and \(q=2\), \(p=1\), write \(f=1+A t+B t^2\), with
\(\deg g\le3\). The coefficient at \(t^4\) on the right of
(2) is
\[
                3(A^2-5B)(3A^2-5B)/625 .
\]
The coefficient at \(t^5\) determines
\[
          c_5=\frac{6A(21A^4-150A^2B+250B^2)}{15625}.
\]
Substitute in the coefficient at \(t^6\) to get
\[
             -\frac{21A^6-90A^4B-150A^2B^2+500B^3}{15625}.
\]
These three coefficients follow by expanding
\(f^{6/5}+c_5t^5 f^{-4/5}\) through \(t^6\).
If \(A\ne0\), the \(t^4\) equation gives \(B=A^2/5\)
or \(3A^2/5\). Then the \(t^6\) coefficient is, respectively,
\(-A^6/15625\) or \(-21A^6/15625\), contrary to
\(\deg g\le3\). If \(A=0\), the \(t^4\) equation gives \(B=0\)
and (2) gives \(g=1\).

For \(e=5\) the two edge polynomials are constant in every case,
contradicting the definition of the first weight. This excludes the aligned
linear-power tuple. The binary conclusion follows from
[SexticKeller655ReducedLowerPencil](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md).

For \(e=4\) and \(h=y^5\), this leaves only \(w=1/2\),
the case where we did not require the highest forms to be binary.
In this case write the minimum-weight part of \(P_\epsilon\) as
\[
                    x^4(1+L\tau+M\tau^2),\qquad
                    \tau=\epsilon/x^2,
\]
where \(L,M\) are homogeneous in \(y,z\) of degrees one and two.
The minimum-weight part of \(Q_\epsilon\) has degree at most
three in \(\tau\) after dividing by \(x^6\). The same
order-by-order weight argument for (1) applies before expressing
the edge polynomials in \(t\): its correction terms at minimum
weight satisfy \(i=5j\). Hence through \(\epsilon^4\) that
part of \(Q_\epsilon\) is given by the binomial expansion
of \(x^6(1+L\tau+M\tau^2)^{3/2}\).
The \(\tau^4\) coefficient is \(3(L^2-4M)^2/128\), so
\(M=L^2/4\). The terms through \(\tau^3\) now give the
cube with the same quadratic parameter. Thus the highest pair
is \(((x^2+L/2)^2,(x^2+L/2)^3)\), where \(L\ne0\) by the
choice of the first weight.

## The exact differential at the pure quintic boundary

The [leading-pair rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
for \(x^e,H\) gives a rational primitive of the relative area
differential of \(K\) on its geometric generic level. Indeed,
fix a nonzero generic value of \(x\) and rescale \(y,z\) by
that value. The derivative of \(x^e\) and the change of scale
only change the differential by a nonzero scalar on each level.
Write \(S\) for the indeterminate level of \(K\). Its affine locus
is smooth by generic smoothness in characteristic zero, and the
relative area differential is nonzero and regular there.

There is a single projective point at infinity, in the chart
\(z=1\) with coordinates \(s,u\) corresponding to \(1/z,y/z\)
in the original plane. Put
\[
 F(s,u)=H(s,u,1)-S s^5,\qquad F(0,u)=u^5.
\]
The differential at infinity is, up to a nonzero sign,
\(s^2\,ds/F_u=-s^2\,du/F_s\). In fact the affine area has
denominator \(s^3\) and the level equation has denominator \(s^5\).
A simple pole is incompatible with a rational primitive. If all
branches at infinity have regular differential, it is regular
on the complete normalization and cannot have a rational
primitive either.

Let \(m\) denote the multiplicity of \(F\) at the origin.
If \(m<5\), the tangent cone contains \(s=0\), since \(F(0,u)=u^5\).
For \(m\le3\) use the blow-up chart \(u=s v\) to study the
other tangent directions. Put \(G=F(s,s v)/s^m\). The differential
is \(s^{3-m}ds/G_v\). A simple root of \(G(0,v)\) gives a
regular differential. A double root can occur only for \(m=3\),
as the tangent cone already has the factor \(s\). At a double
root \(v_0\), if \(G_s(0,v_0)\ne0\), \(v-v_0\) parametrizes a
smooth branch on which \(s\) has order two. Both \(ds\) and \(G_v\)
have order one. If \(G_s(0,v_0)=0\), the quadratic initial
polynomial in \(s,v-v_0\) has two distinct linear factors
transverse to \(s=0\), because the coefficient of \(s^2\)
contains the generic \(-S\) and that of \((v-v_0)^2\) is
nonzero. On both branches \(ds/G_v\) has a simple pole.
Thus every branch off the tangent \(s=0\) is regular or
gives a forbidden simple pole.

If \(s\) has multiplicity one in the tangent cone at \(m\le3\),
the branch in that direction is smooth, with \(u\) a parameter.
The other tangent directions have total contact \(m-1\) with
\(s=0\), so this branch has \(\operatorname{ord}(s)=6-m\).
The other branches contribute order \(m-1\) to \(F_s\) on
this branch: their initial product of tangent forms, evaluated
at \(s=0\), is a nonzero multiple of \(u^{m-1}\), and the
smooth branch factor has unit \(s\)-derivative. Therefore
\(s^2du/F_s\) is regular, giving an exclusion.

The remaining case with \(m<3\) is the double tangent \(s^2\)
at \(m=2\). Weierstrass preparation in \(s\) gives a unit
times \(s^2+A(u)s+B(u)\), with \(\operatorname{ord} B=5\)
and \(\operatorname{ord} A\ge2\).
If the discriminant has order four there are two branches,
with \(s\)-orders two and three in the parameter \(u\) and
\(F_s\)-order two on both. If the discriminant has order five
there is one branch with \(u\)-order two, \(s\)-order five,
and \(F_s\)-order five. In each case \(s^2du/F_s\) is
regular. These are the only discriminant orders because
\(A^2\) and \(4B\) have unequal orders. This finishes the
exclusion of \(m<3\), and at \(m=3\) the tangent \(s=0\)
must occur at least twice.

Consequently \(K\) is at most quadratic in \(z\), with its
quadratic coefficient of degree at most one in \(y\):
if \(m=3\) the degree-three tangent form is a multiple of
\(s^2\), and otherwise the quadratic coefficient vanishes.
The coefficient of \(z\) has degree at most three in \(y\)
because the ordinary degree-five part of \(K\) is \(y^5\).

If \(K=A(y)z+b(y)\) with \(A\ne0\), the level is rational
with parameter \(y\) and the relative differential is \(dy/A(y)\).
Any simple root of \(A\) gives a simple pole.
A polynomial of degree at most three without a simple root
is constant or a single power of degree two or three.
This gives (ii).

If \(K=a(y)z^2+b(y)z+c(y)\) with \(a\ne0\), put
\[
 u_1=2a(y)z+b(y),\qquad
 D(y)=b(y)^2-4a(y)c(y)+4a(y)S.
\]
The level has equation \(u_1^2=D(y)\) in its function field
and differential \(dy/u_1\).
The polynomial \(D\) is squarefree generically: away from
\(a=0\) a repeated root would give a critical point of the
rational function \((b^2-4ac)/a\) at the generic value
\(-4S\), which is impossible. If that rational function is
constant, \(D\) is just a scalar multiple of \(a\).
At the root of a nonconstant \(a\), a persistent root of \(D\)
has nonzero derivative because the coefficient of \(S\) in
that derivative is nonzero. For degree at least three,
\(dy/u_1\) is regular at every point of the complete
normalization. At infinity this follows by taking \(1/y\)
for even degree and \(1/\sqrt y\) for odd degree. Degree
two gives simple poles at infinity. Thus \(\deg D\le1\).
If \(a\) were constant, \(\deg c=5\) and \(\deg b\le3\)
would give \(\deg D\ge5\). So \(a\) is linear nonconstant
and \(b^2-4ac\) is affine linear. Write \(b=2a q+b_0\).
After completing the square the new constant part
\(c-aq^2-b_0q\) is constant in \(y\), since its product
with \(4a\) is \(b_0^2-(b^2-4ac)\). Finally \(\deg q=2\)
by the fifth-degree term of \(K\). This proves (iii).

If \(K\) is independent of \(z\), then \(H/x^5=K(y/x)\)
has degree five as a rational function of \(y/x\), so it
cannot be closed. The displayed necessary forms concern
leading parts of Keller maps.

## The nonpure quintic in the first weighted \((6,5,4)\) tuple

Choose \(y=\ell\), so that \(W=x^2+y\).
The weighted leading-pair slice for \(W^2,h\) forces exactness
of the relative volume differential on their geometric generic
fibre. Fix \(W=s_0\), \(h=b\), where \(s_0,b\) are
algebraically independent. The differential on
\(x^2=s_0-y,\ h(y,z)=b\) is \(1/(2x)\) times the relative
area of the plane level \(h=b\). This is obtained by using
\(W,y,z\) as coordinates. Fixing \(W^2\) instead of \(W\)
only changes the form by a nonzero scalar.

Put \(u=z/y\) and \(f(u)=h(1,u)\), so \(y^5 f(u)=b\).
The differential, up to a nonzero scalar, is
\[
                         \omega=\frac{y^2\,du}{x}.
\]
Indeed \(dy\wedge dz=y\,dy\wedge du\) and the derivative of
\(y^5 f(u)\) in \(y\) on the fixed level is \(5b/y\).
At a finite root of \(f\) of multiplicity \(m\), write \(r\)
for the order of \(u-u_0\) in a uniformizer of the curve.
Then \(y\) has pole order \(rm/5\), \(x\) has half that
pole order, and \(\omega\) has order \(r(1-3m/10)-1\),
which is nonnegative if \(m\le3\).
At \(u=\infty\), let \(m_y\) be the multiplicity of \(y\)
in \(h\). As \(h\) is not pure, \(\deg f=5-m_y>0\).
For \(r\) the pole order of \(u\), \(y\) has zero order
\(r(5-m_y)/5\), \(x\) is a unit, and the form has order
\(r(1-2m_y/5)-1\). This is nonnegative if \(m_y\le2\).
At all other points the level lies in the affine chart, where
the slice differential is regular by generic smoothness.
Thus some linear factor other than \(y\) has multiplicity
four, or \(m_y\ge3\). Otherwise the nonzero differential
would be regular on the complete curve and could not be exact.

If \(m_y=4\), a linear shift in \(z\) and a scalar change
give \(h=y^4 z\). The curve is rational with parameter \(x\)
and its differential is proportional to \(dx/(s_0-x^2)^4\).
At \(x=x_0\), \(x_0^2=s_0\), the nonzero residue comes
from the coefficient of \((x-x_0)^3\) in
\((2x_0+(x-x_0))^{-4}\). This excludes \(m_y=4\).

If \(m_y=3\) and the two residual linear factors are distinct,
complete the square in \(z\) by a linear shift and rescale
to write \(h=y^3(z^2+c y^2)\), \(c\ne0\).
Set \(V=y^2 z\). In the function field after substituting
\(y=s_0-x^2\) one has \(V^2=b y-c y^6\) and
\(\omega\) is proportional to \(dx/(y V)\).
The polynomial on the right in \(x\) has degree twelve and is
squarefree: its six distinct \(y\)-roots lift to distinct
\(x\)-roots for generic \(s_0\). The differential is regular
except for double poles at the two points above \(y=0\).
At these points \(x-x_0\) has order two and \(V\) order one.
At the other finite roots of \(V\), the order of \(dx\)
compensates that of \(V\). At the two points at infinity
\(x,y,V\) have pole orders one, two and six respectively.

If an exact primitive exists, taking its anti-invariant part
under \(V\mapsto -V\) gives a primitive \(V R(x)\) for a
rational \(R\). It has at most simple poles at the two
points above \(y=0\) and is regular elsewhere. Thus \(R\)
has no finite poles except possibly simple poles at
\(x=\pm\sqrt{s_0}\). A pole at another branch point would
give a pole of order at least one in \(VR\), since a pole
of \(R\) has even order there and \(V\) has order one.
Hence the denominator of \(R\) divides \(s_0-x^2\).
If \(R\ne0\), \(VR\) then has a pole at infinity, which
is a contradiction.

If a factor other than \(y\) has multiplicity four,
take it as \(z\) and write \(h=z^4(y+c z)\) up to scalar.
Set \(V=z^2 x\). The generic level has
\(V^2=z^4(s_0+c z)-b\) and differential proportional to
\(dz/(z^2V)\). The polynomial on the right is squarefree
by genericity of \(b\), and has degree four or five.
The only poles of the differential are double poles at the
two points over \(z=0\), where \(V\) is a unit. At the
finite branch points \(dz\) cancels \(V\), and the form
is regular at infinity. An anti-invariant primitive would
be \(VR(z)\) with denominator of \(R\) dividing \(z\),
by the same order comparison, now at unramified points
over \(z=0\). Multiplication by \(V\) gives a pole at
infinity for any nonzero \(R\) of this kind, again a
contradiction.

The only remaining nonpure form is \(h\) proportional to
\(y^3 L^2\) with \(L\) independent of \(y=\ell\). For
example, for \(h=y^3 z^2\) one can take \(y=r^2\) and
\(z=\sqrt b\,r^{-3}\) on the plane level. The resulting
curve \(x^2+r^2=s_0\) has
\(d(x/r)=-s_0\,dr/(x r^2)\), so this relative differential
is exact.
