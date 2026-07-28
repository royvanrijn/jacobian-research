# Hopf-lift classification and minimal \(t^2\)-linear rigidity

## 1. Scope

This note isolates the mechanism behind the adjacent Taylor coefficients in
the two-pair Image-Mathieu counterexample.  It first classifies a parametric
one-profile class, then proves uniqueness in the complete minimum-degree
\(t^2\)-linear numerator class.  Neither result is a search in the full space
\(V_d=\operatorname{End}(\operatorname{Sym}^d)\).

On the Hopf sphere use

\[
 t^2+2xy=1,
 \tag{1.1}
\]

where the phase of \(x\) is uniform and \(t\) is uniform on \([-1,1]\).
Fix

\[
 X=1+ax,\qquad a\ne0.
 \tag{1.2}
\]

The **one-profile Hopf class** consists of Laurent presentations

\[
 p_{r,q,R}
 =c\,x^{-r}X^qR(t^2X^2),
 \qquad c\ne0,\quad r\ge1,\quad q\ge0,
 \tag{1.3}
\]

where \(R\in k[z]\) has a zero of order \(s\) at \(z=1\), with \(s\ge r\).
Here \(k\) has characteristic zero; for the nonvanishing statements below
we use the real integral and assume its displayed values are nonzero.

The restrictions in (1.3) have geometric content.  There is one affine
phase profile \(X\), the \(t\)-dependence sees only the Hopf-invariant
combination \(t^2X^2\), and \(x^{-r}\) records the phase winding.  Arbitrary
linear combinations \(\sum_jB_j(x)t^{2j}\) are not claimed to reduce to
this class.  If \(R(z)=\sum_j\rho_jz^j\), then (1.3) is the explicit
subclass of the general ansatz
\[
 \frac{A(x)}{x^r}\sum_jB_j(x)t^{2j}
 \quad\text{with}\quad
 A(x)=cX^q,\qquad B_j(x)=\rho_jX^{2j}.
\tag{1.4}
\]

## 2. Exact polynomiality criterion

Write

\[
 R(z)=(1-z)^sS(z).
 \tag{2.1}
\]

In the sphere coordinate ring,

\[
\begin{aligned}
 1-t^2X^2
 &=1-t^2(1+ax)^2\\
 &=x\left(2y-t^2(2a+a^2x)\right).
\end{aligned}
\tag{2.2}
\]

Thus \(R(t^2X^2)\) is divisible by \(x^s\) modulo (1.1).  The condition
\(s\ge r\) is therefore a transparent sufficient condition for (1.3) to
represent a polynomial on the sphere.  It is also necessary.

> **Proposition 2.1 (polynomiality equals endpoint contact).**  Suppose
> \(a,c\ne0\), \(q\ge0\), and \(R\ne0\).  Then the Laurent presentation
> \(p_{r,q,R}\) in (1.3), temporarily without assuming \(s\ge r\),
> represents an element of the sphere coordinate
> ring if and only if
> \[
> \operatorname{ord}_{z=1}R(z)\ge r.
> \tag{2.3}
> \]

Indeed, eliminate \(t^2\) by (1.1) and work in \(k[x,y]\).  If
\[
 z=(1-2xy)(1+ax)^2,
\]
then
\[
 1-z=xD(x,y),\qquad D(0,y)=2(y-a)\ne0.
\tag{2.4}
\]
Also \(X^q\equiv1\pmod x\) and \(S(z)\equiv S(1)\ne0\pmod x\).  Hence the
numerator in (1.3) has exact \(x\)-adic order \(s\).  Division by \(x^r\)
is polynomial exactly when \(s\ge r\).

Thus polynomiality and the Taylor cancellation are controlled by the same
integer: contact order with the Hopf endpoint \(t^2X^2=1\).

## 3. The classification theorem

Define the phase-integrated polynomial

\[
 H_{m,q,R}(X)
 =X^{qm}\int_0^1R(v^2X^2)^m\,dv
 \tag{3.1}
\]

and put

\[
 C_m=\int_0^1R(v^2)^m\,dv.
 \tag{3.2}
\]

The letter \(v\) in (3.1) is the uniform Hopf-height variable.  It is
renamed \(w=vX\) in the proof.

> **Theorem 3.1 (one-profile classification).**  Let \(R\) have order
> \(s\ge r\) at \(1\), and suppose \(C_m\ne0\) for every \(m\ge1\).  Then
> the following are equivalent:
>
> 1. for every \(m\ge1\),
>    \[
>    [u^{rm}]H_{m,q,R}(1+u)=0,\qquad
>    [u^{rm-1}]H_{m,q,R}(1+u)\ne0;
>    \tag{3.3}
>    \]
> 2. \(q=r\).
>
> When these conditions hold, the adjacent coefficient is exactly
> \[
> [u^{rm-1}]H_{m,r,R}(1+u)=C_m.
> \tag{3.4}
> \]

For the primitive winding \(r=1\), (3.3) is precisely

\[
 [u^m]H_m(1+u)=0,\qquad
 [u^{m-1}]H_m(1+u)=C_m\ne0.
 \tag{3.5}
\]

### Proof

Changing variables \(w=vX\) in (3.1) gives

\[
 H_{m,q,R}(X)
 =X^{qm-1}J_m(X),\qquad
 J_m(X)=\int_0^X R(w^2)^m\,dw.
 \tag{3.6}
\]

Because \(R(w^2)\) has a zero of order \(s\) at \(w=1\),
\(J_m'(w)=R(w^2)^m\) has a zero of order \(sm\) there.  Hence

\[
 J_m(1+u)=C_m+O(u^{sm+1})
          =C_m+O(u^{rm+1}).
 \tag{3.7}
\]

Through degree \(rm\), the Taylor polynomial of (3.6) is therefore the
Taylor polynomial of

\[
 C_m(1+u)^{qm-1}.
\tag{3.8}
\]

If \(q=0\), its coefficient at degree \(rm\) is
\((-1)^{rm}C_m\ne0\), so (3.3) fails.  We may therefore assume \(q\ge1\).
The adjacent coefficient in (3.3) can be nonzero only if
\(qm-1\ge rm-1\), so \(q\ge r\).  The coefficient at degree \(rm\) can
vanish only if \(qm-1<rm\), so \(q\le r\).  Thus \(q=r\).  Conversely,
\((1+u)^{rm-1}\) has zero coefficient at degree \(rm\) and coefficient one
at degree \(rm-1\), which proves (3.3)--(3.4).  \(\square\)

The proof also shows why the two coefficients are genuinely adjacent:
polynomiality supplies at least \(r\) endpoint zeros, and the two coefficient
conditions pin the profile exponent to the same integer \(r\).

The same proof gives more than the adjacent pair.

> **Corollary 3.2 (full lower-jet ladder).**  Under the hypotheses of
> Theorem 3.1, with \(q=r\), one has for
> \(1\le\ell\le rm\)
> \[
> \boxed{
> [u^{rm-\ell}]H_{m,r,R}(1+u)
> =\binom{rm-1}{\ell-1}C_m.}
> \tag{3.9}
> \]
> The coefficient at degree \(rm\) is zero.

Indeed, all these degrees are at most \(rm\), where (3.7) replaces \(J_m\)
by \(C_m\), and
\[
 [u^{rm-\ell}](1+u)^{rm-1}
 =\binom{rm-1}{\ell-1}.
\tag{3.10}
\]
Thus the adjacent phenomenon is the first step of a complete binomial jet
ladder.

## 4. Phase moments and a parametric family

Before imposing (1.1), phase averaging is constant-term extraction in \(x\).
For \(q=r\) and \(1\le\ell\le rm\),

\[
\begin{aligned}
 \int_{\mathrm{phase}}p_{r,r,R}^m
 &=c^ma^{rm}[u^{rm}]H_{m,r,R}(1+u)=0,\\
 \int_{\mathrm{phase}}x^\ell p_{r,r,R}^m
 &=c^ma^{rm-\ell}[u^{rm-\ell}]H_{m,r,R}(1+u)\\
 &=c^ma^{rm-\ell}\binom{rm-1}{\ell-1}C_m.
\end{aligned}
\tag{4.1}
\]

Thus every \(R\) satisfying the hypotheses gives an all-order parametric
family, and every allowed positive phase multiplier detects it.  A simple
real sufficient condition for \(C_m\ne0\) is

\[
 R(z)=(1-z)^sS(z),\qquad S(z)>0\quad(0\le z\le1).
\tag{4.2}
\]

In particular \(R(z)=(1-z)^s\) gives

\[
 C_m=\int_0^1(1-v^2)^{sm}\,dv
 =\frac{2^{sm}(sm)!}{(2sm+1)!!}\ne0.
\tag{4.3}
\]

This family is larger than taking powers of the primitive example: an
arbitrary factor \(S(t^2X^2)\) satisfying (4.2) may be inserted without
destroying the adjacent jet.

## 5. Degree-minimal uniqueness and Long's polynomial

Now take primitive winding \(r=1\).  The theorem forces \(q=1\), while
polynomiality forces \(R(1)=0\).  Among nonzero \(R\), the minimum possible
degree is one, and then necessarily

\[
 R(z)=\lambda(1-z).
\tag{5.1}
\]

Absorbing \(\lambda\) into \(c\) gives the unique degree-minimal
one-profile lift

\[
 p_{a,c}
 =c\,\frac{1+ax}{x}
   \left(1-t^2(1+ax)^2\right).
\tag{5.2}
\]

Using (2.2), this is the five-term sphere polynomial

\[
\boxed{
 p_{a,c}
 =c(1+ax)\left(2y-(2a+a^2x)t^2\right).
 }
\tag{5.3}
\]

The torus change \(x'=ax,\ y'=a^{-1}y\) preserves (1.1), so all
\(a\ne0\) give one orbit, up to the overall scalar \(c\).  At
\(a=1,c=\tfrac12\),

\[
 p=y+xy-t^2-\frac32xt^2-\frac12x^2t^2,
\tag{5.4}
\]

which is Long's five-term polynomial in Hopf coordinates.  Its
phase-integrated function is

\[
 H_m(X)=X^m\int_0^1(1-v^2X^2)^m\,dv,
\tag{5.5}
\]

and (3.5) recovers the adjacent-jet argument used by the two-pair lift.

> **Corollary 5.1 (minimal uniqueness).**  Up to the sphere torus
> \(x\mapsto ax,\ y\mapsto a^{-1}y\) and nonzero scaling, Long's polynomial
> is the unique member of the one-profile Hopf class with primitive winding
> and linear \(R\).  Equivalently, it is the unique member of minimum
> \(R\)-degree in that class having the adjacent-jet phenomenon.

This is a uniqueness theorem for the Long/Hopf mechanism inside (1.3).
It is not global uniqueness among all low-degree sphere polynomials.  In
particular, it makes no claim about general \(V_4\), about ansatzes with
several independent phase profiles, or about cancellations not explained
by endpoint multiplicity.

## 6. Minimal-degree \(t^2\)-linear classification

The factorization in the preliminary ansatz
\[
 \frac{A(x)}x\left(B_0(x)+B_1(x)t^2\right)
\]
is redundant: only the products
\[
 C(x)=A(x)B_0(x),\qquad D(x)=A(x)B_1(x)
\tag{6.1}
\]
enter the sphere polynomial and its moments.  The invariant low-degree
question should therefore bound \(C,D\), not a chosen factorization.

Consider the minimum numerator degrees containing Long's polynomial:

\[
 p=\frac{C(x)+D(x)t^2}{x},
 \qquad \deg C\le1,\quad\deg D\le3.
\tag{6.2}
\]

Put
\[
 K_m(x)=\int_0^1(C(x)+v^2D(x))^m\,dv.
\tag{6.3}
\]

> **Theorem 6.1 (three-jet rigidity in the minimal linear class).**  For
> (6.2), the following are equivalent, up to multiplying \(p\) by a nonzero
> scalar:
>
> 1. \(p\) is polynomial on the sphere and, for every \(m\ge1\),
>    \[
>    [x^m]K_m(x)=0,\qquad [x^{m-1}]K_m(x)\ne0;
>    \tag{6.4}
>    \]
> 2. \(p\) is polynomial,
>    \[
>    [x]K_1=[x^2]K_2=[x^3]K_3=0,
>    \qquad [x]K_2\ne0;
>    \tag{6.5}
>    \]
> 3. for some \(a\ne0\),
>    \[
>    C(x)=1+ax,\qquad D(x)=-(1+ax)^3.
>    \tag{6.6}
>    \]
>
> Consequently Long's five-term polynomial is the unique polynomial with
> the all-order adjacent-jet phenomenon in the full \(t^2\)-linear class
> (6.2), up to the sphere torus and nonzero scaling.

### Proof

Eliminating \(t^2\) gives
\[
 C(x)+D(x)(1-2xy).
\tag{6.7}
\]
Thus (6.2) is polynomial exactly when \(C(0)+D(0)=0\).  Moreover,
the adjacent condition \([x]K_2\ne0\) forces \(C(0)\ne0\), since otherwise
both \(C,D\) are divisible by \(x\) and \(K_2\) is divisible by \(x^2\).
After scaling, write
\[
 C=1+ax,\qquad
 D=-1+b_1x+b_2x^2+b_3x^3.
\tag{6.8}
\]

Direct integration gives
\[
\begin{aligned}
 [x]K_1
 &=\frac{3a+b_1}{3},\\
 [x^2]K_2
 &=\frac{
 15a^2+10ab_1+3b_1^2+4b_2}{15},\\
 [x^3]K_3
 &=\frac1{35}\bigl(
 35a^3+35a^2b_1+21ab_1^2+28ab_2\\
 &\hspace{31mm}
 {}+5b_1^3+12b_1b_2+8b_3\bigr).
\end{aligned}
\tag{6.9}
\]
They form a triangular system.  Successive vanishing yields
\[
 b_1=-3a,\qquad
 [x^2]K_2=\frac4{15}(b_2+3a^2),\qquad
 [x^3]K_3=\frac8{35}(b_3+a^3),
\tag{6.10}
\]
so
\[
 b_1=-3a,\qquad b_2=-3a^2,\qquad b_3=-a^3.
\tag{6.11}
\]
Therefore \(D=-(1+ax)^3=-C^3\).  Under the first pure equation,
\[
 [x]K_2=\frac{8a}{15},
\tag{6.12}
\]
so the adjacent condition in (6.5) is exactly \(a\ne0\).

Now (6.6) is the primitive profile (5.2), and Theorem 3.1 proves all the
identities in (6.4).  Finally \(x\mapsto ax,\ y\mapsto a^{-1}y\) removes
the remaining nonzero parameter.  \(\square\)

Theorem 6.1 closes the first multi-coefficient island beyond the composite
ansatz (1.3).  Its finite-jet certificate is conceptual and triangular; it
does not enumerate supports or points of \(V_4\).

The triangularity is uniform in the degree.

> **Proposition 6.2 (universal triangular reconstruction).**  Normalize
> \[
> C(x)=1+\sum_{j\ge1}a_jx^j,\qquad
> D(x)=-1+\sum_{j\ge1}b_jx^j.
> \tag{6.13}
> \]
> For every \(m\ge1\), the pure jet \([x^m]K_m\) depends only on
> \(a_1,\ldots,a_m,b_1,\ldots,b_m\), and its coefficient of the new
> variable \(b_m\) is
> \[
> \tau_m
> =m\int_0^1v^2(1-v^2)^{m-1}\,dv
> =\frac{2^{m-1}m!}{(2m+1)!!}\ne0.
> \tag{6.14}
> \]
> Consequently the equations
> \([x^m]K_m=0\), \(1\le m\le N\), solve
> \(b_1,\ldots,b_N\) uniquely and successively in terms of
> \(a_1,\ldots,a_N\).

Indeed, terms of degree greater than \(m\) cannot contribute to
\([x^m]\).  Perturbing \(D\) by \(\epsilon x^m\), the coefficient linear in
\(\epsilon x^m\) is
\[
 m\int_0^1v^2(C(0)+v^2D(0))^{m-1}\,dv
 =m\int_0^1v^2(1-v^2)^{m-1}\,dv,
\tag{6.15}
\]
which evaluates to (6.14).  No term containing \(b_m^2\) can have degree
\(m\).  This proves the claim.

For the rectangle \((\deg C,\deg D)\le(d,d+2)\), the first \(d+2\) pure
jets therefore reconstruct all of \(D\).  After normalizing \(a_1=1\), the
entire remaining problem is an elimination problem in the \(d-1\) weighted
parameters \(a_2,\ldots,a_d\).  Theorems 7.1--7.4 identify those residual
ideals for \(d=2,3,4,5\).

## 7. Continuations, gaps, and obstructions

### 7.1 One further rectangle

The next numerator rectangle is
\[
 \deg C\le2,\qquad \deg D\le4.
\tag{7.1}
\]

It too is rigid.

> **Theorem 7.1 (six-jet rigidity in the \((2,4)\) rectangle).**  Let
> \(p=x^{-1}(C(x)+D(x)t^2)\) satisfy (7.1).  Up to nonzero scaling, the
> following are equivalent:
>
> 1. \(p\) is polynomial and has the all-order adjacent-jet phenomenon;
> 2. \(p\) is polynomial,
>    \[
>    [x^m]K_m=0\quad(1\le m\le6),\qquad [x]K_2\ne0;
>    \tag{7.2}
>    \]
> 3. \(C=1+ax,\ D=-(1+ax)^3\) for some \(a\ne0\).

### Proof

As in Theorem 6.1, polynomiality and the adjacent condition allow the
normalization
\[
\begin{aligned}
 C&=1+ax+ex^2,\\
 D&=-1+b_1x+b_2x^2+b_3x^3+b_4x^4,
\end{aligned}
\tag{7.3}
\]
with \(a\ne0\), since again \([x]K_2=8a/15\) after the first pure equation.
The first four pure jets successively give
\[
\begin{aligned}
 b_1&=-3a,\\
 b_2&=-3a^2-5e,\\
 b_3&=-a(a^2+12e),\\
 b_4&=-3e(3a^2+4e).
\end{aligned}
\tag{7.4}
\]
After these substitutions, the next two pure jets are
\[
\begin{aligned}
 [x^5]K_5
 &=\frac{256}{693}ae(a^2+12e),\\
 [x^6]K_6
 &=\frac{512}{3003}e
   (4a^4+63a^2e+20e^2).
\end{aligned}
\tag{7.5}
\]
Because \(a\ne0\), the fifth jet leaves \(e=0\) or
\(e=-a^2/12\).  On the second branch, the sixth jet is
\[
 [x^6]K_6=\frac{1280}{81081}a^6\ne0.
\tag{7.6}
\]
Hence \(e=0\), and (7.4) reduces to \(D=-(1+ax)^3\).  Theorem 3.1 again
supplies all orders.  \(\square\)

The extra quadratic direction therefore creates a finite-cutoff impostor:
it survives five pure jets but not six.  This illustrates why a bounded
calculation alone is insufficient unless the last surviving branches are
resolved exactly and then connected to an all-order identity.

The next rectangle is also rigid:
\[
 \deg C\le3,\qquad \deg D\le5.
\tag{7.7}
\]

> **Theorem 7.2 (eight-jet rigidity in the \((3,5)\) rectangle).**  Let
> \(p=x^{-1}(C(x)+D(x)t^2)\) satisfy (7.7).  Polynomiality, the pure
> identities
> \[
> [x^m]K_m=0\quad(1\le m\le8),
> \tag{7.8}
> \]
> and \([x]K_2\ne0\) force, up to scaling and the sphere torus,
> \[
> C=1+x,\qquad D=-(1+x)^3.
> \tag{7.9}
> \]
> Consequently these finite conditions are equivalent to the all-order
> adjacent-jet phenomenon throughout (7.7).

### Proof

The first pure equation and the adjacent condition again permit scaling
and the sphere torus to normalize
\[
\begin{aligned}
 C&=1+x+Ex^2+Fx^3,\\
 D&=-1+B_1x+B_2x^2+B_3x^3+B_4x^4+B_5x^5.
\end{aligned}
\tag{7.10}
\]
The first five pure jets solve triangularly:
\[
\begin{aligned}
 B_1&=-3,\\
 B_2&=-3-5E,\\
 B_3&=-1-12E-7F,\\
 B_4&=-3(3E+6F+4E^2),\\
 B_5&=-2E-15F-24E^2-36EF.
\end{aligned}
\tag{7.11}
\]
After removing nonzero rational factors, jets six through eight are
\[
\begin{aligned}
 P_6={}&4F+15E^2+78EF+20E^3+27F^2,\\
 P_7={}&F+4E^2+24EF+8E^3+12F^2+8E^2F,\\
 P_8={}&12F+51E^2+350EF+150E^3+231F^2\\
       &\quad+324E^2F+24E^4+90EF^2.
\end{aligned}
\tag{7.12}
\]
Eliminating \(F\) from \((P_6,P_7)\) and \((P_6,P_8)\) gives
\[
\begin{aligned}
 \operatorname{Res}_F(P_6,P_7)&=3E^2Q_5(E),\\
 \operatorname{Res}_F(P_6,P_8)&=192E^2Q_6(E),
\end{aligned}
\tag{7.13}
\]
where
\[
\begin{aligned}
 Q_5={}&11520E^5-16896E^4+4928E^3-1632E^2-156E+7,\\
 Q_6={}&6912E^6+64512E^5-169440E^4+32592E^3\\
       &\quad-20117E^2-1734E+75.
\end{aligned}
\tag{7.14}
\]
The Euclidean algorithm over \(\mathbb Q[E]\) gives
\[
 \gcd(Q_5,Q_6)=1.
\tag{7.15}
\]
Thus every common zero of \(P_6,P_7,P_8\) has \(E=0\).  Then
\[
 P_6=F(4+27F),\qquad P_7=F(1+12F),
\tag{7.16}
\]
whose only common zero is \(F=0\).  Equations (7.10)--(7.11) now give
\(C=1+x\) and \(D=-(1+x)^3\), and Theorem 3.1 supplies all orders.
\(\square\)

The proof is set-theoretic.  Exact Gröbner reduction strengthens
\(\langle P_6,P_7,P_8\rangle\) to
\(\langle F,E^2\rangle\), so the finite-jet scheme has a nonreduced
transverse direction even though its only point is Long's.  The theorem
does not need that stronger scheme-theoretic statement.

The next rectangle is again rigid, but its shortest current certificate is
an exact Gröbner reduction rather than a hand elimination.

> **Theorem 7.3 (ten-jet rigidity in the \((4,6)\) rectangle).**  Let
> \[
> \deg C\le4,\qquad \deg D\le6.
> \tag{7.17}
> \]
> Polynomiality, the pure identities through order ten, and
> \([x]K_2\ne0\) force Long's family, up to nonzero scaling and the sphere
> torus.  Hence they are equivalent to the all-order adjacent-jet
> phenomenon in this rectangle.

### Exact residual certificate

After six triangular reconstruction jets, normalize
\[
\begin{aligned}
 C&=1+x+Ex^2+Fx^3+Gx^4,\\
 D&=-1+B_1x+\cdots+B_6x^6.
\end{aligned}
\tag{7.18}
\]
Proposition 6.2 uniquely gives
\[
\begin{aligned}
B_1={}&-3,\\
B_2={}&-3-5E,\\
B_3={}&-1-12E-7F,\\
B_4={}&-3(3E+6F+4E^2+3G),\\
B_5={}&-2E-15F-24E^2-24G-36EF,\\
B_6={}&-4F-15E^2-21G-78EF-20E^3-48EG-27F^2.
\end{aligned}
\tag{7.19}
\]
Let \(P_7,\ldots,P_{10}\in\mathbb Q[E,F,G]\) be the remaining four pure
jets after removing their nonzero rational factors.  Exact rational
Gröbner reduction gives the ideal equality
\[
\begin{aligned}
\langle P_7,P_8,P_9,P_{10}\rangle
=\langle\,
 &E^3-4E^2-8G,\ G^2,\ FG,\\
 &-5E^2+3F^2-10G,\\
 &2E^2+EG+4G,\ 2E^2+3EF+4G
\,\rangle.
\end{aligned}
\tag{7.20}
\]
Every common zero has \(G=0\) from \(G^2=0\), then \(E=0\) from
\(2E^2+EG+4G=0\), and finally \(F=0\) from
\(-5E^2+3F^2-10G=0\).  Thus the residual scheme is supported only at
Long's point.  The exact checker verifies (7.19), constructs the four
residual jets, and verifies both containments in (7.20).
Theorem 3.1 supplies all remaining orders.

The scheme is nonreduced: its tangent space is spanned by \(E,F\), exactly
the two invisible directions predicted by Proposition 7.6 below.  This
Gröbner certificate proves the fixed rectangle but is not presented as a
uniform mechanism.

The next rectangle again meets the predicted cutoff.

> **Theorem 7.4 (twelve-jet rigidity in the \((5,7)\) rectangle).**  Let
\[
 \deg C\le5,\qquad \deg D\le7.
\tag{7.21}
\]
> Polynomiality, the pure identities through order twelve, and
> \([x]K_2\ne0\) force Long's family, up to nonzero scaling and the sphere
> torus.  Hence they are equivalent to the all-order adjacent-jet
> phenomenon in this rectangle.

After seven triangular jets, normalize
\[
 C=1+x+Ex^2+Fx^3+Gx^4+Hx^5.
\tag{7.22}
\]
Proposition 6.2 reconstructs
\[
\begin{aligned}
B_1={}&-3,\\
B_2={}&-3-5E,\\
B_3={}&-1-12E-7F,\\
B_4={}&-3(3E+6F+4E^2+3G),\\
B_5={}&-2E-15F-24E^2-24G-36EF-11H,\\
B_6={}&-4F-15E^2-21G-78EF-20E^3-48EG-27F^2-30H,\\
B_7={}&-3(12E^3+32E^2F+E^2+18EF+36EG+20EH\\
       &\hspace{20mm}{}+21F^2+24FG+2G+9H).
\end{aligned}
\tag{7.23}
\]
Let \(P_8,\ldots,P_{12}\) be the five residual pure jets after removing
nonzero rational factors.  Exact grevlex reduction followed by FGLM gives
the lexicographic ideal
\[
\begin{aligned}
\langle P_8,\ldots,P_{12}\rangle=\langle\,
&323471E^5+108000E^4+178200E^3-307800E^2F\\
&\quad+102600EF^2+583200EF-194400F^2+388800H,\\
&16751E^5-2700E^4-43200EF^2+24300G^2,\\
&-5584E^5-7425E^4-8100E^2F+13500EF^2+24300FG,\\
&-53551E^5-15525E^4+32400E^3+81000E^2F\\
&\quad-27000EF^2+145800EG+97200F^2,\\
&13E^5+120F^3,\ E^2(E^3+30F^2),\\
&E^3(21E^2+50F),\ E^6
\,\rangle.
\end{aligned}
\tag{7.24}
\]
The last generator forces \(E=0\) set-theoretically.  Then
\(13E^5+120F^3\) forces \(F=0\), the second generator forces \(G=0\), and
the first forces \(H=0\).  Hence the residual scheme is supported only at
Long's point.  Its tangent space is spanned by \(E,F,G\), of the predicted
dimension three.  The exact checker derives the residual jets and verifies
both ideal containments in (7.24).  Theorem 3.1 supplies all higher orders.

The next unclassified numerator rectangle is
\[
 \deg C\le6,\qquad \deg D\le8.
\tag{7.25}
\]

> **Open problem A.**  Classify (7.25) under polynomiality and the all-order
> adjacent-jet conditions.  More structurally, determine whether the
> rectangles \((d,d+2)\) are rigid for every fixed \(d\), and whether there
> is a uniform conceptual reason for the finite-jet cutoffs.

This remains strictly smaller than a general \(V_4\) classification.

The computed cutoffs suggest a precise statement.

> **Conjecture 7.5 (uniform rectangle rigidity).**  For every \(d\ge2\),
> polynomiality, the adjacent condition \([x]K_2\ne0\), and the pure jets
> \[
> [x^m]K_m=0\qquad(1\le m\le2d+2)
> \tag{7.26}
> \]
> force \(C=1+ax,\ D=-(1+ax)^3\), up to scaling and the sphere torus, in
> the rectangle \((\deg C,\deg D)\le(d,d+2)\).  The residual finite-jet
> scheme is supported at Long's point and has tangent dimension \(d-2\).

Theorems 7.1--7.4 prove the conjecture for \(d=2,3,4,5\); Theorem 6.1 gives
the slightly sharper cutoff three when \(d=1\).  Proposition 6.2 proves the
triangular reconstruction half for every \(d\).  What remains is a uniform
nonlinear argument for the residual \(C\)-parameters.

### 7.2 A possible uniform route

All \(t^2\)-linear height moments are packaged by
\[
\begin{aligned}
 \mathcal G(x,z)
 &=\sum_{m\ge0}K_m(x)z^m\\
 &=\int_0^1\frac{dv}{1-zC(x)-zD(x)v^2}\\
 &=\sum_{k\ge0}
   \frac{(zD(x))^k}{(2k+1)(1-zC(x))^{k+1}}.
\end{aligned}
\tag{7.27}
\]
The pure conditions are the vanishing diagonal
\([x^mz^m]\mathcal G=0\), while the adjacent detectors are the neighboring
diagonal.  The finite rectangle theorems say that, through \(d=3\), this
diagonal gap forces the denominator family to acquire the moving endpoint
factor \(1-t^2(1+ax)^2\).

A uniform proof could therefore avoid expanding individual moments: derive
the differential equation of \(\mathcal G\) in \(x,z\), impose the missing
diagonal, and show that its singular curve must have the repeated endpoint
contact of the one-profile family.  This is the most plausible
simplification of the growing finite-jet calculations.

There is an equivalent inverse-function formulation.  For each height
\(v\), put
\[
 f_v(x)=C(x)+v^2D(x)
\]
and let \(w_v(z)\) be the formal solution of
\[
 w_v=z f_v(w_v).
\]
Lagrange inversion gives the two identities
\[
\begin{aligned}
 1+\sum_{m\ge1}[x^m]f_v(x)^m z^m
 &=\frac{z\,w_v'(z)}{w_v(z)},\\
 \sum_{m\ge1}[x^{m-1}]f_v(x)^m z^m
 &=z\,w_v'(z).
\end{aligned}
\]
After height integration, the pure and adjacent generating series are
therefore
\[
\begin{aligned}
 \mathscr P(z)
 &=\int_0^1\frac{z\,w_v'(z)}{w_v(z)}\,dv,\\
 \mathscr A(z)
 &=\int_0^1z\,w_v'(z)\,dv.
\end{aligned}
\]
Since polynomiality gives \(f_v(0)=1-v^2\) after normalization, all pure
jets vanish if and only if
\[
 \boxed{
 \int_0^1
 \log\!\left(\frac{w_v(z)}{z(1-v^2)}\right)\,dv=0.}
\]
The expression is interpreted coefficientwise; the single endpoint
\(v=1\) is obtained by continuity and does not affect the integral.
Indeed, logarithmic differentiation of this identity gives
\(\mathscr P(z)=1\), and the constant is fixed at \(z=0\).

This packages the entire problem as an averaged inverse-map rigidity
statement.  A uniform proof of Conjecture 7.5 could try to show that the
boxed logarithmic identity forces the inverse branches to come from the
moving endpoint factor \(D=-C^3\).  The adjacent condition is then the
coefficientwise nonvanishing of the simpler average
\(\int_0^1z\,w_v'(z)\,dv\).

There is an important tangent obstruction to any such proof.

> **Proposition 7.6 (eventual tangent invisibility).**  At Long's normalized
> point
> \[
> C=X=1+x,\qquad D=-X^3,
> \tag{7.28}
> \]
> consider a first-order deformation
> \[
> C_\epsilon=X+\epsilon c(x),\qquad
> D_\epsilon=-X^3+\epsilon d(x),
> \tag{7.29}
> \]
> in the rectangle \((d,d+2)\), after fixing the constant and linear
> normalizations, so \(c\in x^2k[x]\).  The eventual linearized pure-jet
> equations force the coefficient of \(x^d\) in \(c\) to vanish, but they
> admit every \(c\) of degree at most \(d-1\) after a unique choice of the
> relevant high coefficients of \(d(x)\).  Thus the eventual system has
> \(d-2\) invisible \(C\)-directions.

Indeed, put
\[
\begin{aligned}
 A_m(X)&=\int_0^X(1-w^2)^{m-1}\,dw,\\
 B_m(X)&=\int_0^Xw^2(1-w^2)^{m-1}\,dw.
\end{aligned}
\tag{7.30}
\]
Their differences from their values at \(X=1\) have order at least \(m\).
Since \(c,d\) have no constant term, differentiating \(K_m\) at
\(\epsilon=0\) shows that its \(x^m\)-coefficient is, up to a nonzero
factor,
\[
 [x^m]X^{m-4}\bigl(d(x)+(2m+1)X^2c(x)\bigr),
\tag{7.31}
\]
where
\[
 \frac{B_m(1)}{A_m(1)}=\frac1{2m+1}.
\tag{7.32}
\]
For \(c_jx^j\), the resulting polynomial in \(m\) is
\[
 (2m+1)c_j\binom{m-2}{j-2},
\tag{7.33}
\]
of degree \(j-1\).  A term \(d_kx^k\) contributes
\[
 d_k\binom{m-4}{k-4},
\tag{7.34}
\]
of degree at most \(k-4\).  Since \(k\le d+2\), no \(d\)-term can cancel
the degree-\((d-1)\) contribution of \(c_dx^d\); hence \(c_d=0\).
Conversely the binomial polynomials in (7.34), for
\(4\le k\le d+2\), form a basis through degree \(d-2\), so they absorb
every contribution from \(c_2,\ldots,c_{d-1}\).

For \(d=3\), the single invisible direction is the nonreduced \(E\)-direction
in Theorem 7.2.  For \(d=4,5\), Theorems 7.3--7.4 exhibit respectively two
and three surviving directions.
Therefore a uniform rigidity proof cannot be a Jacobian-rank or first-order
argument; it must exploit nonlinear endpoint contact.

### 7.3 Nonvanishing of the detecting moments

Endpoint contact proves the pure cancellation independently of the values
\(C_m\).  Detection by \(x^\ell\) fails exactly at orders for which
\[
 C_m=\int_0^1R(v^2)^m\,dv=0.
\tag{7.35}
\]
For real nonzero \(R\), every even \(C_m\) is positive.  Hence only odd
orders can obstruct all-order detection.  Fixed-sign profiles such as
(4.2) avoid the obstruction, but sign-changing profiles may have isolated
odd cancellations.  This is a genuine obstruction, not just a caveat.

> **Proposition 7.7 (unique quadratic obstruction at each odd order).**
> For every odd \(n\ge1\), there is a unique \(c_n\in(0,1)\) such that
> \[
> R_{c_n}(z)=(1-z)(z-c_n)
> \tag{7.36}
> \]
> has the required simple endpoint zero but
> \[
> \int_0^1R_{c_n}(v^2)^n\,dv=0.
> \tag{7.37}
> \]
> Moreover,
> \[
> c_n\longrightarrow3-2\sqrt2
> \qquad(n\longrightarrow\infty,\ n\ {\rm odd}).
> \tag{7.38}
> \]

To see existence, regard the integral in (7.37) as a continuous function
of \(c\).  At \(c=0\), the integrand
\(((1-v^2)v^2)^n\) is nonnegative and not identically zero.  At \(c=1\),
it is \(-(1-v^2)^{2n}\), because \(n\) is odd.  The intermediate value
theorem gives \(c_n\).  Its derivative is
\[
 -n\int_0^1(1-v^2)^n(v^2-c)^{n-1}\,dv<0,
\tag{7.39}
\]
because \(n-1\) is even.  Hence the zero is unique.

At the zero, the positive and negative parts of the integral have equal
\(L^n\)-mass.  Their \(n\)-th roots converge to their respective suprema,
\[
 \max_{0\le z\le c}(1-z)(c-z)=c,\qquad
 \max_{c\le z\le1}(1-z)(z-c)=\frac{(1-c)^2}{4}.
\tag{7.40}
\]
Every accumulation point of \(c_n\) therefore solves
\(c=(1-c)^2/4\).  The unique solution in \((0,1)\) is
\(3-2\sqrt2\), proving (7.38).  For \(n=1\), direct integration gives
\[
 c_1=\frac15.
\tag{7.41}
\]

Consequently the hypothesis \(C_m\ne0\) in Theorem 3.1 cannot be removed
from the mixed-moment statement, even for quadratic \(R\) with the minimum
endpoint contact required by polynomiality.  The pure coefficient still
vanishes at the exceptional order; what disappears is its adjacent
detector.

> **Open problem B.**  Give an algebraic criterion on \(R\) for
> \(C_{2j+1}\ne0\) for every \(j\ge0\), or classify profiles for which one
> of these odd moments vanishes.

This problem concerns a one-variable pushforward measure and is separate
from sphere polynomiality.

### 7.4 Several profiles and odd height dependence

Allowing a sum of terms with different affine profiles destroys the single
substitution in (3.6).  Allowing odd powers of \(t\) also separates the two
Hopf endpoints \(t=1\) and \(t=-1\); polynomiality modulo \(t^2+2xy=1\)
then no longer reduces to one \(x\)-adic contact order.  Either extension
could support mechanisms not equivalent to Long's.

The relevant invariant is likely a contact filtration: the \(x\)-adic
orders of the numerator on the two endpoint branches, together with the
degree of the surviving Taylor polynomial after height integration.  The
one-profile theorem is the case where this filtration has one generator.

### 7.5 Simplified conceptual summary

Within the classified class, the mechanism needs only three facts:

1. polynomiality is endpoint contact \(s\ge r\);
2. height integration turns contact \(sm\) into a Taylor gap of length
   \(sm\);
3. the two adjacent coefficients force the remaining exponent \(q=r\).

The beta evaluation is needed only for the numerical value of \(C_m\), not
for pure-moment cancellation or uniqueness of the exponent.  This separates
the structural Hopf mechanism from the special closed form of Long's mixed
moments.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_hopf_lift_classification.py
```

The exact checker verifies the coefficient identities for several windings,
endpoint multiplicities, and all \(m\le20\), including a non-power profile.
It uses SymPy over \(\mathbb Q\) for the two ideal containments in the fixed
\((4,6)\) Gröbner certificate.  The finite checks are regressions;
Theorem 3.1 is the all-order proof for every surviving family.
