# Hopf-lift classification and minimal \(t^2\)-linear rigidity

## 1. Scope

This note isolates the mechanism behind the adjacent Taylor coefficients in
the two-pair Image-Mathieu counterexample. It classifies a parametric
one-profile class and the complete polynomial winding-one \(t^2\)-linear
numerator class, the latter over every characteristic-zero field. Its
formal reconstruction closes all degree rectangles in that class; it does
not search or classify the full space
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
parameters \(a_2,\ldots,a_d\). Section 8 below solves the unrestricted
formal reconstruction explicitly and classifies every polynomial pair
\(C,D\), with an exact ordinary moment-ideal presentation and sufficient
prefix 3d+2 on c_0!=0. This closes these higher-degree rectangles in
the stated class.

## 7. Real profiles with mixed height powers

**Status and scope.** The results in this section are internal deductions,
with no independent mathematical review, formalization or novelty claim.
They extend the endpoint argument to real polynomials in the signed height.
The real-coefficient hypothesis is essential to the nonvanishing argument
used here; no classification for complex profile coefficients is asserted.
The only proof inputs are the sphere relation (1.1), constant-term phase
averaging, polynomial integration and the radial contraction identity in
[the two-pair note, Sections 2 and 4](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
Neither a bounded search nor a new external theorem is used.

Let
\[
 A=\mathbb R[x,y,t]/(t^2+2xy-1),\qquad X=1+x,
\]
and write \(\Lambda\) for the normalized angular functional: substitute
\(y=(1-t^2)/(2x)\), take the constant term in \(x\), and then integrate
\(t\) against \(dt/2\) on \([-1,1]\). Fix
\[
 r\geq1,\quad q\geq1,\quad \gamma\in\mathbb R\setminus\{0\},\quad
 0\ne P\in\mathbb R[w],\qquad (1-w^2)^r\mid P(w),
\]
where \(r,q\) are integers, and set
\[
 p=\gamma x^{-r}X^qP(tX),\qquad
 C_m=\frac12\int_{-1}^1P(w)^m\,dw.
 \tag{7.1}
\]
This is a single affine phase profile, allowing arbitrary mixtures of
height powers inside \(P\). It represents a polynomial in \(A\): if
\(P(w)=(1-w^2)^rS(w)\), then
\[
 p=\gamma X^q\bigl(2y-(2+x)t^2\bigr)^r S(tX).
 \tag{7.2}
\]
In particular \(C_{2j}>0\) for every \(j\geq1\).

> **Theorem 7.1 (real signed-height profile classification).** For (7.1),
> \[
>  \Lambda(p^m)=\gamma^m\binom{qm-1}{rm}C_m
>  \qquad(m\geq1).
>  \tag{7.3}
> \]
> Thus all pure moments vanish exactly when \(q\leq r\). If \(q=r\),
> \[
>  \Lambda(xp^m)=\gamma^m C_m,
>  \tag{7.4}
> \]
> which is nonzero at every even \(m\). If \(q<r\), then for every fixed
> polynomial marker \(B\in A\),
> \(\Lambda(Bp^m)=0\) for all sufficiently large \(m\).
> The cutoff may depend on \(B\), as specified in Lemma 7.2.

To prove (7.3)--(7.4), use the polynomial
\[
 H_m(X)=X^{qm-1}J_m(X),\qquad
 J_m(X)=\frac12\int_{-X}^{X}P(w)^m\,dw.
 \tag{7.5}
\]
Substitution \(w=tX\) and phase extraction give
\[
 \Lambda(p^m)=\gamma^m[x^{rm}]H_m(1+x),\qquad
 \Lambda(xp^m)=\gamma^m[x^{rm-1}]H_m(1+x).
 \tag{7.6}
\]
Because \(P\) has zeros of order at least \(r\) at both endpoints,
\[
 J_m'(X)=\frac{P(X)^m+P(-X)^m}{2}
\]
is divisible by \((X-1)^{rm}\). Consequently
\[
 H_m(X)\equiv C_mX^{qm-1}\pmod{(X-1)^{rm+1}}.
 \tag{7.7}
\]
This proves (7.3). When \(q\leq r\), the binomial coefficient there is
zero. When \(q>r\), it is positive, and the even moments are nonzero
because \(\gamma^{2j}C_{2j}>0\). For \(q=r\), the coefficient of degree
\(rm-1\) in (7.7) proves (7.4). No assertion that odd marked moments
are nonzero is needed.

> **Lemma 7.2 (all fixed polynomial markers when \(q<r\)).** Suppose
> \(1\leq q<r\). For nonnegative integers \(a,b,h\),
> \[
>  \Lambda(x^ay^bt^h p^m)=0
>  \quad\text{if}\quad
>  m\geq\max\left\{1,
>       \left\lceil\frac{h+2b+1}{q}\right\rceil,
>       \left\lceil\frac{a}{r-q}\right\rceil\right\}.
>  \tag{7.8}
> \]
> Therefore \(m\geq2D+1\) suffices for every marker having a polynomial
> representative of total degree at most \(D\) in \(x,y,t\).

Here is the endpoint estimate needed for arbitrary markers. Put
\[
\begin{aligned}
 N&=rm-a+b,& e&=qm-h-2b-1,\\
 K_m(X)&=\frac12\int_{-X}^{X}
       w^h(X^2-w^2)^bP(w)^m\,dw,\\
 K_m^0(X)&=\frac12\int_{-1}^{1}
       w^h(X^2-w^2)^bP(w)^m\,dw.
\end{aligned}
 \tag{7.9}
\]
The superscript zero means fixed integration endpoints; it is not a
constant in \(X\). It is a polynomial of degree at most \(2b\).
After substituting for \(y\) and changing height variable,
\[
 \Lambda(x^ay^bt^h p^m)
   =\gamma^m2^{-b}[x^N](1+x)^eK_m(1+x).
 \tag{7.10}
\]
The endpoint difference has the stronger divisibility
\[
 K_m(X)-K_m^0(X)\quad\text{is divisible by}\quad
 (X-1)^{rm+b+1}.
 \tag{7.11}
\]
Indeed, split the difference into the integrals from \(1\) to \(X\)
and from \(-X\) to \(-1\). In the first, put
\(X=1+u\), \(w=1+v\). The factor \(P(1+v)^m\) is divisible by
\(v^{rm}\), while
\(X^2-w^2=(u-v)(2+u+v)\). Each term integrated from \(v=0\) to
\(v=u\) therefore has a factor \(u^{rm+b+1}\). The second integral
has the same estimate after putting \(w=-1-v\), using the zero of
\(P\) at \(-1\). These are polynomial integrations, so the estimate
is an exact divisibility statement.

Under (7.8), \(e\geq0\), \(N\geq0\), and \(N\leq rm+b\).
Thus (7.11) permits replacing \(K_m\) by \(K_m^0\) in (7.10).
The resulting polynomial has degree at most
\[
 e+2b=qm-h-1,
\]
whereas
\[
 N-(qm-h-1)=(r-q)m-a+b+h+1>0
\]
by (7.8). Its coefficient of degree \(N\) is zero, proving the lemma.
For \(a+b+h\leq D\), the simpler bound \(m\geq2D+1\) implies
both inequalities in (7.8). A polynomial marker is a finite linear
combination of these monomials. This completes Theorem 7.1.

### Minimal balanced degree and the original Image condition

Let \(n=\deg P\), \(s=n-2r=\deg S\), and write
\(S(w)=\sum_{j=0}^s s_jw^j\), where \(s_s\ne0\). Use the original
bilinears \(R,Z,W,T\), put \(L=R+Z\), and use
\(D=2WR^2-2RT^2-ZT^2\), as in the two-pair note.

> **Proposition 7.3 (primitive homogenization and degree parity).**
> The least balanced degree of a homogeneous lift of (7.1) is
> \[
>  d=2n+q-r.
>  \tag{7.12}
> \]
> Its lift in this degree is
> \[
>  F=\gamma L^qD^r
>       \sum_{j=0}^s s_jR^{2(s-j)}T^jL^j.
>  \tag{7.13}
> \]
> It is \(R\)-primitive. For \(q=r\), this is an Image failure
> witness of even balanced degree \(d=2n\), with marker \(Z\) and
> nonzero marked contractions at every even power. For \(q<r\), it
> obeys the eventual Image conclusion for every fixed original
> polynomial multiplier, including unbalanced multipliers. For a
> multiplier of bidegree \(u,v\), \(m\geq2v+1\) suffices.

Each term of (7.13) has balanced degree
\(q+3r+2s=2n+q-r\), and setting \(R=1\) recovers (7.2).
The highest total-degree term of (7.2), before sphere reduction, is
\[
 \gamma(-1)^rs_s x^{q+n-r}t^n.
\]
Its image is nonzero in the associated graded sphere ring
\(\mathbb R[x,y,t]/(t^2+2xy)\), which is a domain. This proves the
minimal degree assertion: the highest term cannot be removed by the
sphere relation. Equivalently, reduction using \(t^2=1-2xy\) gives
a nonzero term of that same degree in the unique representative of
height degree at most one.

At \((\xi_1,\xi_2,z_1,z_2)=(1,1,1,-1)\), one has
\(R=0\), \(T=2\), \(L=-1\), \(D=4\). Only the \(j=s\) term
in (7.13) remains, and its value is
\(\gamma s_s(-1)^{q+s}2^{s+2r}\ne0\). Thus \(R\nmid F\).
For \(q=r\), radial homogeneity and (7.4) give
\[
 \mathcal E_2(ZF^m)=(dm+2)!\,\gamma^m C_m,
 \tag{7.14}
\]
which is nonzero at every even \(m\), while all pure contractions
vanish. For \(q<r\), Lemma 7.2 first gives eventual zero contraction
against every balanced marker: such a marker is a polynomial in the
four bilinears, and its sphere restriction belongs to \(A\).

To check arbitrary original multipliers, split a multiplier \(Q\) into
its finitely many bihomogeneous parts, of bidegrees \((u,v)\). If
\(u>v\), then \(\mathcal E_2(QF^m)=0\) automatically, since the
number of differentiations exceeds the coordinate degree. If \(u\leq v\),
this contraction is a homogeneous polynomial in \(z\) of degree
\(v-u\). For each multi-index \(\alpha\) with \(|\alpha|=v-u\),
\[
 \partial_z^\alpha\mathcal E_2(QF^m)
     =\mathcal E_2(\xi^\alpha QF^m).
 \tag{7.15}
\]
The marker on the right is balanced. It eventually gives zero by the
preceding paragraph. There are only finitely many such scalar tests;
once all vanish, the homogeneous output polynomial is zero. Taking
the maximum cutoff over the finitely many parts of \(Q\) proves the
full eventual Image conclusion. More explicitly, the balanced markers
in (7.15) have balanced degree \(v\), so they are sums of products of
\(v\) bilinears. Their sphere restrictions have representatives of
total degree at most \(v\). Lemma 7.2 therefore gives the sufficient
cutoff \(m\geq2v+1\).

The classification therefore rules out primitive odd-degree Image
failures within (7.1): the only failing case is \(q=r\), whose least
balanced degree is even. Homogenization in a higher degree introduces
\(R\) padding and loses primitivity. This obstruction allows arbitrary
mixed height powers in a single real profile, but makes no claim about
multiple phase profiles or the full primitive degree-five space.


## 8. Complete complex polynomial quadratic-height classification

Here the coefficient field K is arbitrary of characteristic zero, with
no reality or positivity hypothesis. Use the same angular functional
\(\Lambda\) on \(t^2+2xy=1\). Let C,D be arbitrary polynomials in K[x]
such that

\[
 p=x^{-1}\bigl(C(x)+t^2D(x)\bigr)
 \tag{8.1}
\]

is a polynomial on the quadric. This polynomiality condition is exactly
C(0)+D(0)=0: eliminating t^2 rewrites (8.1) as

\[
 p=\frac{C+D}{x}-2yD.
\]

The condition is also necessary by restriction to x=0, where t^2=1.
This is the full winding-one quadratic-height numerator class, without
a degree bound on either C or D.

> **Theorem 8.1 (complete polynomial winding-one classification).**
> Every pure phase (8.1), meaning \(\Lambda(p^m)=0\) for all m>=1,
> belongs to exactly one of the following branches:
>
> 1. C and D are both divisible by x^2. The phase has strictly
>    positive circular weight and is safe against every fixed
>    polynomial marker.
> 2. \(C=c_0+c_1x\), with c_0!=0, and
>    \[
>    D=-\frac{C^3}{c_0^2}.
>    \tag{8.2}
>    \]
>    If c_1=0 the phase is 2c_0y and is safe. If c_1!=0 it is the
>    retained quartic Hopf profile, up to a nonzero scalar and
>    reciprocal rescaling of x,y. Its fixed marker x has a nonzero
>    moment at every positive power.
>
> Conversely every listed phase is pure. Hence the only unsafe phases
> in this entire class have least balanced degree four. In particular
> it contains no R-primitive Image failure in balanced degree five,
> seven, or any degree greater than four.

### 8.1. The unique formal reconstruction

Fix an arbitrary formal series C in K[[x]], with c_0=C(0)!=0.
The coefficients of a formal D with D(0)=-c_0 are uniquely determined
by the pure equations

\[
 \mu_m(C,D)=[x^m]\int_0^1(C+t^2D)^m\,dt=0,
 \qquad m\ge1.
 \tag{8.3}
\]

Indeed, the m-th equation uses only the coefficients through degree m.
Its coefficient of the new variable d_m=[x^m]D is

\[
 \tau_m(c_0)
 =mc_0^{m-1}\int_0^1t^2(1-t^2)^{m-1}\,dt
 =c_0^{m-1}\frac{2^{m-1}m!}{(2m+1)!!}\ne0.
 \tag{8.4}
\]

No term containing d_m^2 can contribute at degree m. Thus the
equations determine d_m successively for all m. This is the general
c_0 version of Proposition 6.2, and uses only rational polynomial
integration.

The formal solution has the following closed form:

\[
 \boxed{\quad D_*=-\frac{C^3}{(C-xC')^2}.\quad}
 \tag{8.5}
\]

To prove it, write B=C-xC' and X=C/B. These are formal units with
B(0)=c_0 and X(0)=1. Then C+t^2D_*=C(1-t^2X^2). Put

\[
 J_m(v)=\int_0^v(1-u^2)^m\,du,
 \qquad \kappa_m=J_m(1).
\]

Since J'_m has a zero of order m at 1 and X-1 belongs to (x),

\[
 \int_0^1(C+t^2D_*)^m\,dt
  =\frac{C^m}{X}J_m(X)
  \equiv\kappa_m C^{m-1}(C-xC')\pmod{x^{m+1}}.
 \tag{8.6}
\]

These are coefficientwise formal identities: the polynomial primitive
J_m can be composed with the unit series X, and every requested
coefficient depends on finitely many input coefficients. But

\[
 [x^m]C^{m-1}(C-xC')
   =[x^m]\left(C^m-\frac{x}{m}(C^m)'\right)=0.
 \tag{8.7}
\]

So (8.5) satisfies every formal pure equation. Uniqueness from (8.4)
proves that every pure formal pair with c_0!=0 has exactly this D.
No finite pure prefix is asserted sufficient for arbitrary formal C.

### 8.2. Polynomiality forces C to be affine

Now let C,D be polynomials, with C(0)!=0. Equality (8.5) is then an
identity in K(x), since its denominator is a unit at x=0. In particular

\[
                    (C-xC')^2\mid C^3\quad\hbox{in }K[x].
 \tag{8.8}
\]

Suppose h=deg C>=2, and work over an algebraic closure of K.
Then B=C-xC' has degree h: its leading coefficient is (1-h)
times that of C, which is nonzero in characteristic zero. Every root
alpha of C is nonzero. If its multiplicity in C is e, then its
multiplicity in B is exactly e-1, because the leading local term
of -xC' there is -alpha*e times the leading coefficient of C.

If C has k distinct roots, the total multiplicity of B at those roots
is therefore h-k<h. Since B has degree h, it has a root that is not
a root of C. This contradicts (8.8). Thus deg C<=1.
For C=c_0+c_1x, B=c_0 and (8.5) becomes exactly (8.2).
This root-count argument proves an unrestricted degree statement;
observing finitely many nonzero Taylor coefficients of D_* would not
suffice to establish it.

For the remaining branch C(0)=0, polynomiality also gives D(0)=0.
Write C=xC_1,D=xD_1. Then p=C_1(x)+t^2D_1(x) has only nonnegative
circular weights. If c=C_1(0),d=D_1(0), its pure moments are precisely
the moments of c+dt^2. The first two equations give

\[
 c+d/3=0,\qquad c^2+2cd/3+d^2/5=0.
\]

Substitution yields 4d^2/45=0. Hence c=d=0 over every
characteristic-zero field, so C and D are divisible by x^2.
Conversely those divisibilities give strictly positive circular weight,
which proves purity and eventual mixed vanishing directly.

### 8.3. The unsafe branch and the primitive-degree consequence

In branch (8.2) put u=c_1/c_0 and X=1+ux. Then

\[
 p=c_0x^{-1}X(1-t^2X^2).
\]

The endpoint calculation of Section 3 gives

\[
 \Lambda(xp^m)
   =c_0^m u^{m-1}\int_0^1(1-t^2)^m\,dt.
 \tag{8.9}
\]

It is nonzero for every m>=1 when u!=0. For u=0 the phase
is 2c_0y, which has one strict weight sign and is safe. For u!=0,
the highest spatial term of p is -c_0u^3 x^2t^2, nonzero in
the associated graded quadric ring, so its least balanced lift degree
is exactly four. Its degree-four lift is a scalar/rescaling of the
retained quartic; every higher-degree homogeneous lift introduces
a factor of R. Thus no unsafe phase in this class supplies a new
R-primitive higher-degree witness.

The formal reconstruction, its polynomial obstruction, and the zero
branch are self-contained algebraic arguments. The coefficient and
endpoint manipulations are the same classical polynomial-integration
mechanism used in Sections 3 and 6. The Image-space interpretation uses
the already stated original contraction and radial identities; no
new external classification theorem is imported.

### 8.4. Exact moment ideals, with c_0 inverted

The triangular reconstruction also gives an equality of ordinary ideals,
not just a classification of field-valued points. Fix degree bounds
deg C<=h and deg D<=k, and work over the universal coefficient ring

\[
 \mathscr R=\mathbb Q[c_0,c_0^{-1},c_1,\ldots,c_h,d_1,\ldots,d_k],
 \quad C=\sum_{i=0}^hc_ix^i,\quad
 D=-c_0+\sum_{i=1}^kd_ix^i.
\]

Write f_m=mu_m(C,D) for (8.3), B=C-xC', and

\[
 \delta_j=[x^j](D-D_*),\qquad
 E(x)=B(x)^2D(x)+C(x)^3,\qquad E_j=[x^j]E(x).
\]

The formal coefficients of D_* belong to \(\mathscr R\), since the
constant coefficient of B is the unit c_0. Also delta_0=E_0=0.

> **Theorem 8.2 (ordinary ideal presentation and finite prefix).**
> For every N>=1, as ordinary ideals of \(\mathscr R\),
> \[
> (f_1,\ldots,f_N)
>   =(\delta_1,\ldots,\delta_N)=(E_1,\ldots,E_N).
> \tag{8.10}
> \]
> In particular, for \(M=\max(2h+k,3h)\),
> \[
> (f_m:m\ge1)=(f_1,\ldots,f_M)=(E_1,\ldots,E_M).
> \tag{8.11}
> \]
> When M=0 all three ideals are zero. The equalities persist under
> every coefficient-algebra base change in which c_0 stays invertible.

To prove the first equality, substitute D=D_*+(D-D_*) in f_m.
The moment depends only on d_1,...,d_m, and its d_m coefficient is
the unit tau_m(c_0) in (8.4). Its specialization with every delta_j
zero is f_m(C,D_*)=0 by the universal formal identity (8.6)--(8.7).
Thus

\[
 f_m-\tau_m(c_0)\delta_m
                 \in(\delta_1,\ldots,\delta_{m-1}).
\]

One may first form this polynomial identity with independent coefficients
d_1,...,d_m, then specialize d_j=0 for j>k (and c_j=0 for j>h).
Thus the same identity holds when m exceeds the degree bound k;
delta_m is then a prescribed coefficient of -D_* rather than an
independent variable. Induction, and invertibility of tau_m(c_0),
prove equality of every finite prefix. This argument is in the coefficient ring itself; it
does not replace an ideal by its radical or infer membership from a
set-theoretic calculation.

For the second equality use the exact formal identity
E=B^2(D-D_*). Its coefficients satisfy

\[
 E_m=c_0^2\delta_m+
       \sum_{j=1}^{m-1}([x^j]B^2)\delta_{m-j}.
\]

The same triangular induction applies, since c_0^2 is a unit.
Finally E is a polynomial of degree at most M, so all E_j above
M are zero. Taking the union of (8.10) over all N proves (8.11),
including the automatic infinite tail. For the earlier rectangle
deg C<=d, deg D<=d+2, the explicit sufficient prefix is M=3d+2.

This theorem is restricted to the open c_0!=0 coefficient branch.
It preserves possible nonreduced structure and does not assert that
the coefficient ideal of E is radical, compute its nilpotence index,
or give this same prefix on the c_0=0 branch. It supplies the exact
ordinary ideal presentation of the higher residual equations left
by Proposition 6.2, while Theorem 8.1 gives their field-point
classification.

The [exact standard-library checker](verify_hopf_quadratic_height_classification.py)
constructs the formal solution and independently solves the triangular
moment equations on finite rational fixtures, checks both outputs and
all requested moments, verifies exact polynomial division in affine and
nonlinear examples, and includes zero-branch and corrupted-coefficient
controls. It also reconstructs moment polynomials in independent delta
variables and verifies explicit ordinary-ideal triangular certificates
on separate finite coefficient fixtures. Its
[receipt](hopf_quadratic_height_classification_verification.json)
records the checked finite orders. The all-degree conclusion follows
from the written root argument, and the ideal presentation and infinite
tail follow from the unit triangular identities, not Taylor-prefix
experiments.
The theorem received a separate internal reconstruction of (8.4)--(8.8)
and both branches. External human review, formal verification and
worldwide priority remain unestablished.

## Verification scope

The written classifications are the all-order arguments. The preexisting
[rank-one Gaussian checker](../scripts/verify_three_real_gmc_rank_one_classification.py)
is associated with the earlier classification material; it is not claimed
to verify Section 7. That section rests on its displayed elementary proof.
Its endpoint estimates, degree calculation and arbitrary-marker bridge
received focused internal AI proof checks. No retained computational
checker, independent replay, formal verification or external mathematical
review is claimed for Section 7. Finite checks do not replace the
universal proof.

Section 8 has its separately named checker and receipt, with the exact
finite scope described there. Those computations do not upgrade the
assurance claimed for Section 7 or classify phases outside (8.1).
