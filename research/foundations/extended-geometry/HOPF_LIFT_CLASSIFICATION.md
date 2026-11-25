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


Section 9 classifies nonlinear single-profile dilations over arbitrary characteristic-zero fields, using a classical interval-moment theorem for the converse. Sections 10–11 examine affine multiples of two known seeds; their finite necessary equations are paired with independent all-order converses.

## 9. Nonlinear single-profile dilations over characteristic-zero fields

This is a written classification of the specified family, over every field of characteristic zero. It does not classify arbitrary phases on the affine quadric or settle arbitrary primitive Image degree-five existence. All integrals below are polynomial integration with rational denominators. No finite list of experimental moments is used as an all-orders certificate.

### 9.1. The family and the conclusions

Let K be a field of characteristic zero and
\[
 A=K[x,y,t]/(t^2+2xy-1),\qquad U=1+x.
\]
Use the normalized functional
\[
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
 \tag{9.1}
\]
Fix integers r,q,k>=1, a nonzero scalar gamma in K, and a nonzero polynomial H in K[w] such that (1-w^2)^r divides H. Put
\[
 p=\gamma x^{-r}U^qH(tU^k),\qquad
 C_m=\frac12\int_{-1}^1H(w)^m\,dw.                 \tag{9.2}
\]
The displayed Laurent expression represents an element of A. Write H=(1-w^2)^r S, n=deg H, and s=deg S=n-2r. Indeed
\[
 B_k=\frac{1-t^2U^{2k}}x
 =2y-t^2\sum_{j=1}^{2k}\binom{2k}j x^{j-1},
 \qquad p=\gamma U^q B_k^r S(tU^k).                 \tag{9.3}
\]
Here a phase is **pure** if Lambda(p^m)=0 for every integer m>=1. It is **eventually pure** if this holds for all sufficiently large m. It is **safe against all fixed markers** if, for every B in A, Lambda(Bp^m)=0 for all sufficiently large m, with the cutoff allowed to depend on B.

**Theorem 9.1 (nonlinear-profile purity, safety and least degree).** For the nonzero family (9.2):

1. For every m>=1,
   \[
    \Lambda(p^m)=\gamma^m\binom{qm-k}{rm}C_m.       \tag{9.4}
   \]
   The binomial coefficient has its usual polynomial definition for an arbitrary integer upper argument.
2. Eventual purity holds exactly when q<=r. Exact purity holds exactly when q<=r and
   \[
     C_m=0\quad\text{for every }m\geq1\text{ with }qm<k.
                                                               \tag{9.5}
   \]
   Thus (9.5) is an exact finite reduction justified by the unrestricted zero tail (9.4).
3. Safety against all fixed markers holds exactly when q<r. In that branch a monomial marker x^a y^b t^h is zero at every
   \[
    m\geq\max\left\{1,
       \left\lceil\frac{k(h+2b+1)}q\right\rceil,
       \left\lceil\frac a{r-q}\right\rceil\right\}.
                                                               \tag{9.6}
   \]
4. If q=r, the single fixed marker x^k is unsafe. More precisely,
   \[
    \Lambda(x^kp^m)=\gamma^m C_m\quad(rm\geq k),     \tag{9.7}
   \]
   and the left side is zero when rm<k. Infinitely many of the values (9.7) are nonzero.
5. The least total degree of p in A, equivalently the least balanced degree of a homogeneous two-pair lift, is
   \[
     D=(k+1)n+q-r.                                \tag{9.8}
   \]
   Its degree-D lift is radially primitive. Every pure unsafe phase in this family has q=r and therefore composite least balanced degree
   \[
     D=(k+1)n,\qquad k+1\geq2,\quad n\geq2.        \tag{9.9}
   \]

Consequently this family has no radially primitive pure unsafe phase in any prime balanced degree, in particular five or seven. This conclusion allows arbitrary complex coefficients and arbitrary mixtures of height powers within the one fixed profile H. No nonpower assertion is imposed or needed for this exclusion; some members are proper powers.

### 9.2. The one classical nonvanishing input

For a nonzero complex polynomial Q, Müger--Tuset, *On the moments of a polynomial in one variable*, arXiv:1910.12148v1, Theorem 2.1, proves
\[
 \limsup_{m\to\infty}\left|\int_0^1Q(v)^m\,dv\right|^{1/m}>0.
                                                               \tag{9.10}
\]
In particular its moments cannot vanish eventually. Their Remark 2.2 credits the already known eventual-vanishing conclusion to François--Pakovich--Yomdin--Zhao, Corollary 4.1. Only that non-eventual-vanishing consequence is needed here.

It applies to H(2v-1), whose moments are C_m. For arbitrary K, let K_0 be the subfield generated over Q by the finitely many coefficients of H. This finitely generated characteristic-zero field embeds into C: send a transcendence basis to algebraically independent complex numbers and extend across its finite algebraic extension. Polynomial integration commutes with the embedding. If C_m vanished eventually in K, the embedded nonzero complex polynomial would contradict (9.10). Hence
\[
       C_m\ne0\quad\text{for infinitely many }m.    \tag{9.11}
\]
No ordering or positivity is asserted for K.

### 9.3. Exact endpoint formula and the purity directions

Define the polynomial
\[
 J_m(V)=\frac12\int_{-V}^VH(w)^m\,dw.
\]
Since H vanishes to order at least r at both endpoints,
\[
 J_m'(V)=\frac{H(V)^m+H(-V)^m}2
\]
is divisible by (V-1)^{rm}. Therefore
\[
 J_m(V)-C_m\quad\text{is divisible by }(V-1)^{rm+1}. \tag{9.12}
\]
Changing variables w=tU^k in polynomial integration gives
\[
 \Lambda(p^m)=\gamma^m[x^{rm}]U^{qm-k}J_m(U^k).     \tag{9.13}
\]
This identity can be read in K[[x]], since U is a unit there. In particular it remains valid when qm-k is negative. Since U^k-1 has x-order one, (9.12) says that replacing J_m(U^k) by C_m does not change the coefficient in (9.13). This proves (9.4) for every m, including all small orders.

If q<=r and qm>=k, then 0<=qm-k<rm, so the binomial in (9.4) is zero. If qm<k its upper argument is a negative integer and
\[
 \binom{qm-k}{rm}=(-1)^{rm}\binom{k-qm+rm-1}{rm}\ne0.
\]

a) Thus every q<=r phase is eventually pure, and it is exactly pure if and only if (9.5) holds.

b) If q>r, then for all sufficiently large m one has qm-k>=rm. The binomial is then a nonzero positive integer, interpreted in K. Eventual purity would force eventual vanishing of C_m, contrary to (9.11). This proves both excluded directions for q>r.

When q=r and rm>=k, multiplying by x^k changes coefficient extraction in (9.13) to [x^{rm-k}]. The endpoint error still starts above that degree, and
\[
 [x^{rm-k}]U^{rm-k}=1.
\]
This proves (9.7). If rm<k the Laurent integrand before coefficient extraction has a positive power of x times a power series in x, so its constant term is zero. Equation (9.11) now proves unsafety against x^k.

### 9.4. All fixed markers in the strict branch

Let q<r and let a,b,h be nonnegative integers. Set
\[
 N=rm-a+b,\qquad e=qm-k(h+2b+1),
\]
\[
 \begin{split}
 K_m(V)&=\frac12\int_{-V}^V
       w^h(V^2-w^2)^bH(w)^m\,dw,\\
 K_m^0(V)&=\frac12\int_{-1}^1
       w^h(V^2-w^2)^bH(w)^m\,dw.
 \end{split}                                      \tag{9.14}
\]
The fixed-endpoint polynomial K_m^0 has degree at most 2b in V. Substituting for y and then setting w=tU^k gives exactly
\[
 \Lambda(x^ay^bt^hp^m)
   =\gamma^m2^{-b}[x^N]U^eK_m(U^k).                \tag{9.15}
\]
The stronger endpoint estimate is
\[
 K_m(V)-K_m^0(V)\quad\text{divisible by }(V-1)^{rm+b+1}. \tag{9.16}
\]
For a direct proof, write V=1+u in the integral between 1 and V, and put w=1+v. Its integrand contains v^{rm}(u-v)^b times a polynomial, since H(1+v)^m is divisible by v^{rm} and V^2-w^2=(u-v)(2+u+v). Each integrated term contains u^{rm+b+1}. The other endpoint has the same estimate using w=-1-v. This is exact polynomial integration over any characteristic-zero field.

For m in (9.6), e>=0 and (r-q)m>=a, hence N>=0. Also N<=rm+b. We may therefore replace K_m by K_m^0 in (9.15), by (9.16) composed with U^k. The resulting polynomial has degree at most
\[
 e+2kb=qm-k(h+1),
\]
while its requested coefficient has degree N and
\[
 N-(qm-k(h+1))=(r-q)m-a+b+k(h+1)>0.
\]
Thus (9.15) is zero. A fixed polynomial marker is a finite sum of monomial markers, so taking the maximum of their cutoffs proves the forward safety direction. Conversely safety implies eventual purity by the marker 1, excluding q>r, and q=r was excluded by (9.7). This proves the exact safety statement of Theorem 9.1.

For a marker with a representative of total degree at most v, the simpler uniform cutoff m>=k(2v+1) suffices in the strict branch.

### 9.5. Least degree, explicit homogeneous lift and radial primitivity

Use the original two-pair bilinears
\[
 R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad
 W=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2,
 \qquad T^2+2ZW=R^2.
\]
Every bilinear has bidegree (1,1). Put L=R+Z and
\[
 \mathcal B_k=2WR^{2k}
 -T^2\sum_{j=1}^{2k}\binom{2k}j R^{2k-j}Z^{j-1}.
\]
Writing S(w)=sum_{j=0}^s s_j w^j, with s_s nonzero, define
\[
 F=\gamma L^q\mathcal B_k^r
       \sum_{j=0}^s s_jR^{(k+1)(s-j)}T^jL^{kj}.     \tag{9.17}
\]
Each term has balanced degree q+(2k+1)r+(k+1)s=D, and R=1 recovers (9.3).

The highest total-degree term of (9.3) is the single nonzero monomial
\[
 \gamma(-1)^rs_s x^{kn+q-r}t^n.
\]
It remains nonzero in the associated graded ring K[x,y,t]/(t^2+2xy): the irreducible quadratic t^2+2xy cannot divide a nonzero monomial. Hence the sphere relation cannot lower that degree, proving (9.8). This also proves the least balanced degree assertion by homogeneous dehomogenization.

At (xi_1,xi_2,z_1,z_2)=(1,1,1,-1) one has R=0,Z=-1,T=2,W=2. Modulo R,
\[
 \mathcal B_k=-T^2Z^{2k-1},\qquad
 F=\gamma(-1)^rs_s T^nZ^{kn+q-r}.
\]
The value at this point is nonzero. Thus R does not divide F: (9.17) is radially primitive.

For q=r, (9.8) becomes (9.9). The factor n is at least 2r>=2, and k+1>=2, so D cannot be prime. Raising the lift degree by homogeneous radial padding introduces a factor R and loses radial primitivity. This proves the prime-degree exclusion precisely within this family.

### 9.6. Original Image contractions and scope

For the usual two-pair contraction
\[
 \mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta,
\]
the normalized radial/angular identity for a balanced polynomial G of degree e is
\[
 \mathcal E_2(G)=(e+1)!\Lambda(G|_{R=1}).          \tag{9.18}
\]
This is the standard identity proved in *Two-pair Image Mathieu counterexample*, Sections 2 and 4. It can also be checked termwise on the bilinears above. Equations (9.4), (9.7) and (9.18) show that a pure q=r phase has
\[
 \mathcal E_2(F^m)=0\quad(m\geq1),\qquad
 \mathcal E_2(Z^kF^m)=(Dm+k+1)!\gamma^mC_m\quad(rm\geq k).
                                                               \tag{9.19}
\]
The fixed marker Z^k has infinitely many nonzero contractions. Unlike the real-coefficient subfamily, the argument need not give nonzero values at every even order.

For q<r, the fixed-marker conclusion extends to every original polynomial marker, including unbalanced ones. A bihomogeneous marker Q of bidegree (u,v) gives zero contraction automatically if u>v. If u<=v, the output is homogeneous in z of degree v-u. For each multi-index alpha with |alpha|=v-u,
\[
 \partial_z^\alpha\mathcal E_2(QF^m)
       =\mathcal E_2(\xi^\alpha QF^m).
\]
The marker on the right is balanced of degree v, so (9.6), or m>=k(2v+1), eventually makes every such scalar test zero. A homogeneous polynomial of degree v-u whose derivatives of that full order vanish is zero in characteristic zero. Taking a maximum over finitely many bihomogeneous parts proves the full eventual Image conclusion.

The family includes the prior k=1 real-profile calculation, now allowing arbitrary characteristic-zero coefficients, as well as the nonlinear dilation used to realize odd composite degrees. The known unsafe quartic and sextic profiles are examples. It does not contain every cubic-height coefficient family, every hyperbolic-coordinate family, or every polynomial on the quadric. In particular it does not contradict the prime-degree witnesses at every prime at least eleven, which use different coefficient architectures. It gives an exclusion sector for the remaining degree-five problem. The separate [joint-Hopf construction](PRIMITIVE_IMAGE_DEGREE_SEVEN_AND_HOPF_SUBALGEBRA.md) realizes degree seven outside this fixed-profile class.

The nonpower restriction is a separate condition. For example r=q=2,k=1,H=(1-w^2)^2 gives a proper square. The composite-degree exclusion above already holds without restricting to nonpowers.

### 9.7. Mathematical sources

- M. Müger and L. Tuset, [*On the moments of a polynomial in one variable*, arXiv:1910.12148v1](https://arxiv.org/abs/1910.12148), Theorem 2.1 and its proof, printed pages 2--4. Only non-eventual vanishing of the interval power moments is imported. Its Remark 2.2 gives the earlier FPYZ credit.
- [*Hopf lift classification*](HOPF_LIFT_CLASSIFICATION.md), Section 7, contains the k=1 real-profile endpoint and full-marker calculation. The argument here extends that calculation to arbitrary k and characteristic-zero coefficients, using the preceding classical interval theorem for the reverse directions.
- [*Odd primitive two-pair Image witnesses*](ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md), Sections 2--4, uses the q=r,k=2r odd-profile specialization. Its real even-moment positivity is stronger than (9.11) for that specialization.
- [*Two-pair Image Mathieu counterexample*](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md), Sections 2 and 4, supplies the bilinear radial/angular and original-marker contraction conventions used in (9.18).

## 10. Affine multiples of the quartic seed and primitive degree five

Work over any characteristic-zero field K, with
\[
 A_K=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
\]
Integration here means the rational polynomial functional. Define the known quartic seed
\[
 p_0=\frac{1+x}{x}\bigl(1-(1+x)^2t^2\bigr)
     =(1+x)\bigl(2y-(2+x)t^2\bigr).               \tag{10.1}
\]
For scalars a,b,c,d in K let
\[
 \ell=a+bx+cy+dt,\qquad f=p_0\ell.
\]
This is a concrete multiplicative two-profile attempt to leave the fixed-profile quartic mechanism. It allows the full affine-linear space on the sphere, including y and t, and usually has least degree five.

**Theorem 10.1 (four moments decide the quartic-affine family).** The following conditions are equivalent:

1. Lambda(f^m)=0 for every m>=1.
2. Lambda(f^m)=0 for 1<=m<=4.
3. b=c=d=0.

If a is nonzero, the resulting phase is the known unsafe quartic a*p_0, with a fixed marker x surviving at every positive order. If a=0 it is the zero phase. In particular no nonconstant affine multiplier of p_0 yields a pure phase, and no member of this degree-five presentation is simultaneously pure and radially primitive.

The equivalence is proved by a finite necessary-equation certificate and a separate unrestricted converse. It is not inferred from vanishing finitely many moments of an unexplained candidate. It says nothing about arbitrary quintics not divisible by p_0 on the sphere, and it does not classify arbitrary quintics or septics.

### 10.1. Exact moments used in the certificate

For nonnegative i,j,h define
\[
 M_m(i,j,h)=\Lambda(x^iy^jt^h p_0^m).
\]
If h is odd this is zero by height parity. If h is even, direct binomial expansion of (10.1), substitution for y, and polynomial integration give
\[
 M_m(i,j,h)=2^{-j}\sum_{v=0}^{m}(-1)^v\binom mv
 \binom{m+2v}{m-i+j}
 \sum_{w=0}^{j}\frac{(-1)^w\binom jw}{2v+h+2w+1}. \tag{10.2}
\]
A binomial with lower argument outside [0,m+2v] is interpreted as zero. Formula (10.2) is valid for every positive integer m; it has no bounded-search hypothesis. Multinomial expansion gives the homogeneous polynomial
\[
 \mathcal M_m=\Lambda(f^m)
 =\sum_{i+j+h\leq m}
 \frac{m!}{(m-i-j-h)!i!j!h!}
 a^{m-i-j-h}b^ic^jd^h M_m(i,j,h).                 \tag{10.3}
\]
Its first two instances are
\[
 \mathcal M_1=\frac{10b-3c}{15},                 \tag{10.4}
\]
\[
 \mathcal M_2=\frac2{105}
 (56ab+4ac+28b^2-4bc+3c^2+4d^2).                \tag{10.5}
\]
All quantities in (10.2)--(10.5) have rational coefficients, so the identities hold after scalar extension to every K under consideration.

Suppose the first four moments vanish. Equation (10.4) gives b=3c/10. If c=0, then b=0, and (10.5) gives 8d^2/105=0; hence d=0. There is no restriction on a.

It remains to exclude c!=0. Write z=a/c and v=d^2/c^2. Substituting b/c=3/10 into (10.5) gives
\[
  \mathcal M_2/c^2=\frac8{2625}(25v+130z+27),
  \qquad v=-\frac{130z+27}{25}.                  \tag{10.6}
\]
Expanding (10.2)--(10.3) for m=3,4 and making precisely these substitutions gives
\[
 \mathcal M_3/c^3=\frac8{625625}P(z),\qquad
 \mathcal M_4/c^4=\frac{256}{31906875}Q(z),       \tag{10.7}
\]
where
\[
 P(z)=32175z^2+56680z+14254,
\]
\[
 Q(z)=60775z^3+79220z^2+1802z-643.
\]
These small formulas can be derived using the finite sums (10.2)--(10.3), without computer algebra. The exact checker described in Section 12 independently compares these formulas with literal Laurent expansion.

### 10.2. The explicit inconsistency certificate

Set
\[
 U(z)=-1531111624900z^2-1247587548025z+290136844604,
\]
\[
 V(z)=810588507300z+1031834055265.
\]
Integer polynomial multiplication gives the exact Bezout identity
\[
             U(z)P(z)+V(z)Q(z)=3472141285450021. \tag{10.8}
\]
The nonzero integer on the right is invertible in every characteristic-zero field. Thus P and Q cannot both vanish at any z in K, or even in an algebraic closure. Equations (10.6)--(10.7) exclude c!=0. This proves that four vanishing moments force b=c=d=0.

For the converse let b=c=d=0. Put
\[
 c_m=\frac12\int_{-1}^1(1-w^2)^m dw
     =\frac{4^m(m!)^2}{(2m+1)!}\ne0.
\]
An endpoint calculation gives, for every polynomial H(x),
\[
 \Lambda(H(x)p_0^m)=c_m[x^m]H(x)(1+x)^{m-1}.    \tag{10.9}
\]
Indeed, change height variable w=(1+x)t. The endpoint primitive differs from its value at 1 by a multiple of x^{m+1}, so it does not contribute to the indicated coefficient. Taking H=1 and H=x gives
\[
 \Lambda((a p_0)^m)=0,\qquad
 \Lambda(x(a p_0)^m)=a^m c_m.                    \tag{10.10}
\]
This is the unrestricted pure/marked converse and completes the theorem.

### 10.3. Degree-five primitivity and the original Image statement

The highest filtration term of p_0 is -x^2t^2 of degree four. If (b,c,d) is nonzero, the highest part of ell is bx+cy+dt, a nonzero linear form. The associated graded sphere ring is
\[
       K[x,y,t]/(t^2+2xy),
\]
which is a domain. Hence the product of these highest terms is nonzero, and f has least degree exactly five. Its degree-five homogeneous lift is radially primitive.

Explicitly use the bilinears R,Z,W,T with T^2+2ZW=R^2 and define
\[
 F_4=(R+Z)\bigl(2WR^2-(2R+Z)T^2\bigr),\qquad
 L_1=aR+bZ+cW+dT,\qquad F_5=F_4L_1.             \tag{10.11}
\]
The first is the degree-four lift of p_0. If (b,c,d) is nonzero, neither F_4 nor L_1 is divisible by R, and R is irreducible in the original four-variable polynomial ring. Thus R does not divide F_5. Equivalently their nonzero leading sphere classes multiply nontrivially in the graded domain above.

The exact radial/angular contraction law gives
\[
 \mathcal E_2(F_5^m)=(5m+1)!\Lambda(f^m).
\]
Therefore the theorem excludes pure radially primitive Image phases of the form (10.11). If b=c=d=0 and a!=0, F_5=aR F_4 is pure but contains radial padding R; its primitive representative has degree four. The proper-power and fixed-marker conditions do not create an exception to this purity obstruction.

### 10.4. Source and proof boundaries

Long's actual arXiv:2609.16178v1, Section 3, Theorem 3.1/Corollary 3.2, proves the known quartic Hopf coefficient law (10.9), after the harmless bilinear normalization uv=2xy. Section 8 changes to square-root-free abelian coordinates; it does not produce a new primitive degree-five or degree-seven construction. The July paper arXiv:2607.19012v1, Theorem 2.1/Proposition 3.1, proves the related interval--Laurent cancellation. No priority claim is made for that mechanism or for the quartic seed.

Müger--Tuset's *The Mathieu conjecture for SU(2) reduced to an abelian conjecture*, DOI10.1016/j.indag.2023.10.001, supplies an integration formula and a conditional reduction (Theorem 2.5). Its mixed-support reverse statement is Conjecture 2.3; the matrix-coefficient reverse statement is explicitly a conjecture in Section 3. These are not mathematical premises here. Remark 2.6 records a generating-function continuation obstruction. No mixed purity-to-support theorem is being assumed.

The only calculations specific to the new restricted exclusion are the complete finite expressions (10.2)--(10.8), all displayed above. The separate exact-arithmetic verification is a check of those rational polynomial identities, not a formal proof system or a universal moment-prefix theorem.

## 11. Affine multiples of the sextic seed and primitive degree seven

Let K be any field of characteristic zero, and use
\[
 A_K=K[x,y,t]/(t^2+2xy-1),\qquad U=1+x,
\]
\[
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
\]
All integration is the rational polynomial functional. The known quartic and odd sextic seeds are
\[
 p_0=\frac Ux(1-U^2t^2),\qquad
 p_6=Ut p_0=\frac{U^2t}{x}(1-U^2t^2).
                                                        \tag{11.1}
\]
These represent polynomials on the sphere. Consider the full affine-linear multiplier
\[
     \ell=a+bx+cy+dt,\qquad f=p_6\ell,
     \qquad a,b,c,d\in K.                               \tag{11.2}
\]

**Theorem 11.1 (four orders decide the sextic-affine family).** The following are equivalent:

1. Lambda(f^m)=0 for every positive integer m.
2. Lambda(f^m)=0 at the four orders m=1,2,4,6.
3. b=c=d=0.

For a!=0 the remaining phase a*p_6 is the known pure unsafe sextic, with the fixed marker x nonzero at every even positive order. For a=0 the phase is zero. Every nonconstant affine ell gives f least degree seven and a radially primitive degree-seven lift, but such a phase is never pure. Hence this entire affine-multiplier construction supplies no primitive septic Image counterexample. Arbitrary degree-seven phases outside (11.2), as well as the global degree-five problem, remain outside the theorem.

The proof uses four necessary moment equations with an exact coprimeness certificate, followed by a separately proved all-orders converse. No finite moment prefix is extrapolated into an unrestricted vanishing claim.

### 11.1. Exact moment equations and parity

For i,j,h>=0 define T_m(i,j,h)=Lambda(x^i y^j t^h p_6^m). It is zero if m+h is odd. Otherwise elementary binomial expansion gives the all-orders formula
\[
 T_m(i,j,h)=2^{-j}\sum_{v=0}^{m}(-1)^v\binom mv
  \binom{2m+2v}{m-i+j}
  \sum_{w=0}^{j}\frac{(-1)^w\binom jw}{m+2v+h+2w+1}.
                                                        \tag{11.3}
\]
A binomial coefficient whose lower argument is outside its ordinary range is zero. To check (11.3), expand (1-U^2t^2)^m, substitute y=(1-t^2)/(2x), and extract x^{m-i+j} from U^{2m+2v}; polynomial integration of the remaining even height powers gives its last denominator.

Write M_m=Lambda(f^m). The usual multinomial formula with (11.3) computes each M_m as a rational homogeneous polynomial in a,b,c,d. The first equation is
\[
                   M_1=-\frac{2d}{15}.                  \tag{11.4}
\]
Thus every pure phase has d=0. Once d=0, f is odd in t, so every odd pure moment vanishes identically. We henceforth impose d=0. The second moment is
\[
 M_2=\frac2{3465}
 (264ab+44ac+132b^2-44bc+53c^2).                         \tag{11.5}
\]
Formula (11.3) determines the complete fourth and sixth moment polynomials explicitly; the exact checker in Section 12 independently compares them with literal Laurent expansion. The reductions needed for the proof are displayed next.

### 11.2. The c=0 branch and the denominator guard

If c=0, (11.5) becomes
\[
                 M_2=\frac8{105}b(2a+b).
\]
Either b=0, as required, or b=-2a. In the latter branch the fourth moment computed from (11.3) is
\[
                  M_4=-\frac{2048}{15015}a^4.
                                                        \tag{11.6}
\]
Its vanishing forces a=0 and hence b=0. Thus c=0 is completely settled.

Suppose instead c!=0 and set
\[
       z=b/c,\qquad \alpha=a/c,\qquad w=6z+1.
\]
Equation (11.5) gives
\[
       44(6z+1)\alpha+132z^2-44z+53=0.                 \tag{11.7}
\]
If w=0, then z=-1/6 and the final three terms in (11.7) equal 64, a contradiction in characteristic zero. Thus w!=0, and necessarily
\[
        \alpha=-\frac{132z^2-44z+53}{44(6z+1)}.
                                                        \tag{11.8}
\]
This accounts explicitly for the only denominator introduced in the elimination.

### 11.3. Two coprime univariate necessary equations

Substitute b/c=z and (11.8) into the homogeneous fourth and sixth moments from (11.3). Exact rational polynomial expansion gives
\[
 \frac{M_4}{c^4}=-\frac{4P(z)}{6455143695(6z+1)^3},\qquad
 \frac{M_6}{c^6}=\frac{2Q(z)}{48450701243889945(6z+1)^5},
                                                        \tag{11.9}
\]
where
\[
\begin{split}
 P(z)={}&2971558656z^7-5200227648z^6-13319044992z^5\\
        &+7148449616z^4+152847024z^3-485630508z^2\\
        &+167187968z+6731307,
\end{split}
\]
\[
\begin{split}
 Q(z)={}&596608834727854080z^{11}-2119204298356231680z^{10}\\
 &+804612135850617600z^9+4002560438537396736z^8\\
 &-1083906565862605056z^7-786969297304179264z^6\\
 &+423512030238845856z^5-13510803549968736z^4\\
 &-19554621670633104z^3+2571951562119874z^2\\
 &-446879603799373z+5677522608208.
\end{split}
\]
A concise explicit coprimeness certificate is available modulo 101. Coefficients in the following lists are in ascending degree order, starting with the constant term. Define
\[
 \overline U=[24,81,76,19,51,70,47,30,93,1,69],
 \qquad \overline V=[68,96,31,6,17,98,9].
\]
Direct integer polynomial multiplication and reduction modulo 101 gives
\[
              \overline U P+\overline V Q\equiv1\pmod{101}.
                                                        \tag{11.10}
\]
The leading coefficients of P and Q reduce respectively to 84 and 63, both nonzero. Thus (11.10) proves coprimeness over Q: if P and Q had a nonconstant common factor over Q, Gauss's lemma would give a primitive integer common factor. Its leading coefficient would be a 101-unit because it divides the leading coefficient of P. Its reduction would therefore remain nonconstant and divide both reductions, contradicting (11.10). Equivalently their integer Sylvester determinant is nonzero; its residue is 15 modulo 101.

Coprimeness over Q supplies a rational Bezout identity, which remains an identity over every characteristic-zero field. Hence P and Q have no common root in K or an algebraic closure. Equations (11.9) contradict M_4=M_6=0 when c!=0. Combined with Sections 11.1--11.2, this proves that the four necessary moment equations force b=c=d=0.

The displayed small certificate (11.10) and the leading-coefficient guard are sufficient for this proof.

### 11.4. The separate all-orders converse and original Image degrees

Put H(w)=w(1-w^2) and
\[
                  \kappa_m=\frac12\int_{-1}^1H(w)^m\,dw.
\]
Changing variable w=Ut in (11.1) gives
\[
 \Lambda(p_6^m)=[x^m]U^{m-1}J_m(U),\qquad
 J_m(V)=\frac12\int_{-V}^{V}H(w)^m\,dw.
\]
Since H has zeros at both endpoints, J_m(V)-kappa_m is divisible by (V-1)^{m+1}. Therefore, for every m>=1,
\[
 \Lambda(p_6^m)=\kappa_m[x^m]U^{m-1}=0,\qquad
 \Lambda(xp_6^m)=\kappa_m.                              \tag{11.11}
\]
Odd kappa_m vanish, and every even kappa_m is a nonzero rational number: its real polynomial integrand is nonnegative and not identically zero. This is a proof of rational nonvanishing, not an ordered-field hypothesis on K. Thus a*p_6 is pure for all orders, and unsafe for a!=0. This completes the theorem's unrestricted converse.

The highest filtration term of p_6 is -x^3t^3, of degree six. If (b,c,d) is nonzero, the highest term of ell is the nonzero linear form bx+cy+dt. Their product is nonzero in the associated graded domain K[x,y,t]/(t^2+2xy). Consequently f has least degree exactly seven.

Use the standard bilinears R,Z,W,T with T^2+2ZW=R^2. The explicit lift is
\[
 F_6=T(R+Z)^2\bigl(2R^2W-(2R+Z)T^2\bigr),\qquad
 F_7=F_6(aR+bZ+cW+dT).                                \tag{11.12}
\]
The nonzero degree-seven graded class shows that F_7 is radially primitive whenever (b,c,d) is nonzero. Alternatively R is irreducible in the original four-variable polynomial ring and divides neither factor in (11.12) in this case. The exact two-pair contraction law gives
\[
                   \mathcal E_2(F_7^m)=(7m+1)!M_m.
\]
The four necessary scalar equations therefore exclude pure primitive lifts of this form. When b=c=d=0, the degree-seven presentation aR F_6 has radial padding; its nonzero primitive representative has degree six. No proper-power condition can evade an obstruction to purity itself.

### 11.5. Existing classification and source boundaries

The seed (11.1), its all-orders even marked moments and primitive nonpower sextic lift are already retained in *Two-pair Image Mathieu counterexample*, Theorem 5.5, and *Odd cubic height: rational reconstruction and the sharp unsafe degree-six threshold*, Section 7. The latter's theorem addresses winding-one odd cubic-height phases, not arbitrary septic phases.

After d=0, a nonzero c in (11.2) introduces an x^-2 term of height five on the angular chart, so the present full affine-multiplier calculation is not merely the old winding-one odd cubic classification. The quartic-affine quintic exclusion concerns a different seed and different parity. The endpoint converse (11.11) is written above in full, and the new necessary-equation exclusion is exactly (11.3)--(11.10). No external support-converse conjecture is used.

The arithmetic certificate concerns the displayed rational polynomial identities and coprimeness; the unrestricted converse is the endpoint proof above. The separate [joint-Hopf construction](PRIMITIVE_IMAGE_DEGREE_SEVEN_AND_HOPF_SUBALGEBRA.md) realizes primitive degree seven; degree five remains open.

## 12. Verification scope of the profile and affine-multiplier extensions

Section 9 is a written proof. Its non-eventual interval-moment input is Müger–Tuset, Theorem 2.1, in the retained primary PDF; the field-transfer, endpoint divisibility, marker cutoff and degree calculations are given in full. The old real-profile and quadratic-height claims retain their earlier scopes and proof inputs.

The [exact affine-multiplier checker](verify_hopf_affine_multiplier_obstructions.py) constructs the quartic and sextic phases by literal Laurent multiplication, computes their relevant moment polynomials with rational arithmetic, and compares every coefficient with an independently expanded binomial formula. It checks the reduced necessary polynomials and the quartic integer and sextic modulo-101 Bézout certificates, including rejection of corrupted certificates. Its [receipt](hopf_affine_multiplier_obstructions_verification.json) records only these finite polynomial-identity checks. The all-order conclusions come from the printed endpoint converses and the necessary-equation implications.

These are internal mathematical and arithmetic checks. Formal verification, external mathematical review, independent replay assurance and worldwide priority are unclaimed. No classification outside the stated profile and affine-multiplier families, and no arbitrary degree-five exclusion, is inferred. Primitive degree seven is realized by the separate joint-Hopf construction.
