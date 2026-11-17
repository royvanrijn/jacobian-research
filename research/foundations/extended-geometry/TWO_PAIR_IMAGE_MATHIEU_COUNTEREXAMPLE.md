# A bidegree-\((4,4)\) counterexample in two contraction pairs

## 1. Statement

Work over a characteristic-zero field and use the contraction pairs

\[
 (\xi _1,z_1),\qquad(\xi _2,z_2).
\]

Put

\[
 \mathcal M_2
 =(\partial_{z_1}-\xi _1)k[\xi _1,\xi _2,z_1,z_2]
  +(\partial_{z_2}-\xi _2)k[\xi _1,\xi _2,z_1,z_2],
 \tag{1.1}
\]

and define the four bilinears

\[
\begin{aligned}
 R&=\xi _1z_1+\xi _2z_2,&
 Z&=\xi _1z_2,\\
 W&=2\xi _2z_1,&
 T&=\xi _1z_1-\xi _2z_2.
\end{aligned}
\tag{1.2}
\]

They satisfy the rank-one-quadric identity

\[
 T^2=R^2-2ZW.
\tag{1.3}
\]

Set

\[
 \boxed{
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right),
 \qquad Q=Z.}
\tag{1.4}
\]

> **Theorem 1.1.** For every \(m\geq1\),
> \[
>  \boxed{
>  \mathcal E_2(F^m)=0,\qquad
>  \mathcal E_2(QF^m)
>   =\frac{(4m+2)!\,m!}{(2m+1)!!}\ne0.}
> \tag{1.5}
> \]
> Consequently \(F^m\in\mathcal M_2\) but
> \(QF^m\notin\mathcal M_2\) for every \(m\geq1\).
> Thus \(\operatorname{SIC}(2)\) is false.

The one-pair Special Image Conjecture is known, so the minimum failing pair
dimension is exactly two.  The witness \(F\) is bihomogeneous of bidegree
\((4,4)\), has ordinary total degree eight, and has sixteen expanded terms.
The multiplier \(Q\) has bidegree \((1,1)\).

In the monomial bases
\(\xi _1^i\xi _2^{4-i}\) and \(z_1^jz_2^{4-j}\), the coefficient matrix of
\(F\), with \(i,j=0,\ldots,4\), is
\[
 \begin{pmatrix}
 -1&2&0&0&0\\
 -3/2&2&6&0&0\\
 -1/2&3/2&6&6&0\\
 0&1&3/2&2&2\\
 0&0&-1/2&-3/2&-1
 \end{pmatrix},
\qquad \det=48.
\tag{1.6}
\]
In particular \(F\) has tensor rank five, so this is a genuinely
nonseparable two-pair witness rather than a constant-coefficient GVC point
on the rank-one Segre cone.

## 2. Contraction as circular Gaussian expectation

For a monomial, let

\[
 \mathcal E_2(\xi^\alpha z^\beta)
 =\partial_z^\alpha z^\beta.
\tag{2.1}
\]

Zhao's image-kernel identity gives

\[
 \mathcal M_2=\ker\mathcal E_2.
\tag{2.2}
\]

Here the identity has a direct proof over the stated field. For every
polynomial \(H\),
\(\mathcal E_2(\xi_iH)=\partial_{z_i}\mathcal E_2(H)
=\mathcal E_2(\partial_{z_i}H)\), so \(\mathcal E_2\) annihilates
\(\mathcal M_2\). Conversely, reducing \(\xi_iH\) to
\(\partial_{z_i}H\) modulo the vector space \(\mathcal M_2\),
and inducting on total \(\xi\)-degree, gives
\(P-\mathcal E_2(P)\in\mathcal M_2\) for every \(P\). Thus a
polynomial with zero contraction lies in \(\mathcal M_2\), proving
(2.2) without an additional field-extension assumption.
The same reduction proves the image-kernel identity with any finite
number of contraction pairs.

Both \(F^m\) and \(QF^m\) have equal dual and coordinate degrees.  Their
contractions are therefore scalars.  If \(G_1,G_2\) are independent
standard circular complex Gaussians, the Wick rule gives, for every
balanced polynomial \(H\),

\[
 \mathcal E_2(H)
 =\mathbb E\,
 H(\overline G_1,\overline G_2,G_1,G_2).
\tag{2.3}
\]

All coefficients in (1.4) are rational, and both sides of (1.5) are
rational scalars obtained by formal differentiation.  It therefore
suffices to prove the identities over \(\mathbb C\) using (2.3); the
resulting rational identities then extend to every characteristic-zero
field.

Write

\[
 S=|G_1|^2+|G_2|^2,\qquad
 (G_1,G_2)=\sqrt S\,(U_1,U_2),
\tag{2.4}
\]

where \(U\) is uniform on the unit sphere in \(\mathbb C^2\) and is
independent of \(S\).  The variable \(S\) has the gamma distribution of
shape two, hence

\[
 \mathbb E(S^n)=(n+1)!.
\tag{2.5}
\]

On the sphere put

\[
 x=\overline U_1U_2,\qquad
 y=2\overline U_2U_1,\qquad
 t=|U_1|^2-|U_2|^2.
\tag{2.6}
\]

Then \(t^2+2xy=1\), and radial homogeneity gives

\[
 F=S^4p,\qquad Q=Sx,
\quad
 p=(1+x)\left(y-\frac12(2+x)t^2\right).
\tag{2.7}
\]

It remains to calculate two angular moments.

## 3. The all-order angular identity

The Hopf coordinates may be chosen so that \(t\) is uniform on
\([-1,1]\), the phase of \(x\) is uniform, and

\[
 xy=\frac{1-t^2}{2}.
\tag{3.1}
\]

Phase averaging is therefore constant-term extraction after substituting
\(y=(1-t^2)/(2x)\).  Equation (2.7) becomes

\[
 p=\frac{1+x}{2x}
 \left(1-t^2(1+x)^2\right).
\tag{3.2}
\]

Since the integrand is even in \(t\), define

\[
 H_m(X)
 =X^m\int_0^1(1-s^2X^2)^m\,ds.
\tag{3.3}
\]

The phase coefficient and the \(t\)-integral give

\[
\begin{aligned}
 \mathbb E_U(p^m)
 &=2^{-m}[u^m]H_m(1+u),\\
 \mathbb E_U(xp^m)
 &=2^{-m}[u^{m-1}]H_m(1+u).
\end{aligned}
\tag{3.4}
\]

After the change of variables \(v=sX\),

\[
 H_m(X)
 =X^{m-1}J_m(X),
\qquad
 J_m(X)=\int_0^X(1-v^2)^m\,dv.
\tag{3.5}
\]

Now \(J_m'(X)=(1-X^2)^m\) has a zero of order \(m\) at \(X=1\).
Consequently the Taylor coefficients through order \(m\) of
\(H_m(X)\) at \(X=1\) are the same as those of
\(J_m(1)X^{m-1}\).  Equation (3.4) immediately yields

\[
 \mathbb E_U(p^m)=0,\qquad
 \mathbb E_U(xp^m)=2^{-m}J_m(1).
\tag{3.6}
\]

The remaining beta integral is

\[
\begin{aligned}
 J_m(1)
 &=\int_0^1(1-v^2)^m\,dv\\
 &=\frac12B\left(\frac12,m+1\right)
 =\frac{2^m m!}{(2m+1)!!}.
\end{aligned}
\tag{3.7}
\]

Thus

\[
 \mathbb E_U(p^m)=0,\qquad
 \mathbb E_U(xp^m)=\frac{m!}{(2m+1)!!}.
\tag{3.8}
\]

The endpoint-zero mechanism in this calculation is classified in
[`HOPF_LIFT_CLASSIFICATION.md`](HOPF_LIFT_CLASSIFICATION.md).  Inside the
one-profile class treated there, the phase winding forces the profile
exponent, and the present five-term polynomial is the unique member of
minimum profile degree, up to the sphere torus and overall scaling.

Combining (2.3), (2.5), (2.7), and (3.8) gives

\[
\begin{aligned}
 \mathcal E_2(F^m)
 &=(4m+1)!\,\mathbb E_U(p^m)=0,\\
 \mathcal E_2(QF^m)
 &=(4m+2)!\,\mathbb E_U(xp^m)
 =\frac{(4m+2)!\,m!}{(2m+1)!!},
\end{aligned}
\]

which proves Theorem 1.1.

## 4. A second proof by direct coefficient extraction

This section gives an algebraic proof independent of the circular-Gaussian
and Hopf-coordinate argument above.  In particular, it uses neither a
Gaussian measure nor radial decomposition nor the Taylor-jet argument in
(3.5)--(3.6).

Let \(I:\mathbb Q[t]\to\mathbb Q\) be the linear functional
\[
 I(t^j)=\frac1{j+1}.
\]
It satisfies the formal fundamental theorem
\(I(P')=P(1)-P(0)\), first on monomials and hence on every polynomial.
For \(a\geq1\), applying this to
\(P=t^a(1-t)^{b+1}\), whose two endpoint values vanish, gives
\[
 (b+1)I\!\left(t^a(1-t)^b\right)
 =aI\!\left(t^{a-1}(1-t)^{b+1}\right).
\]
Induction on \(a\), starting from
\(I((1-t)^b)=1/(b+1)\), proves the entirely algebraic beta identity
\[
 I\!\left(t^a(1-t)^b\right)=\frac{a!\,b!}{(a+b+1)!}.
 \tag{4.1}
\]
We continue to write \(I(P)=\int_0^1P(t)\,dt\), but no analytic integration
is involved.  For a polynomial \(H\) bihomogeneous of bidegree \((n,n)\),
direct monomial inspection now gives the coefficient-extraction identity
\[
 \boxed{\displaystyle
 \mathcal E_2(H)
  =(n+1)!\operatorname {CT}_u
    \int_0^1H\left(1,u,t,\frac{1-t}{u}\right)\,dt.}
 \tag{4.2}
\]
Indeed, for
\[
 H=\xi _1^a\xi _2^{n-a}z_1^bz_2^{n-b},
\]
the constant term in \(u\) is zero unless \(a=b\); when \(a=b\), (4.1)
turns the right-hand side of (4.2) into \(a!(n-a)!\), exactly its
contraction.

Under the substitution in (4.2), put
\[
 v=2t-1,\qquad x=Z=\frac{1-v}{2u}.
\]
Then
\[
 R=1,\qquad T=v,\qquad ZW=\frac{1-v^2}{2},
\]
and hence
\[
 F=\frac{1+x}{2x}\left(1-v^2(1+x)^2\right),\qquad Q=x.
 \tag{4.3}
\]
Here is the precise constant-term change.  Before applying \(I\), work over
the coefficient field \(\mathbb Q(v)\).  If \(c=(1-v)/2\), the substitution
\(x=c/u\) sends \(u^j\) to \(c^jx^{-j}\).  Thus exponent zero is preserved
and no nonzero exponent can contribute to it:
\(\operatorname {CT}_u=\operatorname {CT}_x\).  After simplification,
(4.3) lies in \(\mathbb Q[v][x,x^{-1}]\), so the formal functional \(I\)
applies coefficientwise.  The expression is even in \(v\).  Expanding it
in (4.2) is especially short:
\[
\begin{aligned}
 \operatorname {CT}_x F^m
 &=2^{-m}\sum_{k=0}^m(-1)^k\binom mk
   \operatorname {CT}_x\!
    \left(x^{-m}(1+x)^{m+2k}\right)v^{2k}\\
 &=2^{-m}\sum_{k=0}^m(-1)^k\binom mk
   \binom{m+2k}{m}v^{2k},\\
 \operatorname {CT}_x(xF^m)
 &=2^{-m}\sum_{k=0}^m(-1)^k\binom mk
   \binom{m+2k}{m-1}v^{2k}.
\end{aligned}
\]
Moreover the formal fundamental theorem already used above gives
\[
 I\!\left((2t-1)^{2k}\right)
 =I\!\left(\left(\frac{(2t-1)^{2k+1}}
 {2(2k+1)}\right)'\right)
 =\frac1{2k+1}.
\]
Since \(F^m\) and \(QF^m\) have bidegrees \((4m,4m)\) and
\((4m+1,4m+1)\), respectively, (4.2) now gives, for every \(m\geq1\),
\[
\begin{aligned}
 \frac{\mathcal E_2(F^m)}{(4m+1)!}
  &=2^{-m}\sum_{k=0}^m
    \frac{(-1)^k\binom mk\binom{m+2k}{m}}{2k+1},
    \tag{4.4}\\
 \frac{\mathcal E_2(QF^m)}{(4m+2)!}
  &=2^{-m}\sum_{k=0}^m
    \frac{(-1)^k\binom mk\binom{m+2k}{m-1}}{2k+1}.
    \tag{4.5}
\end{aligned}
\]
These formulas came directly from formal contraction, not from Wick
expectation or phase averaging.

We now evaluate both sums by finite differences.  If
\(\Delta A(X)=A(X+1)-A(X)\), induction on \(m\) gives
\[
 \sum_{k=0}^m(-1)^k\binom mk A(k)=(-1)^m\Delta^mA(0).
\]
Each application of \(\Delta\) lowers the degree of a nonconstant
polynomial by one.  Consequently, for every polynomial \(A(X)\) of degree
less than \(m\),
\[
 \sum_{k=0}^m(-1)^k\binom mk A(k)=0.
 \tag{4.6}
\]
For the pure sum,
\[
 \frac{\binom{m+2X}{m}}{2X+1}
   =\frac1{m!}\prod_{j=2}^{m}(2X+j)
 \tag{4.7}
\]
is a polynomial of degree \(m-1\).  Equations (4.4), (4.6), and (4.7)
immediately prove \(\mathcal E_2(F^m)=0\).

The same observation gives a useful general residue principle.  If
\(\deg A<m\), polynomial division at \(X=-1/2\) gives
\[
 A(X)=A(-1/2)+(2X+1)D(X),\qquad \deg D<m.
\]
Consequently
\[
 \boxed{\displaystyle
 \sum_{k=0}^m(-1)^k\binom mk\frac{A(k)}{2k+1}
 =A(-1/2)
   \sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+1}.}
 \tag{4.8}
\]
Thus the alternating quotient transform on all polynomials of degree
less than \(m\) has rank one: it is evaluation at \(-1/2\), multiplied by
one universal scalar.

There is a useful higher-degree-denominator form of the same argument.  Fix
a polynomial \(L\) such that \(L(k)\ne0\) for \(0\leq k\leq m\), and write
\[
 {\cal T}_{m,L}(A)
 =\sum_{k=0}^m(-1)^k\binom mk\frac{A(k)}{L(k)}.
\]
If
\[
 A=R+LD,\qquad \deg D<m,
\]
then (4.6), applied to \(D\), gives the exact invariance
\[
 \boxed{{\cal T}_{m,L}(A)={\cal T}_{m,L}(R).}
 \tag{4.8a}
\]
In particular, if \(L\) is normalized to be monic of degree \(r\geq1\)
and \(\deg A<m+r\), polynomial division gives \(\deg D<m\) and
\[
 {\cal T}_{m,L}(A)={\cal T}_{m,L}(A\bmod L).
 \tag{4.8b}
\]
Thus, on this degree range, the transform factors through the
\(r\)-dimensional quotient \(\mathbb Q[X]/(L)\): it sees only the remainder
modulo \(L\).  (As a single scalar-valued functional its linear-map rank is
of course at most one.)  Formula (4.8) is the case \(r=1\).  For a repeated
denominator
\(L=(X+\tfrac12)^r\), the remainder is equivalently determined by the
first \(r\) Taylor coefficients of \(A\) at \(-1/2\).  This jet
interpretation is a structural extension of the finite-sum certificate;
by itself it does not assert the existence of further SIC witnesses.

The repeated-pole case admits an explicit all-order evaluation.  For
\(s\geq0\), put
\[
 B_{m,s}=\sum_{k=0}^m
 \frac{(-1)^k\binom mk}{(2k+1)^s}.
\]
Thus \(B_{m,0}=0\) for \(m\geq1\), and the sum \(B_m\) in (4.10) is
\(B_{m,1}\).  For \(m,s\geq1\), the termwise identity
\[
 \frac{2m+1}{(2k+1)^s}
 =\frac1{(2k+1)^{s-1}}
  +\frac{2(m-k)}{(2k+1)^s}
\]
and \((m-k)\binom mk=m\binom{m-1}k\) give the triangular recurrence
\[
 \boxed{(2m+1)B_{m,s}=B_{m,s-1}+2mB_{m-1,s}.}
 \tag{4.8c}
\]
The boundary values are \(B_{0,s}=1\) and \(B_{m,0}=0\) for \(m\geq1\).
Hence (4.8c) evaluates every repeated-pole residue using only rational
arithmetic.

More explicitly, suppose \(r\geq1\), \(\deg A<m+r\), and expand the
remainder in the local coordinate \(\lambda=2X+1\):
\[
 A(X)=\sum_{j=0}^{r-1}c_j\lambda^j+\lambda^rD(X),
 \qquad \deg D<m.
\]
Then (4.6) gives the exact jet formula
\[
 \boxed{\displaystyle
 \sum_{k=0}^m(-1)^k\binom mk
 \frac{A(k)}{(2k+1)^r}
 =\sum_{j=0}^{r-1}c_jB_{m,r-j}
 =\sum_{j=0}^{r-1}
 \frac{A^{(j)}(-1/2)}{2^j j!}\,B_{m,r-j}.}
 \tag{4.8d}
\]
There is also a compact closed description.  With an auxiliary variable
\(a\), the same first-order recurrence as (4.11) gives
\[
 G_m(a):=\sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+a}
 =\frac{2^m m!}{\prod_{\ell=0}^m(a+2\ell)}.
\]
Formal differentiation in \(\mathbb Q(a)\) yields
\[
 B_{m,s}
 =\frac{(-1)^{s-1}}{(s-1)!}
 \left.\frac{d^{\,s-1}}{da^{s-1}}G_m(a)\right|_{a=1}
 \qquad(s\geq1).
 \tag{4.8e}
\]
If \(h_q\) denotes the complete homogeneous symmetric polynomial of degree
\(q\), expansion at \(a=1\) makes this still more concrete:
\[
 \boxed{\displaystyle
 B_{m,s}
 =\frac{2^m m!}{(2m+1)!!}\,
 h_{s-1}\left(1,\frac13,\ldots,\frac1{2m+1}\right)>0.}
 \tag{4.8f}
\]
Indeed,
\[
 \frac{G_m(1+y)}{G_m(1)}
 =\prod_{\ell=0}^m
  \left(1+\frac{y}{2\ell+1}\right)^{-1},
\]
whose coefficient of \(y^q\) is
\((-1)^qh_q(1,1/3,\ldots,1/(2m+1))\).  Consequently the jet in (4.8d)
cannot vanish when its local coefficients \(c_j\) are all nonnegative and
at least one is positive.  Equations (4.8c)--(4.8f) make the surviving
\(r\)-jet effective; realizing such a jet by an admissible SIC multiplier
remains a separate problem.

For the mixed sum, take
\[
 A_m(X)=\binom{m+2X}{m-1}
       =\frac1{(m-1)!}\prod_{j=2}^{m}(2X+j).
 \tag{4.9}
\]
It has degree \(m-1\) and \(A_m(-1/2)=1\).  Formula (4.8) therefore
reduces the mixed sum to
\[
 B_m=\sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+1}.
 \tag{4.10}
\]
This last sum has a purely finite recurrence.  Since
\[
 m\binom{m-1}{k}=(m-k)\binom mk,
\]
termwise subtraction gives
\[
\begin{aligned}
 (2m+1)B_m-2mB_{m-1}
 &=\sum_{k=0}^m(-1)^k\binom mk\\
 &=0.
\end{aligned}
\tag{4.11}
\]
Together with \(B_0=1\), this proves
\[
 B_m=\prod_{j=1}^m\frac{2j}{2j+1}
     =\frac{2^m m!}{(2m+1)!!}.
 \tag{4.12}
\]
Substitution in (4.5) yields
\[
 \mathcal E_2(QF^m)
  =(4m+2)!\,2^{-m}B_m
  =\frac{(4m+2)!\,m!}{(2m+1)!!}.
\]
This proves both identities in Theorem 1.1 for all \(m\) a second time.
The only general inputs are the monomial contraction rule, the formal beta
identity (4.1), and the elementary finite-difference identity (4.6).

## 5. Propagation to every balanced degree \(d\geq4\)

The quartic witness propagates without changing its angular cancellation.
The mechanism is multiplication by powers of the invariant \(R\).

> **Lemma 5.1 (radial shift).** If \(H\) is bihomogeneous of bidegree
> \((n,n)\), then for every \(a\geq0\),
> \[
>  \mathcal E_2(R^aH)
>  =\frac{(n+a+1)!}{(n+1)!}\mathcal E_2(H).
>  \tag{5.1}
> \]

There are again two proofs.  On circular Gaussians, \(R=S\), so radial
independence gives (5.1).  Directly and algebraically, the substitution in
(4.2) sends \(R\) to \(1\).  Formula (4.2), applied in bidegrees \(n\) and
\(n+a\), therefore gives
\[
 \mathcal E_2(R^aH)
 =(n+a+1)!\operatorname {CT}_u\int_0^1
 H\left(1,u,t,\frac{1-t}{u}\right)\,dt
 =\frac{(n+a+1)!}{(n+1)!}\mathcal E_2(H).
\]
Thus the propagation also has a proof over an arbitrary
characteristic-zero field with no Gaussian input.

For \(d\geq4\), put
\[
 \boxed{F_d=R^{d-4}F,\qquad Q=Z.}
 \tag{5.2}
\]
Then \(F_d\in V_d\).  Apply Lemma 5.1 first to \(H=F^m\), which has
bidegree \((4m,4m)\), with \(a=(d-4)m\).  Apply it again to \(H=QF^m\),
which has bidegree \((4m+1,4m+1)\).  Theorem 1.1 gives
\[
\begin{aligned}
 \mathcal E_2(F_d^m)&=0,\\
 \mathcal E_2(QF_d^m)
 &=\frac{(dm+2)!}{(4m+2)!}\mathcal E_2(QF^m)\\
 &=\boxed{\frac{(dm+2)!\,m!}{(2m+1)!!}}\ne0
 \qquad(m\geq1).
\end{aligned}
\tag{5.3}
\]

**Geometric verification boundary.** The nullcone and semistability
interpretation in the next statement and Section 6 still requires a
specified group action and an audited implication from its nullcone to
eventual fixed-marker vanishing. This is recorded as a verification
obligation in the dependency graph. That geometric interpretation is
conditional on this input; the contraction identities and Image failures
proved here do not depend on it.

> **Theorem 5.2 (geometric input pending verification).** The moment--nullcone assertion
> \(\mathrm{MN}_d\) fails for every \(d\geq4\).

The fixed multiplier \(Q\) is independent of \(d\), and every mixed moment
in (5.3) is nonzero.  Since a one-sided nullcone point has eventual mixed
vanishing for every fixed multiplier, \(F_d\) is outside the nullcone.
Thus (5.3) supplies a moment-zero semistable point in every balanced degree
\(d\geq4\).

The radial order can be bounded independently of \(d\) by combining
ordinary powers with invariant multiplication.  More generally, let
\[
 d=4r+k,\qquad r\geq1,\quad k\geq0,
\]
and put
\[
 G_{r,k}=R^kF^r.
 \tag{5.4}
\]
The same calculation, now starting from \(F^{rm}\) and \(QF^{rm}\), gives
\[
\boxed{
 \mathcal E_2(G_{r,k}^m)=0,\qquad
 \mathcal E_2(QG_{r,k}^m)
 =\frac{((4r+k)m+2)!\,(rm)!}{(2rm+1)!!}\ne0.}
 \tag{5.5}
\]
Thus for \(r=\lfloor d/4\rfloor\) and \(k=d-4r\), every degree \(d\geq4\)
has a witness whose \(R\)-adic order is at most three.  Moreover \(F\) is
not divisible by \(R\): at
\[
 (\xi _1,\xi _2,z_1,z_2)=(1,1,1,-1)
\]
one has \(R=0\) and \(F=-2\).  Since \(R\) is irreducible, the exact
\(R\)-adic order of \(G_{r,k}\) is \(k\).  In particular:

> **Corollary 5.3.** Every degree divisible by four has an
> \(R\)-primitive counterexample \(F^{d/4}\), and every other degree
> \(d\geq4\) has one divisible by no more than \(R^3\).

The primitive family in degrees divisible by four need not consist of
proper powers.  Put
\[
 D=2WR^2-2T^2R-ZT^2
 \tag{5.6}
\]
and, for \(h\geq1\), define
\[
 \Phi_h=(R+Z)D
 \sum_{j=0}^{h-1}\binom{h-1}{j}
 R^{4(h-1-j)}T^{2j}(R+Z)^{2j}.
 \tag{5.7}
\]
Every summand has balanced degree \(4h\).  On the Hopf sphere this is the
one-profile polynomial from the
[Hopf-lift classification](HOPF_LIFT_CLASSIFICATION.md)
\[
 x^{-1}(1+x)\rho_h\!\left(t^2(1+x)^2\right),
 \qquad
 \rho_h(z)=(1-z)(1+z)^{h-1}.
 \tag{5.8}
\]
The endpoint-contact argument gives
\[
\begin{aligned}
 \mathcal E_2(\Phi_h^m)&=0,\\
 \mathcal E_2(Q\Phi_h^m)
 &=(4hm+2)!\,C_{h,m},\\
 C_{h,m}
 &=\int_0^1(1-v^2)^m(1+v^2)^{(h-1)m}\,dv>0.
\end{aligned}
\tag{5.9}
\]
The integral is rational by direct polynomial integration.  At the point
used above, \(R=0\), \(Z=-1\), and \(T=2\); only the \(j=h-1\) summand
survives and is nonzero.  Hence \(\Phi_h\) is \(R\)-primitive.  The factor
\(R+Z\) has multiplicity exactly one: neither \(D\) nor the sum in (5.7)
vanishes identically on \(R+Z=0\).  Thus \(\Phi_h\) is not a proper power.
Notice that \(\Phi_1=2F\).

> **Theorem 5.4.** Every balanced degree divisible by four admits an
> explicit \(R\)-primitive, non-proper-power failure of
> \(\mathrm{MN}_d\).

The congruence restriction is intrinsic to the specified one-profile class
depending on height through \(t^2\), not just to (5.8).  If the profile
polynomial has degree \(h\), its leading
term in
\[
 c\,x^{-r}(1+ax)^r\rho\!\left(t^2(1+ax)^2\right)
\]
is a nonzero scalar multiple of \(x^{2h}t^{2h}\), of total angular degree
\(4h\).  Endpoint contact makes the full Laurent expression polynomial and
introduces no term of larger degree.  Hence its minimal balanced
homogenization has degree \(4h\); placing it in any larger \(V_d\) adds
the radial factor \(R^{d-4h}\).  Therefore no lift in this even-height class can
produce an \(R\)-primitive witness in degrees \(1,2,\) or \(3\) modulo
four.  Odd height dependence does reach the even class omitted by this
restriction, as the following construction shows.

### Odd height factors and all even primitive degrees

**Status.** This extension is an internal deduction, checked by a focused
AI proof audit and exact bounded computations. It has no independent
mathematical review, formal verification or novelty claim. The all-order
argument below supplies the deduction; finite checks do not prove it.

For each integer \(k\geq0\), use \(D\) from (5.6) and put
\[
 \boxed{F_k=T^k(R+Z)^{k+1}D,\qquad d=2k+4.}
 \tag{5.11}
\]
Here \(d\) is balanced degree: \(F_k\) has bidegree \((d,d)\) and
ordinary total degree \(2d\). The notation \(F_k\) in this subsection
denotes the height family (5.11), distinct from the radially padded
\(F_d\) in (5.2). In particular \(F_0=2F\).

> **Theorem 5.5 (internal deduction: all even primitive degrees).**
> Over every characteristic-zero field, for all integers \(k\geq0\)
> and \(m\geq1\),
> \[
>  \mathcal E_2(F_k^m)=0,
>  \qquad
>  \mathcal E_2(ZF_k^m)=
>  \begin{cases}
>   (dm+2)!\displaystyle\int_0^1w^{km}(1-w^2)^m\,dw,&km\text{ even},\\
>   0,&km\text{ odd}.
>  \end{cases}
>  \tag{5.12}
> \]
> In the even case the displayed rational scalar is strictly positive
> over \(\mathbb R\), hence nonzero in every characteristic-zero field.
> Moreover \(R\nmid F_k\), and \(F_k\) is not a proper power, even up
> to a nonzero scalar. Thus every even balanced degree \(d\geq4\)
> admits an \(R\)-primitive, non-proper-power two-pair Image witness,
> with the fixed marker \(Z\).

**Proof of the moment identities.** On the Hopf sphere \(R=1\),
\(Z=x\), \(W=y\), \(T=t\), and \(t^2+2xy=1\), so
\[
 D=2y-(2+x)t^2=x^{-1}\bigl(1-t^2(1+x)^2\bigr).
\]
Consequently the angular polynomial is
\[
 p_k=x^{-1}t^k(1+x)^{k+1}
       \bigl(1-t^2(1+x)^2\bigr).
 \tag{5.13}
\]
The Laurent expression represents a polynomial on the sphere because it
was obtained from (5.11). Phase averaging takes its constant term in
\(x\). If \(km\) is odd, both constant terms, for \(p_k^m\) and
\(xp_k^m\), are odd polynomials in the uniform height \(t\); their
averages vanish.

Suppose \(km\) is even. Symmetry changes the height average on
\([-1,1]\) into integration on \([0,1]\). Set
\[
\begin{aligned}
 H_{k,m}(X)
 &=X^{(k+1)m}\int_0^1t^{km}(1-t^2X^2)^m\,dt\\
 &=X^{m-1}J_{k,m}(X),\qquad
 J_{k,m}(X)=\int_0^Xw^{km}(1-w^2)^m\,dw.
\end{aligned}
 \tag{5.14}
\]
These are polynomial identities, including at \(X=0\). Constant-term
extraction gives
\[
 \mathbb E_U(p_k^m)=[x^m]H_{k,m}(1+x),\qquad
 \mathbb E_U(xp_k^m)=[x^{m-1}]H_{k,m}(1+x).
 \tag{5.15}
\]
Now \(J_{k,m}'(X)=X^{km}(1-X^2)^m\) vanishes to order \(m\) at
\(X=1\). Hence \(J_{k,m}(X)-J_{k,m}(1)\) is divisible by
\((X-1)^{m+1}\), and
\[
 H_{k,m}(X)\equiv J_{k,m}(1)X^{m-1}
                       \pmod{(X-1)^{m+1}}.
 \tag{5.16}
\]
The coefficients in (5.15) are therefore \(0\) and \(J_{k,m}(1)\),
respectively. The radial degrees are \(dm\) and \(dm+1\), so (2.5)
supplies the factors \((dm+1)!\) and \((dm+2)!\). This proves (5.12).

The integral is positive because its integrand is positive on \((0,1)\).
It also has the exact rational forms
\[
 J_{k,m}(1)=\sum_{j=0}^m
             \frac{(-1)^j\binom mj}{km+2j+1}
          =\frac{2^m m!}{\prod_{j=0}^m(km+2j+1)}.
 \tag{5.17}
\]
For the product formula, integration of the derivative of
\(w^{a+1}(1-w^2)^m\), with \(a=km\), gives
\((a+2m+1)I_{a,m}=2mI_{a,m-1}\) for
\(I_{a,m}=\int_0^1w^a(1-w^2)^m\,dw\). Iterate with \(a\) held
fixed, starting from \(I_{a,0}=1/(a+1)\), then set \(a=km\).

There is no analytic assumption in the characteristic-zero extension:
one may instead use the monomial coefficient formula (4.2) with
\(R=1\), \(T=2t-1\), and \(x=(1-t)/u\). Its constant term is
unchanged by this substitution over \(\mathbb Q(t)\). The resulting
height average is exactly the one above, and polynomial integration is
the rational functional used in Section 4. All identities therefore
hold over \(\mathbb Q\) and extend to every characteristic-zero field.

**Proof of primitivity and the power assertion.** At
\((\xi_1,\xi_2,z_1,z_2)=(1,1,1,-1)\),
\[
 R=0,\quad T=2,\quad R+Z=-1,\quad D=4,\qquad
 F_k=(-1)^{k+1}2^{k+2}\ne0.
 \tag{5.18}
\]
Thus \(R\nmid F_k\). Put \(L=R+Z\). The bilinear forms \(T\)
and \(L\) have coefficient matrices
\[
 \begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
 \begin{pmatrix}1&1\\0&1\end{pmatrix}.
\]
Both have rank two, so each is irreducible over every field extension:
a nontrivial factorization of a form of bidegree \((1,1)\) would be a
product of one linear form in \(\xi\) and one in \(z\), giving rank
one. They are distinct irreducibles. At \((1,1,1,1)\) one has
\(T=0,D=16\); at \((1,1,2,-1)\) one has \(L=0,D=-1\).
Thus neither divides \(D\), and their exact multiplicities in \(F_k\)
are \(k\) and \(k+1\). A proper-power exponent greater than one would
divide both consecutive integers, which is impossible. For \(k=0\),
the multiplicity one of \(L\) already suffices.

For even \(k\), the marker in (5.12) survives at every \(m\geq1\).
For odd \(k\), it survives exactly at even \(m\). Either case gives
infinitely many nonzero marked contractions and therefore violates the
eventual Image/Mathieu conclusion. The argument needs infinite survival;
survival at every power is unnecessary.

The earlier quartic identities are exact specializations: for \(k=0\),
(5.17) is \(J_{0,m}(1)=2^m m!/(2m+1)!!\) and \(F_0=2F\).
Dividing the contractions in (5.12) by \(2^m\) recovers Theorem 1.1.
Combining this specialization with Lemma 5.1 recovers (5.3) in every
balanced degree \(d\geq4\). These implications use the stated
normalizations; no cubic positivity theorem is needed for the quartic
counterexample or its padding.

### Invariant padding and the remaining primitive frontier

There is a sharp obstruction internal to the invariant-multiplier
strategy.  The first fundamental theorem for one vector and one covector
gives
\[
 k[\xi _1,\xi _2,z_1,z_2]^{\mathrm{SL}_2}=k[R].
 \tag{5.10}
\]
Consequently every homogeneous invariant multiplier of balanced degree
\((k,k)\) is a scalar multiple of \(R^k\).  No alternative invariant
multiplier can remove the residual factors \(R,R^2,R^3\) in the three
nonzero congruence classes.  Equivalently, invariant-power propagation is
zero in the primitive quotient \(V_d/RV_{d-1}\) whenever it uses radial
padding.

This also records exactly what the two successful operations contribute.
Ordinary powers, the simplest plethystic operation here, produce the
\(R\)-primitive congruence class \(d\equiv0\pmod4\).  Tensoring with the
invariant line spanned by \(R^k\) fills the other degrees.  General Cartan
projections, plethysms other than ordinary powers, polarizations, block
embeddings, and differential intertwiners have no automatic scalar
radial-shift identity, so no moment-preservation claim is made for them.
Theorem 5.5 supplies primitive witnesses in every even degree \(d\geq4\).
The remaining primitive degree frontier is odd \(d\geq5\), beginning
with degree five. Require \(F\bmod R\ne0\) from the start when seeking
such a witness. Mixed-height or multiple-profile constructions are
possible research directions; no cancellation theorem for them is claimed.
Invariant padding cannot solve this primitive question. A Cartan or
trace-free projection would have to solve new moment identities rather
than inherit them from the radial-shift lemma.

### The affine-multiplier quintic family is excluded

The first quintic deformation of the retained quartic is already rigid.
In this subsection \(a,b,c,e\) are arbitrary coefficients in a
characteristic-zero field; the letter \(e\) is used for the coefficient
of \(W\) to keep it distinct from balanced degree. Put
\[
 G_{a,b,c,e}=(aR+bT+cZ+eW)(R+Z)D.
 \tag{5.19}
\]
This is the complete four-parameter family obtained by multiplying the
retained quartic by a bilinear form: \(R,T,Z,W\) span all four
bilinears. Its elements have balanced degree five.

> **Theorem 5.6 (internal deduction: affine quintic obstruction).**
> If \(\mathcal E_2(G_{a,b,c,e}^{\,j})=0\) for \(j=1,2,3,4\),
> then \(b=c=e=0\). Conversely, for every \(a\), the radial multiple
> \(G_{a,0,0,0}=aR(R+Z)D\) has every pure contraction zero.
> Consequently this entire affine-multiplier family contains no
> \(R\)-primitive quintic Image witness. The unrestricted primitive
> degree-five problem remains open.

**Proof.** On the sphere write
\[
 p_0=(1+x)\bigl(2y-(2+x)t^2\bigr),\qquad
 \ell=a+bt+cx+ey,
\]
and let \(M_j\) be the normalized sphere moment of \((\ell p_0)^j\).
The contraction formula (4.2) gives
\(\mathcal E_2(G^{\,j})=(5j+1)!M_j\). All calculations below are
identities over \(\mathbb Q[a,b,c,e]\).

For clarity, the exact monomial rule used to obtain the moments is
\[
 L(x^i y^j t^h)=0\quad(i\ne j\text{ or }h\text{ odd}),\qquad
 L(x^i y^i t^{2s})
   =\frac1{2^i}\sum_{v=0}^i
       \frac{(-1)^v\binom iv}{2s+2v+1}.
 \tag{5.20}
\]
It follows by substituting \(xy=(1-t^2)/2\), extracting the phase
constant term and integrating the even height. The first two moments are
\[
 M_1=\frac{10c-3e}{15},\qquad
 M_2=\frac{2(56ac+4ae+4b^2+28c^2-4ce+3e^2)}{105}.
 \tag{5.21}
\]
Thus \(M_1=M_2=0\) forces
\[
 c=\frac{3e}{10},\qquad
 b^2=-\frac{130ae+27e^2}{25}.
 \tag{5.22}
\]
Odd powers of \(b\) vanish in each \(M_j\) by height parity. After
these substitutions, the next two moments reduce to
\[
\begin{aligned}
 M_3&=\frac{8e(32175a^2+56680ae+14254e^2)}{625625},\\
 M_4&=\frac{256e(60775a^3+79220a^2e+1802ae^2-643e^3)}{31906875}.
\end{aligned}
 \tag{5.23}
\]
If \(e\ne0\), the ratio \(u=a/e\) would be a common root of
\[
 P(u)=32175u^2+56680u+14254,\qquad
 Q(u)=60775u^3+79220u^2+1802u-643.
\]
The following integer identity rules this out over every
characteristic-zero field:
\[
\begin{aligned}
 A(u)P(u)+B(u)Q(u)&=3472141285450021,\\
 A(u)&=-1531111624900u^2-1247587548025u+290136844604,\\
 B(u)&=810588507300u+1031834055265.
\end{aligned}
 \tag{5.24}
\]
Therefore \(e=0\), and (5.22) gives \(c=0\) and \(b=0\). For the
converse, (5.19) becomes a radial multiple of \(F_0=2F\), so Theorem
5.5 and the radial-shift lemma give all pure contractions zero. Such a
multiple is zero or divisible by \(R\), proving the stated obstruction.

The [quintic verifier](../scripts/verify_two_pair_image_quintic_affine_multiplier.py)
and [exact receipt](two_pair_image_quintic_affine_multiplier_verification.json)
reconstruct (5.21)--(5.24), both from sphere moments and from the
original four-variable contractions. Four moments suffice for this
restricted family because the printed elimination proves the conclusion;
this is not a finite-prefix claim for arbitrary quintics. The proof and
computations are internal, with independent mathematical review pending.

A broader obstruction is proved in the
[real mixed-height profile classification](HOPF_LIFT_CLASSIFICATION.md#7-real-profiles-with-mixed-height-powers).
For \(p=\gamma x^{-r}(1+x)^qP(t(1+x))\), with real coefficients,
\(r,q\geq1\), and \((1-w^2)^r\mid P(w)\), all pure moments vanish
exactly when \(q\leq r\). The case \(q<r\) satisfies eventual mixed
vanishing for every fixed original multiplier; the failing case
\(q=r\) has least balanced degree \(2\deg P\), which is even.
Thus arbitrary mixed height powers inside this single real phase profile
do not produce a primitive odd-degree witness. Theorem 5.6 separately
excludes all complex affine multiples of the retained quartic. Multiple
phase profiles and the rest of the full quintic coefficient space remain
outside these obstructions.

## 6. Frontier consequences

The family begins in balanced bidegree \((4,4)\). The complete
[cubic positivity theorem](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md) proves
the Image implication for the full bidegree-\((3,3)\) stratum, and for
each lower balanced bidegree. Thus four is the exact least failing
balanced degree among bihomogeneous forms in two pairs.
Since \(F_0\) is \(R\)-primitive and not a proper power, four is also
the exact least failing primitive balanced degree. The lower bound uses
the cubic and lower-degree positivity theorem; the even-degree
construction alone supplies the upper bound.

The family also falsifies the proposed two-pair moment--nullcone equality
in every degree \(d\geq4\): the pure moments of \(F_d\) all vanish, while
the fixed multiplier \(Q\) proves that \(F_d\) cannot be one-sided.
The cubic source separately gives a safe, moment-zero cubic outside the
orthogonal nullcone. Consequently moment--nullcone failure in degree
three does not contradict the cubic Image theorem.

Because the coefficient matrix has full rank, the witness does not
contradict the split-symbol GVC theorem.  That theorem concerns rank-one
forms \(A(\xi)P(z)\), whereas (1.4) is nonseparable.  The two-variable GVC
frontier and the ordinary-Laplacian polarization problem therefore remain
open in exactly their previously stated nonhomogeneous or degree-raising
forms.

## Reproduction

Run

```bash
python3 research/foundations/scripts/verify_two_pair_image_mathieu_counterexample.py
python3 research/foundations/scripts/audit_two_pair_image_mathieu_coefficient_extraction.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_even_primitive.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.SIC2C4FiniteSum
```

The dependency-free checker builds \(F\) from (1.2), verifies (1.3), checks
the sixteen-term expansion and coefficient-matrix determinant, and performs
exact sparse contraction through \(m=8\).  It also checks (5.3) for
\(4\leq d\leq10\) and \(1\leq m\leq4\), checks the bounded-radial-order
family (5.4)--(5.5) for \(4\leq d\leq15\) and \(1\leq m\leq3\), checks
the non-power profiles (5.7)--(5.9) for \(1\leq h\leq4\) and
\(1\leq m\leq3\), and replays through \(m=99\) the two finite sums
obtained by expanding (3.3).  The arguments in Sections 2--5, rather than
any finite cutoff, prove (1.5), (5.3), (5.5), and (5.9) for every allowed
parameter.

The second dependency-free audit is deliberately separate.  It checks the
formal coefficient-extraction substitution, the polynomial-division
residue identity (4.8), the finite-difference cancellations, and the recurrence
(4.11) using exact rational arithmetic.  Its finite cutoff is only a
regression; the displayed degree bounds, divisibility at \(X=-1/2\), and
termwise recurrence are the independently checkable all-order certificate.
The Lean module formalizes the general finite-difference cancellation, the
rank-one residue identity (4.8), the specialized normalized products
(4.7)--(4.9), and the recurrence and explicit double-factorial evaluation
of \(B_m\).  It also formalizes the coefficient functional and beta
recurrence, the linear assembly of (4.2), the literal Laurent witness (4.3),
its all-order binomial expansion, both constant-term extractions, and the
closed values in (4.4)--(4.5).  It now also defines the original
four-variable \(F,Q\) as `MvPolynomial` objects and identifies their
displayed substitution with the Laurent witness, including all pure and
mixed powers.  The only remaining Lean integration seam is a wrapper
identifying the generic balanced coefficient array in the contraction
theorem with `MvPolynomial.coeff`; the algebraic proof itself is formalized
on both sides of that representation boundary.

No term-count minimality or literature-priority claim is made here.

The [even primitive verifier](../scripts/verify_two_pair_image_even_primitive.py)
and its [receipt](two_pair_image_even_primitive_verification.json) record
the focused checks for Theorem 5.5. Original four-variable contractions are
checked for \(k=0,1,2,3\), \(m=1,2,3\), and degree-six sphere
calculations are checked through \(m=8\), including both parity cases.
These are internal bounded checks with no independent-review or
formalization claim. The retained receipt identifies the actual local
checker and run performed while recording the supplied observation.

## Source for the Image framework

A. van den Essen, D. Wright, and W. Zhao,
[*On the Image Conjecture*](https://arxiv.org/abs/1008.3962),
J. Algebra 340 (2011), 211--224, supplies the contraction
image-kernel identity and the one-pair positive theorem used for the sharp
dimension conclusion. Definition 1.1 requires eventual membership for each
fixed multiplier, so infinitely many nonzero marked contractions suffice
to contradict that conclusion.

For the one-pair lower bound, the precise input is Theorem 2.8 of the
[retained primary edition](sources/van_den_essen_wright_zhao_1008_3962v2.pdf),
arXiv:1008.3962v2, printed page 5, with its proof on pages 5--6.
Apply it with \(A=k[\xi]\), \(a=\xi\), and \(B=A[z]\).
The characteristic-zero assumption makes \(A\) a \(\mathbb Q\)-algebra;
\(\xi\) is a non-zero-divisor and \((\xi)\) is radical because
\(A/(\xi)=k\) is a field. Therefore
\((\partial_z-\xi)k[\xi,z]\) is Mathieu. Together with the two-pair
failure established above, this proves the minimum failing pair
dimension is exactly two. The imported theorem uses the one-variable
Factorial theorem, Theorem 4.9 in the same source; its hypotheses are not
an extra assumption on the present construction. This records the exact
source application, without asserting independent review of that source.

The failure also persists in every pair dimension \(r\geq2\): view
\(F\) and \(Z\) as polynomials in the first two pairs of the larger
ring. The monomial contraction rule gives
\(\mathcal E_r(F^m)=\mathcal E_2(F^m)\) and
\(\mathcal E_r(ZF^m)=\mathcal E_2(ZF^m)\), since all exponents in
the additional pairs are zero. The image-kernel identity therefore gives
failure for every \(r\geq2\), alongside the positive one-pair theorem.
