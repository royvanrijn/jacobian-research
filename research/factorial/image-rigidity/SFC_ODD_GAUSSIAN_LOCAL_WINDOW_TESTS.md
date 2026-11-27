# Local tests for the odd-start quadratic Gaussian window

For a positive order \(q\), write
\[
 G_q(A,B,C)=\mathbb E(A+B X_q+C X_q^2)^q,\qquad
                        X_q\sim\mathcal N(0,1/q).
\]
These are rational homogeneous polynomials of degree \(q\).
Nonvanishing of the three consecutive forms on projective
coefficient space gives the
[large-height exclusion and fixed-support decision](SFC_INDEPENDENT_COMMON_HEIGHT_DECISION.md#3-three-bivariate-columns-of-equal-degree)
for homogeneous bivariate trinomials at the same start. The
[even-start theorem](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md)
handles all even starts. Here \(n\) is even and the window
consists of \(n-1,n,n+1\).

## 1. Degree and the odd-pair resultant

**Theorem.** For every positive odd \(j\), the projective curves
\(G_j=0\) and \(G_{j+2}=0\) have no common component. In the
chart \(B=1\), the resultant in \(x=A/B\) of the two forms,
as a polynomial in \(h=C/B\), has a zero of exact order
three at \(h=0\). Every closed point of the odd-pair
intersection has degree at most \(j(j+2)\) over \(\mathbb Q\).

Write \(F_j(x,h)=\mathbb E(X_j+x+hX_j^2)^j\). It is monic in
\(x\). At \(h=0\) its roots are \(\mathrm i\) times the scaled
Hermite nodes. Each is simple, and the two odd orders share
only the zero node by
[scaled Hermite isolation](SFC_SCALED_HERMITE_ROOT_ISOLATION.md).
At that node the derivative in \(x\) is
\(j\mathbb E X_j^{j-1}\ne0\).
The ratio of the derivatives in \(h\) and \(x\) is one,
so the simple root has an expansion
\[
                x_j(h)=-h+\frac{2(j-1)}{3j}h^3+O(h^5).
\]
Indeed the parity identity \(F_j(-x,-h)=-F_j(x,h)\) rules
out even powers. For the cubic coefficient substitute
\(x=-h+\gamma h^3\) and compare the terms of order three
after dividing by \(j\mathbb E X_j^{j-1}\). The cubic term
from \(h(X_j^2-1)\) is
\[
 \frac{\binom j3}{j}
 \frac{\mathbb E X_j^{j-3}(X_j^2-1)^3}{\mathbb E X_j^{j-1}}
                  =-\frac{2(j-1)}{3j}.
\]
For \(j\ge3\), under the tilt by \(X_j^{j-3}\) the consecutive
first three moments of \(X_j^2\) are
\((j-2)/j,(j-2)/j,(j^2-4)/j^2\), giving the displayed value.
For \(j=1\) the cubic term is absent and the root is \(-h\).

The cubic coefficients at \(j\) and \(j+2\) differ by
\(4/(3j(j+2))\). The monic resultant is the product of the
cross-order root differences. All pairs other than the two
zero branches have nonzero differences at \(h=0\).
This proves exact order three and excludes a common component
meeting \(B\ne0\). The line \(B=0\) is not a component because
the first form is nonzero at \([A:B:C]=[1:0:0]\).
The degree bound is plane Bézout, applied to the proper
intersection and then to the full Galois orbit of any point.

## 2. The pure quadratic direction at intermediate primes

**Theorem.** At every prime \(n+1<p<2n-1\), a common zero
of \(G_{n-1},G_n,G_{n+1}\) over \(\overline{\mathbb Q}_p\)
cannot reduce to \([0:0:1]\).

More generally put \(P=(p-1)/2\), and consider at least three
distinct orders in \(P<q<p\). If a common zero reduces to
\([0:0:1]\), all the orders other than the least one must
take the same value under \(q\mapsto 4q^2+q\bmod p\).
Consequently four orders in this interval already rule out
that reduction, and for three orders the two upper orders
\(q_1,q_2\) must satisfy \(4(q_1+q_2)+1=0\bmod p\).

Normalize \(C=1\), put \(T=B/2\), \(R=A-B^2/4\), and use
\(v(p)=1\). Then \(\alpha=v(R)>0,\ \beta=v(T)>0\), with
infinity allowed. Completing the square gives the coefficient
of \(R^aT^{2b}\) at order \(q\) as
\[
 c_{ab,q}=\binom qa\binom{2(q-a)}{2b}
              \frac{(2(q-a-b)-1)!!}{q^{q-a-b}} .
\]
For \(a+b\le q\), its valuation is
\(\mathbf 1_{\{q-a>P\}}-\mathbf 1_{\{b>P\}}\).
Indeed the binomial times the double factorial equals
\((2(q-a))!/((2b)!\,2^{q-a-b}(q-a-b)!)\), and all other
factors are units. The two indicators cannot give a
negative valuation.

Set \(D_q=q-P\). The only candidates for minimum valuation
are \(1,R^{D_q},T^{p+1}\) with valuations
\(1,D_q\alpha,(p+1)\beta\), ignoring zero terms.
The constant term is uniquely smallest among the
coefficients of valuation one. In the valuation zero region
either \(a\ge D_q\) or \(b\ge P+1\), giving the two remaining
candidates. The two constant bounds \(1,(p+1)\beta\) must
agree. Otherwise, a tie with \(D_q\alpha\) is possible
for at most one \(q\), and no tie is possible if \(R=0\).
Let \(q_0\) be the least order. We obtain
\[
 \beta=\frac{1}{p+1},\qquad D_{q_0}\alpha\ge1 .
 \tag{1}
\]
In particular \(\alpha>2\beta\), since \(D_{q_0}<P+1\).
At an order \(q>q_0\), the terms with \(a\ge D_q\)
have valuation at least \(D_q\alpha\ge1+\alpha\)
when \(R\ne0\), so they lie strictly above \(1+2\beta\).
All other terms involving \(R\) also lie above
\(1+2\beta\): their coefficient has valuation one or
they have \(b\ge P+1\).

Define \(V=T^{p+1}/p\), a unit, and put
\(k_q=p c_{0,P+1,q}/c_{00,q}\).
At all \(q>q_0\), the terms of valuation at most \(1+2\beta\)
in \(G_q\), after dividing by \(c_{00,q}\), give
\[
 0=1+k_q V+T^2(a_q+k_q V b_q)+E_q,\qquad v(E_q)>2\beta ,
 \tag{2}
\]
where \(a_q=c_{01,q}/c_{00,q}\) and
\(b_q=c_{0,P+2,q}/c_{0,P+1,q}\). These are successive
coefficient ratios, since
\[
 \frac{c_{0,b+1,q}}{c_{0,b,q}}
         =\frac{2q(q-b)}{(2b+2)(2b+1)}.
\]
Hence \(a_q=q^2\) and
\(b_q\equiv q(2q-1)/6\bmod p\). There are no terms
omitted at the valuation of the correction in (2):
for \(a=0\) the next candidates have \(b\ge2\) in
the valuation one region or \(b\ge P+3\) in the
valuation zero region.

The constant residue in (2) gives \(\overline{k_q V}=-1\)
for all upper orders. Thus any difference of \(k_q\)'s
has valuation at least one, since these are rational
\(p\)-adic units. Subtract (2) at two upper orders, divide
by \(T^2\), and reduce. The coefficient
\(\overline{a_q-b_q}=(4q^2+q)/6\) must be the same for all
upper orders. There are at most two distinct upper orders
with the same quadratic value. For two such orders this is
equivalent to \(4(q_1+q_2)+1=0\bmod p\).

For the consecutive window we can also use the constant
residue in (2) in its explicit form. Wilson's congruence
gives
\[
 \overline{k_q}=-(2q)^{P+1}\,q(q-1)\cdots(q-P).
\]
Equality at \(q=n,n+1\) requires
\[
 1=\left(\frac{n+1}{n}\right)^{P+1}\frac{n+1}{n-P}
   =\chi\left(\frac{n+1}{n}\right)
                  \frac{2(n+1)^2}{n(2n+1)} ,
\]
with \(\chi\) the quadratic character. The two signs imply
\(p\mid (3n+2)(4n^2+5n+2)\). The second residue condition
requires \(8n+5=0\bmod p\). At \(n=-5/8\) the two factors
reduce to \(1/8\) and \(7/16\). The only possible prime is
therefore \(7\), which would require \(n=5\) from the interval
and fails \(8n+5=0\bmod7\). This proves the theorem.

## 3. Exact tests over finite fields

Multiplying \(G_q\) by \(q^q\) gives an integer polynomial,
with the coefficient of \(A^a B^{2b} C^{q-a-2b}\) equal to
\[
 \binom{q}{a,\,2b,\,q-a-2b}\,
      (2(q-a-b)-1)!!\;q^{a+b}.
\]
A common projective zero over \(\overline{\mathbb Q}\) would
reduce to a common projective zero over the algebraic closure
of every residue field. A modular test must therefore keep
the homogeneous directions at infinity.

There are two tests. At a prime \(p>2n+2\), check the chart
\(C=1\) by a Gröbner basis in \(A,B^2\); over an algebraic
closure replacing \(B\) by \(B^2\) loses no coefficient
directions. On \(C=0\), the point \(B=0\) is excluded by the
pure \(A\) term. Test \(C=0,\ B=1\) by a univariate gcd.
Unit ideals in both charts give an exact exclusion.

At an intermediate prime of Section 2, it suffices that the
reductions in \(A=1\) and in \(A=0,\ B=1\) have unit ideals.
The only direction outside those charts is \([0:0:1]\),
which cannot lift by Section 2. This combines exact residue
calculations with a common-lift obstruction. The pure
quadratic direction always survives in the residue system
since \(2(n-1)-1\ge p\).

### Lifting test at a non-pure residue point

In \(A=1\), put \(U=B^2\), \(V=C\), and
\(F_q(U,V)=q^q G_q(1,B,C)\). For a rational residue point,
choose integer representatives \(u,v\). The three rows
\[
 \left(\partial_U F_q(u,v),\ \partial_V F_q(u,v),\
                         F_q(u,v)/p\right)\bmod p
 \tag{3}
\]
give an exclusion if their determinant is nonzero.
The first two columns then contain an invertible two-row minor.
Any lift to the algebraic closure has both coordinate corrections
of valuation at least one. Otherwise their minimum valuation
would lie in \((0,1)\). For the two rows of the invertible minor,
some linear term would have that valuation, while the constant
term has valuation at least one and all nonlinear terms have
strictly higher valuation. Now divide the coordinate corrections
by \(p\), apply the three first-order equations and reduce. The
result is a linear system whose augmented matrix is (3), so it
has no solution. The same argument applies at a point over a
finite extension of the residue field using representatives
in its unramified extension of \(\mathbb Q_p\).

For example, at start \(71\) and prime \(109\), the
\(A=1\) residue ideal has reduced basis \(V-43,U+17\).
Thus the non-pure point \((U,V)=(92,43)\) survives the
modular test, and the \(A=0,B=1\) chart has unit ideal.
The rows in (3) are \((67,68,105),(34,32,90),(98,39,56)\),
with determinant \(106\bmod109\). This supplies a
certificate at \(109\) by combining Section 2 with the
first-order lift test at the other point.

### The constant-zero chart at fixed distance from the prime

For \(q=p-d\) in the intermediate interval calculation, put
\(D=q-(p-1)/2\). The restriction of \(q^qG_q\) modulo \(p\)
to \(A=0,\ B=1\) is a nonzero scalar times
\[
 J_d(C)=\sum_{t=0}^{\lfloor(d-1)/2\rfloor}
  \frac{(d-1)_{\underline{2t}}}{(2d-2)_{\underline{2t}}}
         \frac{(d/2)^t}{t!}\,C^{d-1-2t},
 \tag{2}
\]
where the underlined exponents are falling factorials and
the displayed formula is valid when \(p>2d\). To see this,
the surviving indices \(b\) of the \(B^{2b}\) terms run from
\(D\) to \(\lfloor q/2\rfloor\), because the double factorial
index is \(q-b\). Set \(b=D+t\) and divide by the coefficient
at \(t=0\). The multinomial factor gives
\((d-1)_{\underline{2t}}/(2d-2)_{\underline{2t}}\) after
reducing \(2D\equiv 1-2d\bmod p\). The other factors give
\[
 q^t\,\frac{(2(P-t)-1)!!}{(2P-1)!!}
             \equiv \frac{(d/2)^t}{t!}\pmod p.
\]
This proves (2).

In particular, take \(n=p-3\) for a prime \(p\ge11\).
The two upper orders have \(d=3,2\), so the polynomials are
\(C^2+1/4\) and \(C\), with no common root.
The residue direction \([0:0:1]\) cannot lift by
Section 2, including at \(p=23\). Thus for these starts any
algebraic common zero must have \(A\) a unit at every place
over \(p\) in projective integral coordinates. To prove
exclusion at such a start by this route, a unit-ideal test
in \(A=1\) is enough.

## 4. A verified prefix of odd windows

**Theorem.** The quadratic Gaussian window has no common
nonzero projective coefficient direction for every odd
start from 1 through 155 over any characteristic-zero field.

The [exact checker](verify_sfc_odd_gaussian_local_windows.py)
uses the two large-prime charts of Section 3 at \(65537\)
for all 64 odd starts through 127. The remaining 14 starts
are checked at the intermediate prime \(163\), where the
new local comparison applies even at start \(155\), when
the constant residues alone agree. At orders 156 and 157
both residues of \(k_q\) are 92, and the two residues of
\(a_q-b_q\) are 113 and 23. The commands, run from
the repository root with Python 3 and Singular, are

```sh
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_local_windows.py
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_local_windows.py \
    --first 129 --last 155 --prime 163 --intermediate
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_local_windows.py \
    --first 71 --last 71 --prime 109 --intermediate --residue-point 92 43
```

Combining the odd prefix with the even-start gate proves
the large-height exclusion and the finite independent-height
decision for every fixed homogeneous bivariate trinomial
support and every prescribed start through 156. The cutoff
can depend on the prescribed start. The unresolved modular
condition is to find, for every remaining odd start, an
intermediate prime where all non-pure residues are absent
or have an obstruction to lifting.

An independent sufficient test uses the odd-pair resultant
of Section 1. The
[dyadic reduction](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md#odd-starts-and-the-remaining-valuation-pair)
forces \(v_2(C/B)=v_2(n)/2\) at every genuine common zero
of the triple. Hence it suffices to exclude a Newton polygon
slope \(-v_2(n)/2\) in the resultant in \(h=C/B\).
If that slope occurs, the resultant supplies an odd-pair zero
at the specified valuation. The even-order equation must
still vanish at the same coefficient direction.

## 5. Odd starts with a large prime-power divisor in the middle order

**Theorem.** Let the middle order be \(n=2k p^e\),
with \(k,e\ge1\) and \(p\) an odd prime satisfying
\(p>4k\). Define
\[
 P_k(Z)=\sum_{j=0}^{k}\frac{k_{\underline j}}{(2j)!}Z^j,
 \qquad
 Q_k(Z)=\sum_{j=0}^{k-1}
             \frac{(k-1)_{\underline j}}{(2j+1)!}(-Z)^j,
 \qquad
 R_k=\operatorname{Res}_Z((2k)!\,P_k,(2k-1)!\,Q_k).
\]
The integer \(R_k\) is nonzero. If \(p\nmid R_k\),
the moments \(G_{n-1},G_n,G_{n+1}\) have no common
projective zero over any characteristic-zero field.
In particular, for each fixed \(k\), the theorem covers
all exponents \(e\) at all sufficiently large primes.

The proof uses \(p\)-adic reduction of the middle
equation to locate the only surviving chart and to
determine one valuation there. The odd equations on
that chart either separate at the cubic term or reduce
to \(P_k,Q_k\). The degree-\(k\) test is independent of
\(e\).

### The middle equation forces a valuation

Assume a common projective zero. By the Nullstellensatz,
we may take algebraic coefficients and embed them in
\(\overline{\mathbb Q}_p\). Write \(v=v_p\) and choose
projective integral coordinates with at least one unit.
Use
\[
 V(r)=v((2r-1)!!),\quad V(0)=0,\qquad
 L=\sum_{d=1}^e p^{-d}=\frac{1-p^{-e}}{p-1}.
\]
Since \(2n<p^{e+1}\), no odd factor through \(2n-1\)
has valuation greater than \(e\). In the order-\(n\)
integer form \(n^nG_n\) the coefficient indexed by
\(i\) constant factors and \(b\) linear factors, with
\(b\) even, has valuation
\[
 e(i+b/2)+V(n-i-b/2)+v\binom n{i,b,n-i-b}.
\]
For \(l=i+b/2>0\), the \(l\) odd factors removed from
\((2n-1)!!\) have total valuation strictly less than
\(el\) because the first removed factor \(2n-1\) is a
unit. Thus \(p^{-V(n)} n^n G_n\) reduces to a nonzero
multiple of \(C^n\). Consequently \(\bar C=0\).

In \(G_{n+1}\), the order is a unit congruent to \(1\),
and all Gaussian moments of degrees at least \(p\)
vanish modulo \(p\). In the truncated polynomial
ring modulo \(X^p\) in characteristic \(p\),
\((A+BX+CX^2)^{n+1}=A^n(A+BX+CX^2)\).
The residue of the equation is \(A^n(A+C)=0\).
It forces \(\bar A=0\). Hence \(B\) is a unit.
Set \(x=A/B,\ h=C/B\), so both have positive valuation
with infinity allowed at zero.

Rescale the middle equation to
\[
 \mathbb E(u+Y+wY^2)^n=0,\qquad
 u=\sqrt n\, x,\quad w=h/\sqrt n,\quad Y\sim\mathcal N(0,1).
\]
A nonzero term with \(i,b,j\) constant, linear, and
quadratic factors has \(b\) even, \(i+b+j=n\), and a
Gaussian factor of valuation \(V(r)\) with
\(r=(n-i+j)/2\). At \(v(w)=-L/2\), its valuation is
\[
 \frac{nL}{2}+i\bigl(v(u)-L/2\bigr)
       +V(r)-rL+v\binom n{i,b,j}.
\]
For \(d\le e\), the contribution to \(V(r)\) is
\(\lfloor(r+(p^d-1)/2)/p^d\rfloor\). It equals \(r/p^d\)
if all three counts are divisible by \(p^d\), and it
exceeds \(r/p^d-1/2\) otherwise. In the latter case
the contribution to the multinomial valuation at that
level is at least one because the counts sum to a
multiple of \(p^d\). Thus the last three terms beyond
\(nL/2\) are nonnegative in total. Equality requires
\(i=0\) and all counts divisible by \(p^e\), since
\(v(u)>e/2>L/2\). The pure linear and pure quadratic
terms both achieve equality. If \(v(w)>-L/2\), the
pure linear term is uniquely least, including at \(w=0\).
If \(v(w)<-L/2\), the pure quadratic term is uniquely
least. The middle equation therefore forces
\[
                            c:=v(h)=e/2-L/2 .
 \tag{4}
\]

### The odd equations when \(v(x)\ge v(h)\)

Use the two odd orders \(m=n-1,n+1\). For either of
them divide the chart polynomial by
\(m\mu_{m-1}\), where \(\mu_j=\mathbb E X_m^j\).
The term of degree one in \(x,h\) is \(x+h\).
Only odd total degrees \(\ell\) occur. The valuation
of any coefficient at degree \(\ell\) is bounded below
by minus the valuation of
\(\prod_{s=1}^{(\ell-1)/2}(m-2s)\). Indeed this is the
denominator of \(\mu_{m-\ell}/\mu_{m-1}\) up to units,
and additional factors from powers of \(h\) have
nonnegative valuation.

There is no loss for \(\ell<p\). For \(\ell\ge p\),
the loss is at most \(\ell/(p-1)\): at each level
\(p^d\le\ell\) at most
\(\lfloor(\ell+p^d)/(2p^d)\rfloor\le\ell/p^d\) factors
are divisible by \(p^d\), since \(m\equiv\pm1\bmod p^d\)
for \(1\le d\le e\). Higher levels cannot occur.
For \(p\ge7\) we have \(c\ge(p-2)/(2(p-1))\), so
\((\ell-3)c>\ell/(p-1)\) for \(\ell\ge p\).
For \(p=5\) the same inequality holds when \(e\ge2\).
The sole remaining case \(p=5,e=1\) here has \(k=1\),
hence \(n=10\), already covered by Section 4.

Suppose \(v(x)\ge c\). All degrees at least five
then have valuation strictly larger than \(3c\).
The degree-three coefficients are integral at \(p\).
A zero thus requires \(v(x+h)\ge3c\). Write
\(x=-h+y\), so terms in the cubic involving \(y\)
also have valuation larger than \(3c\).
The cubic at \(x=-h\), divided by \(h^3\), is
\(-2(m-1)/(3m)\), by the tilt computation in Section 1.
The two equations therefore give
\[
  0=y-\frac{2(m-1)}{3m}h^3 + E_m,
             \qquad v(E_m)>3c,\qquad m=n-1,n+1.
\]
The two cubic coefficients differ by a \(p\)-adic
unit. Their equations cannot both hold.

### The odd equations when \(v(x)<v(h)\)

Suppose \(v(x)<c\). In each total degree \(\ell\)
the pure \(x^\ell\) term has strictly smaller
valuation than every term involving \(h\): replacing
an \(x\) by an \(h\) increases the parameter valuation,
and the corresponding Gaussian moment ratio is
\(p\)-integral. Thus a least-valuation comparison in
the two odd equations can use the pure \(x\)-polynomials
alone.

Let \(P=p^e,\ N=n/2=kP,\ \lambda=P L\). After division
by \(m\mu_{m-1}x\), the coefficient of \(x^{2r}\)
in the pure polynomial is
\[
 U^+_r=\frac{(2(n+1))^r N_{\underline r}}{(2r+1)!}
   \quad (0\le r\le N),\qquad
 U^-_r=\frac{(2(n-1))^r (N-1)_{\underline r}}{(2r+1)!}
   \quad (0\le r<N).
\]
This follows by cancelling the factors
\(m-2,m-4,\ldots,m-2r\) against the binomial
coefficient. Since \(k<p/4\),
\(v(N-i)=v(i)\) for \(0<i<N\). Therefore, with
\(W(r)=v((2r+1)!)-v(r!)\), we have
\(v(U^-_r)=-W(r)\) and
\(v(U^+_r)=-W(r)+e-v(r)\) for \(r>0\).
Here \(U^+_0=1\).

If \(s_d=r\bmod p^d\), the \(d\)-th contribution to
\(W(r)\) is
\(\lfloor r/p^d\rfloor+\mathbf 1_{s_d\ge(p^d-1)/2}\).
It follows that \(v(U^+_r)\ge-rL\), with equality
precisely at multiples of \(P\): the correction
\(e-v(r)\) counts the levels where \(s_d\ne0\).
In particular the endpoints are on this line.
The upper odd equation forces \(v(x)=L/2\), since
otherwise one endpoint is uniquely least. At that
valuation the plus polynomial retains \(r=tP\),
\(0\le t\le k\). In the minus polynomial the minimum
of \(v(U^-_r)+rL\) occurs exactly when
\(s_d=(p^d-1)/2\) at every level. These indices are
\(r=(P-1)/2+tP\), \(0\le t<k\).

Put \(Z=(-1)^\lambda 2 x^{2P}/p^\lambda\), a unit.
For an integer \(j\) with base-\(p\) digits \(j_d\),
the unit part of \(j!\) reduces to
\((-1)^{v(j!)}\prod_d j_d!\) by grouping the nonmultiples
of \(p\) in blocks and using Wilson's congruence.
For \(r=tP\) the high digits of \(N,N-r,2r+1\)
are respectively \(k,k-t,2t\), with the lower digit
factorials contributing one. The factor
\((2(n+1))^{tP}\) reduces to \(2^t\), so the plus
residue is \(P_k(\bar Z)\). For
\(r=(P-1)/2+tP\) divide by the first retained term.
The lower digits of \(N-1-r,2r+1\) are independent of
\(t\), and their high digits are \(k-1-t,2t+1\).
The additional factor \((2(n-1))^{tP}\) reduces to
\((-2)^t\). Thus the minus residue is \(Q_k(\bar Z)\).
The assumption on \(R_k\) excludes a common residue.

For completeness, \(P_k(-2X^2)\) is proportional to
the probabilists' Hermite polynomial of degree \(2k\),
and \(Q_k(2X^2)\) to the degree-\((2k-1)\) polynomial
divided by \(X\). This follows from the coefficient
formula obtained from \(\exp(Xt-t^2/2)\).
Hermite orthogonality for the positive Gaussian
weight gives simple real roots. Consequently all
roots of \(P_k\) are strictly negative and all roots
of \(Q_k\) are strictly positive, so \(R_k\ne0\).

For \(k=1\) the second polynomial is constant and
\(R_1=1\). For \(k=2\) the only root of \(Q_2\) is \(6\),
and \(P_2(6)=10\). For \(k=3\) the integer polynomials
may be taken to be \(Z^3+30Z^2+180Z+120\) and
\(Z^2-20Z+60\). Dividing the former by the latter
leaves \(160(7Z-18)\). Substitution of \(18/7\) in
the latter gives \(24\cdot31/49\). Thus the theorem
covers \(n=2p^e\) and \(4p^e\) under the displayed prime
bound, and \(6p^e\) for \(p>12\) with \(p\ne31\).

## 6. Exact tests at prime divisors of the prescribed orders

For fixed even \(n\) and \(p\ge5\) dividing \(n\),
write \(e=v_p(n)\) and \(V(r)=v_p((2r-1)!!)\),
with \(V(0)=0\). Compute the rational bounds
\[
 b=\min_{1\le l\le n}\frac{el+V(n-l)-V(n)}{2l},
 \qquad
 b_*=\max_{\substack{m=n-1,n+1\\ 5\le j\le m,\ j\ {\rm odd}}}
       \frac{V((m-1)/2)-V((m-j)/2)}{j-3}.
\]
If \(b>\max(0,b_*)\), every common direction would
reduce to the pure linear direction. Indeed the strict
positivity of the numerators defining \(b\) gives
the pure quadratic reduction of the middle form
after dividing by its minimum coefficient valuation,
as in Section 5. The reduction of \(G_{n+1}\) then
forces \(\bar A=0\). In \(B=1\) put \(x=A/B\), \(h=C/B\).
Relative to the pure quadratic term, a term of
\(n^nG_n\) with \(i\) constant factors and \(j\) linear
factors has valuation at least
\[
 el+V(n-l)-V(n)-(i+j)v(h),
 \qquad l=i+j/2 .
\]
If \(0<v(h)<b\), this is strictly positive for
\(l>0\), since \(i+j\le2l\). Thus \(v(h)\ge b\)
with infinity allowed. If \(v(x)\ge v(h)\), the bound
\(b_*\) makes all terms above the cubic have strictly
greater valuation than \(3v(h)\). Cubic comparison as
in Section 5 gives a contradiction. If \(h=0\) and
\(v(x)\ge v(h)\), the middle pure linear term already
excludes that direction.

It remains to test \(v(x)<v(h)\). Within each total
degree in the odd equations, retaining only the pure
\(x\) term gives the coefficient of least valuation.
Compute the two Newton polygons from the pure
coefficients \(U^\pm_r\) of Section 5, using their
factorial formulas with \(N=n/2\), without the
restriction \(n=2kp^e\) with \(p>4k\).
A common zero would require a common positive root
valuation of the pure polynomials in \(x\). For each
such valuation, scale \(x\) by the corresponding
power of a uniformizer in a ramified extension and
keep the initial terms of both polynomials. Remove
their initial powers of \(x\). The resulting reductions,
viewed as polynomials in the square of the scaled
variable, must then have a common nonzero root.
The [prime-divisor checker](verify_sfc_odd_gaussian_prime_divisors.py)
uses these exact comparisons and gcds over \(\mathbb F_p\).
The unit parts of factorials needed for the reductions
are computed by stripping factors of \(p\) from each
integer in the product.

### The first and last orders

The same quadratic lower bound applies at a prime
\(p\ge3\) dividing either outer order \(m\), with
\(n\) replaced by \(m\) in the definition of \(b\).
If \(b>0\) the residue of \(C\) is zero. For a prime
dividing the first order, the next order has residue
\(A^m(A+C)\), yielding the pure linear chart.
The lower bound \(v(h)\ge b\) then applies.
If \(b\) exceeds every positive root valuation of
the pure even Hermite polynomial, the even equation
requires \(v(x)<v(h)\). A comparison of the two
remaining pure Hermite polynomials and their initial
forms tests this case exactly.

For a prime dividing the last order, the two
preceding orders are \(p\)-integral. On \(C=0,B=1\)
their reductions are just the scaled linear
Gaussian moment polynomials in \(A\). The chart
\(C=B=0\) is excluded by the pure \(A\) term. If
the reductions have no common nonzero root, then
any surviving residue must have \(A=0\).
The same lower bound and even/odd initial-form
comparison apply in that chart. These endpoint
tests use only factorial valuations and univariate
gcds in addition to the middle-order test.

### Exact coverage of odd starts

Together with Section 4, the three prime-divisor
tests certify every odd start through \(12001\).
The combined test certifies all 5923 odd starts
from \(157\) through \(12001\). For comparison, its
middle-order branch alone certifies 382 of the 423
odd starts from \(157\) through \(1001\).

~~~sh
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_prime_divisors.py 157 5001
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_prime_divisors.py 5003 12001
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_prime_divisors.py \
    157 1001 --middle-only --expect-unresolved \
    159 161 179 185 191 215 223 239 251 255 269 287 319 323 349 359 383 \
    419 431 447 479 485 499 503 511 539 569 575 599 629 639 647 719 \
    767 791 839 863 899 959 971
~~~

Internal replay used Python 3.12.3. The underlying
exclusions hold over characteristic zero. Combined
with the even-start theorem, this gives the
large-height exclusion for fixed homogeneous
bivariate trinomial supports at every start through
\(12002\), as well as the families of Sections 5 and 7.
For every start in this proved domain, the
[explicit cutoff transfer](SFC_EXPLICIT_GAUSSIAN_NONLIFTING_CUTOFF.md)
now supplies a numerical corner bound from the exponent span and
maximum prescribed order. This does not enlarge the all-order scope
of the Gaussian exclusion.

## 7. Large prime powers dividing the initial order

**Theorem.** Suppose the odd start is \(q=(2k+1)p^e\),
where \(k\ge0\), \(e\ge1\), \(p\) is odd prime and
\(p>4k+2\). Define
\[
 S_k(Z)=\sum_{j=0}^k\frac{k_{\underline j}}{(2j+1)!}Z^j,
 \qquad
 T_k=\operatorname{Res}_Z((2k+1)!\,S_k(Z),(2k+1)!\,S_k(2Z)),
 \qquad T_0=1 .
\]
If \(p\nmid T_k\) and either \(e=1\) or \(2\) is a
quadratic nonresidue modulo \(p\), the Gaussian gate
closes over characteristic zero. In particular this
covers all odd prime starts and all powers of every odd
prime at which \(2\) is a nonresidue.

Let \(P=p^e\) and \(L=\sum_{d=1}^e p^{-d}\).
Since \(2q<p^{e+1}\), the reduction argument of
Section 5 applied to \(q^qG_q\) forces \(C\) to
vanish in the residue field. The reduction of
\(G_{q+1}\) is \(A^q(A+C)\), so \(A\) vanishes too
and we may use \(B=1,\ x=A/B,\ h=C/B\).
Both \(x,h\) now have positive valuation.

In the rescaled order-\(q\) equation put
\(u=\sqrt q\,x,\ w=h/\sqrt q\). The term comparison
of Section 5 uses the valuation of the pure
quadratic term at \(v(w)=-L/2\) as a lower bound
for all terms. At smaller \(v(w)\) its quadratic
count is uniquely maximal, so it dominates.
Thus \(v(h)\ge e/2-L/2\), with infinity allowed.

Set \(N=(q+1)/2\). In the even order \(q+1\),
divide the pure polynomial by its constant
Gaussian moment. In the odd order \(q+2\),
divide it by \((q+2)\mu_{q+1}x\). The two
coefficient lists in \(x^2\) are
\[
 E_r=\frac{(2(q+1))^r N_{\underline r}}{(2r)!},
 \qquad
 O_r=\frac{(2(q+2))^r N_{\underline r}}{(2r+1)!},
 \qquad 0\le r\le N .
\]
Write \(M_d=(p^d-1)/2,\ N_d=(p^d+1)/2\), and
\(s_d=r\bmod p^d\). The ordinates of the two
Newton polygons are the negatives of
\[
 v(r!)+\sum_{d=1}^e\mathbf 1_{s_d=N_d},
 \qquad
 v(r!)+\sum_{d=1}^e\mathbf 1_{s_d\in\{M_d,N_d\}},
\]
respectively. To see this, use \(N\bmod p^d=N_d\)
in the Legendre formula for the numerator and
denominator. All factorial arguments are less than
\(p^{e+1}\).

For \(p\ge5\), the even polygon through \(N_e\)
has vertices \(0,N_1,\ldots,N_e\). The odd polygon
through \(M_e\) has vertices \(0,M_1,\ldots,M_e\).
For an edge after the first one ending at \(N_d\)
or \(M_d\), the width is \(D_d=(p-1)p^{d-1}/2\),
and the slope is
\[
             -\sum_{j=1}^{d-1}p^{-j}-1/D_d.
\]
At the left endpoint the offset of each of the
first \(d-1\) level contributions from its linear
function \(r/p^j\) is maximal. These contributions
increase by at most the linear increment at
interior points. At the right endpoint they attain
that increment and the \(d\)-th level adds one.
All interior points therefore lie strictly above
the chord. The first edges have slopes \(-1/N_1\)
and \(-1/M_1\), respectively, with earlier ordinates
zero. The slopes increase strictly. For \(p=3\)
the prime bound forces \(k=0\). When \(e\ge2\) the
even hull starts with the single edge \(0,N_2=5\)
of slope \(-3/5\): the ordinates at \(1,2,3,4\)
are \(0,-1,-1,-1\), above that edge. Subsequent
edges and the odd polygon follow the same argument.
When \(p=3,e=1\), both first edges are as in the
original lists.

After \(N_e\) in the even polygon and \(M_e\) in
the odd polygon there is an edge of slope \(-L\)
and width \(kP\) if \(k>0\). Its initial points
are exactly \(N_e+tP\) or \(M_e+tP\) for
\(0\le t\le k\), since these are the indices
at which all \(e\) level offsets are simultaneously
maximal. The odd polygon then has a final flat
edge of width one.

The positive root valuations of the even pure
polynomial in \(x\) are at most \(3/10\), from
these slopes and the factor two for \(x^2\).
Since \(e/2-L/2>3/10\), the even equation cannot
vanish with \(v(x)\ge v(h)\). All positive
degrees would have strictly greater valuation
than its constant, because coefficients of terms
involving \(h\) are the pure coefficient times
integral Gaussian moment ratios and binomial
coefficients. Thus \(v(x)<v(h)\), and both
equations can be tested by the initial forms of
their pure polynomials, as in Section 6.

For a shared edge of width \(D_d\), the initial
forms are binomials after removing their initial
monomials. At an odd endpoint \(r=M_j\), the even
endpoint is \(r+1=N_j\), and
\[
 \frac{E_{r+1}}{O_r}
 =\frac{(q+1)(N-r)}{r+1}
         \left(\frac{q+1}{q+2}\right)^r
 \equiv 2\cdot2^{-r}\pmod p.
\]
Hence the ratio of the normalized final
coefficients is \(2^{-D_d}\equiv
2^{-(p-1)/2}\pmod p\). A nonresidue gives opposite
coefficients and excludes a common root. Such edges
are absent for \(e=1\). When \(2\) is a residue and
these slopes occur, their initial binomials agree.

For the last shared slope \(-L\), put
\(\lambda=PL\) and \(Z=(-1)^\lambda 2x^{2P}/p^\lambda\).
The unit-factorial computation from Section 5 gives
the normalized even and odd residues \(S_k(\bar Z)\)
and \(S_k(2\bar Z)\). Indeed, the top digits of
\(N-r\) and of the corresponding factorial argument
\(2r\) or \(2r+1\) are \(k-t\) and \(2t+1\).
The lower digits are independent of \(t\). The
extra power \((2m)^{tP}\) gives \(2^t\) at the
even order and \(4^t\) at the odd order.
The assumption \(p\nmid T_k\) excludes this slope.

For \(k>0\), \(S_k(2X^2)\) is proportional to the
odd Hermite polynomial of degree \(2k+1\) divided
by \(X\). Schur's irreducibility of the squared
nonzero Hermite polynomial, used in the
[Hermite isolation proof](SFC_SCALED_HERMITE_ROOT_ISOLATION.md#square-node-groups),
implies that \(S_k\) is irreducible. Its constant
coefficient is one and its linear coefficient
changes under scaling by two. Thus \(S_k(Z)\)
and \(S_k(2Z)\) are coprime over \(\mathbb Q\)
and \(T_k\ne0\). For each fixed \(k\) only
finitely many primes are excluded by this
resultant condition.
